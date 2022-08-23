
# packages ----------------------------------------------------------------
library(splatter)
library(SCPA)
library(testSctpa)
library(tidyverse)
library(magrittr)
library(msigdbr)
library(Seurat)
library(ggnewscale)

# functions ---------------------------------------------------------------
get_pvals <- function(input) {
  
  split_input <- SplitObject(input, split.by = "group")
  
  pval_output <- list()
  for (i in 1:nrow(input)) {
    pval_output[[as.character(i)]] <- wilcox.test(split_input[[1]]@assays$PAS@data[i, ], 
                                                  split_input[[2]]@assays$PAS@data[i, ])$p.value
  }
  
  comp_output <- unlist(pval_output) %>%
    data.frame() %>% 
    set_colnames("pval") %>%
    mutate(pathway = rownames(input))
  
  return(comp_output)
  
}


# import counts -----------------------------------------------------------
load("int_data_N8.RData")


# define a gene set to permute --------------------------------------------
comp_genes <- msigdbr("Homo sapiens", "H") %>%
  filter(gs_name == "HALLMARK_COMPLEMENT") %>%
  pull(gene_symbol) %>%
  unique()

`%notin%` <- Negate(`%in%`)

non_comp_genes <- rownames(int_data_N8)[rownames(int_data_N8) %notin% comp_genes]


# get data to estimate ----------------------------------------------------
counts <- int_data_N8@assays$RNA@counts[1:length(non_comp_genes), 1:1000] %>% as.matrix()
comp_counts <- int_data_N8@assays$RNA@counts[1:length(comp_genes), 1:1000] %>% as.matrix()

# estimate parameters -----------------------------------------------------
params <- splatEstimate(counts)
params_comp <- splatEstimate(comp_counts)

# simulate background data ------------------------------------------------
params <- setParams(params, update = list("group.prob" = c(0.5, 0.5),
                                          "de.prob" = 0.2,
                                          "de.facLoc" = 0.2))

background_expression <- splatSimulate(params, method = "groups")

background_expression <- background_expression@assays@data$counts %>% as.matrix()
rownames(background_expression) <- non_comp_genes

# 1. defacloc data -----------------------------------------------------------
# loop over complement differential expression ----------------------------
de_prob <- seq(0, 1, 0.05)
de_data <- list()
for (i in de_prob) {
  
  cat("Generating data with", i, "de prob", "\n")
  
  variable_params <- setParams(params_comp, update = list("group.prob" = c(0.5, 0.5),
                                                          "de.prob" = 0.3,
                                                          "de.facLoc" = i))
  
  complement_sim <- splatSimulate(variable_params, method = "groups")
  complement_expression <- complement_sim@assays@data$counts
  rownames(complement_expression) <- comp_genes
  
  combined_df <- rbind(complement_expression, background_expression)
  
  combined_df <- CreateSeuratObject(combined_df) %>%
    NormalizeData()
  combined_df$group <- complement_sim$Group
  Idents(combined_df) <- "group"
  
  de_data[[as.character(i)]] <- combined_df
  
}

# test with scpa ----------------------------------------------------------
scpa_out <- list()
for (i in names(de_data)) {
  
  group1 <- seurat_extract(de_data[[i]], meta1 = "group", value_meta1 = "Group1")
  group2 <- seurat_extract(de_data[[i]], meta1 = "group", value_meta1 = "Group2")
  
  print("testing with scpa")
  pathways <- msigdbr("Homo sapiens", "H") %>%
    filter(gs_name == "HALLMARK_COMPLEMENT") %>%
    format_pathways()
  
  scpa_out[[i]] <- compare_pathways(samples = list(group1, group2),
                                    pathways = pathways) %>%
    select(Pval, Pathway) %>%
    set_colnames(c("pval", "pathway")) %>%
    mutate(de = i) %>%
    mutate(method = "scpa") %>%
    mutate(pathway = gsub(pattern = "_", replacement = "-", x = pathway))
  
}

scpa_out <- bind_rows(scpa_out)

# sctpa -------------------------------------------------------------------
### test all methods across all de datasets
sctpa_out <- list()
temp <- list()
method <- c("AUCell", "Vision", "GSVA", "ssGSEA", "zscore")

for (c in method) {
  print(c)
  
  for (i in names(de_data)) {
    
    print(i)
    
    temp[[i]] <- cal_PAS(seurat_object = de_data[[i]],
                         tool = c,
                         normalize = "log",
                         species = "human",
                         pathway = "hallmarker",
                         n_cores = 4) %>%
      get_pvals() %>%
      filter(pathway == "HALLMARK-COMPLEMENT") %>%
      mutate(de = i) %>%
      mutate(method = c)
    
  }
  
  sctpa_out[[as.character(c)]] <- bind_rows(temp)
  
}


de_facloc <- rbind(bind_rows(sctpa_out), scpa_out)

saveRDS(de_data, "defacloc_loop_expression.rds")
saveRDS(de_facloc, "defacloc_results.rds")


# 2. de prob data ------------------------------------------------------------
# loop over complement differential expression ----------------------------
de_prob <- seq(0, 1, 0.05)
de_data <- list()
for (i in de_prob) {
  
  cat("Generating data with", i, "de prob", "\n")
  
  variable_params <- setParams(params_comp, update = list("group.prob" = c(0.5, 0.5),
                                                          "de.prob" = i,
                                                          "de.facLoc" = 0.3))
  
  complement_sim <- splatSimulate(variable_params, method = "groups")
  complement_expression <- complement_sim@assays@data$counts
  rownames(complement_expression) <- comp_genes
  
  combined_df <- rbind(complement_expression, background_expression)
  
  combined_df <- CreateSeuratObject(combined_df) %>%
    NormalizeData()
  combined_df$group <- complement_sim$Group
  Idents(combined_df) <- "group"
  
  de_data[[as.character(i)]] <- combined_df
  
}


# test with scpa ----------------------------------------------------------
scpa_out <- list()
for (i in names(de_data)) {
  
  group1 <- seurat_extract(de_data[[i]], meta1 = "group", value_meta1 = "Group1")
  group2 <- seurat_extract(de_data[[i]], meta1 = "group", value_meta1 = "Group2")
  
  print("testing with scpa")
  pathways <- msigdbr("Homo sapiens", "H") %>%
    filter(gs_name == "HALLMARK_COMPLEMENT") %>%
    format_pathways()
  
  scpa_out[[i]] <- compare_pathways(samples = list(group1, group2),
                                    pathways = pathways) %>%
    select(Pval, Pathway) %>%
    set_colnames(c("pval", "pathway")) %>%
    mutate(de = i) %>%
    mutate(method = "scpa") %>%
    mutate(pathway = gsub(pattern = "_", replacement = "-", x = pathway))
  
}

scpa_out <- bind_rows(scpa_out)


# sctpa -------------------------------------------------------------------
### test all methods across all de datasets
sctpa_out <- list()
temp <- list()
method <- c("AUCell", "Vision", "GSVA", "ssGSEA", "zscore")

for (c in method) {
  print(c)
  
  for (i in names(de_data)) {
    
    print(i)
    
    temp[[i]] <- cal_PAS(seurat_object = de_data[[i]],
                         tool = c,
                         normalize = "log",
                         species = "human",
                         pathway = "hallmarker",
                         n_cores = 4) %>%
      get_pvals() %>%
      filter(pathway == "HALLMARK-COMPLEMENT") %>%
      mutate(de = i) %>%
      mutate(method = c)
    
  }
  
  sctpa_out[[as.character(c)]] <- bind_rows(temp)
  
}


de_prob <- rbind(bind_rows(sctpa_out), scpa_out)

saveRDS(de_data, "deprob_expression.rds")
saveRDS(de_prob, "deprob_results.rds")


# 3. merge all data ------------------------------------------------------
# idea data
idea_deprob <- readRDS("/Volumes/data/benchmarking/idea/sim_data/idea_deprob.rds")
idea_facloc <- readRDS("/Volumes/data/benchmarking/idea/sim_data/idea_facloc.rds")

idea_deprob <- lapply(idea_deprob, function(x) x@gsea %>%
                        filter(annot_id == "HALLMARK_COMPLEMENT")) %>%
  bind_rows() %>%
  select(annot_id, pvalue_louis) %>%
  rename("pathway" = "annot_id", "pval" = "pvalue_louis") %>%
  mutate(method = "iDEA") %>%
  mutate(de = seq(0, 1, 0.05)) %>%
  select(pval, pathway, de, method)

idea_facloc <- lapply(idea_facloc, function(x) x@gsea %>%
                        filter(annot_id == "HALLMARK_COMPLEMENT")) %>%
  bind_rows() %>%
  select(annot_id, pvalue_louis) %>%
  rename("pathway" = "annot_id", "pval" = "pvalue_louis") %>%
  mutate(method = "iDEA") %>%
  mutate(de = seq(0, 1, 0.05)) %>%
  select(pval, pathway, de, method)

# fgseadata
fgsea_deprob <- readRDS("/Volumes/data/benchmarking/fgsea/fgsea_deprob.rds")
fgsea_deprob <- lapply(fgsea_deprob, function(x) x %>%
                         filter(pathway == "HALLMARK_COMPLEMENT")) %>%
  bind_rows() %>%
  select(pathway, pval) %>%
  mutate(method = "fgsea") %>%
  mutate(de = seq(0, 1, 0.05)) %>%
  data.frame() %>%
  select(pval, pathway, de, method)

fgsea_facloc <- readRDS("/Volumes/data/benchmarking/fgsea/fgsea_defacloc.rds")
fgsea_facloc <- lapply(fgsea_facloc, function(x) x %>%
                         filter(pathway == "HALLMARK_COMPLEMENT")) %>%
  bind_rows() %>%
  select(pathway, pval) %>%
  mutate(method = "fgsea") %>%
  mutate(de = seq(0, 1, 0.05)) %>%
  data.frame() %>%
  select(pval, pathway, de, method)


# merge data
facloc <- readRDS("defacloc_results.rds") %>%
  mutate(method = gsub(pattern = "scpa", replacement = "SCPA", x = method)) %>%
  rbind(idea_facloc, fgsea_facloc) %>%
  mutate(method = gsub(pattern = "fgsea", replacement = "fGSEA", x = method))

deprob <- readRDS("deprob_results.rds") %>%
  mutate(method = gsub(pattern = "scpa", replacement = "SCPA", x = method)) %>%
  rbind(idea_deprob, fgsea_deprob) %>%
  mutate(method = gsub(pattern = "fgsea", replacement = "fGSEA", x = method))


# 4.  plot data -----------------------------------------------------------

ggplot(facloc, aes(de, -log10(pval))) +
  geom_smooth(aes(group = method, col = method), se = F, lwd = 0.6) +
  geom_point(shape = 21, size = 2.5, alpha = 0.8, aes(fill = method)) +
  scale_x_discrete(breaks = seq(0, 1, 0.2), expand = c(0.05, 0.05)) +
  facet_wrap(~method) +
  labs(x = "Size of DE factors", title = "Increasing DE size") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "none")
ggsave("facloc_results.png", device = "png", dpi = 600,
       width = 5, height = 4.5, units = "in")

ggplot(deprob, aes(de, -log10(pval))) +
  geom_smooth(aes(group = method, col = method), se = F, lwd = 0.6) +
  geom_point(shape = 21, size = 2.5, alpha = 0.8, aes(fill = method)) +
  scale_x_discrete(breaks = seq(0, 1, 0.2)) +
  facet_wrap(~method) +
  labs(x = "DE probability", title = "Increasing DE probability") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "none")
ggsave("deprob_results.png", device = "png", dpi = 600,
       width = 5, height = 4.5, units = "in")









