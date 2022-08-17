
# packages ----------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(glmGamPoi)
library(SCPA)
library(msigdbr)
library(magrittr)
library(ComplexHeatmap)
library(AnnotationDbi)
library(hgu133plus2.db)
library(janitor)



# 1.  process tissue data -------------------------------------------------
# metadata
meta <- read.csv("Metadata.csv", stringsAsFactors = F)

# processing functions
proc_data <- function(filepath) {
  df <- read.table(filepath, header = T, stringsAsFactors = F) %>%
    mutate(Gene = make.names(Gene, unique = T)) %>%
    column_to_rownames("Gene") %>%
    select(-Accession)
}

seu_obj <- function(mat) {
  df <- CreateSeuratObject(mat, min.cells = 3, min.features = 300) %>%
    PercentageFeatureSet(., pattern = "MT.", col.name = "percent_mt") %>%
    subset(., subset = percent_mt < 10)
}

# load in and process files
seq_files <- grep(list.files(), pattern = "^GSM", value = T)
cd4 <- list()
for (i in seq_files) {
  print(i)
  cd4[[i]] <- proc_data(i)
}

names(cd4) <- meta$sample_id

cd4_seu <- list()
for (i in names(cd4)) {
  print(i)
  cd4_seu[[i]] <- seu_obj(cd4[[i]])
}

# in a previous round of processing, non T cells were identified by a lack of cd3
bad_cells <- read.csv("exclude_cells.csv")$exclude
bad_cells <- data.frame(cell = gsub(bad_cells, pattern = "_[0-9]+", replacement = ""),
                        sample = str_extract(string = bad_cells, pattern = "[0-9]+") %>% as.numeric()) %>%
  group_split(sample)

`%notin%` <- Negate(`%in%`)

for (i in 1:length(cd4_seu)) {
  cd4_seu[[i]] <- cd4_seu[[i]][, colnames(cd4_seu[[i]]) %notin% bad_cells[[i]]$cell]
}


# sctransform
cd4_seu <- lapply(cd4_seu, function(x) {
  x <- SCTransform(x, method = "glmGamPoi")
})


# annotate data
tissue_info <- meta$tissue
names(tissue_info) <- meta$sample_id

stim_info <- meta$stimulation
names(stim_info) <- meta$sample_id

for (i in names(cd4_seu)) {
  cd4_seu[[i]]$stimulation <- stim_info[i]
  cd4_seu[[i]]$tissue <- tissue_info[i]
}


# integrate data
features <- SelectIntegrationFeatures(cd4_seu)
cd4_seu <- PrepSCTIntegration(cd4_seu, anchor.features = features)

cd4_seu <- lapply(cd4_seu, function(x) {
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(cd4_seu, 
                                  reduction = "rpca", 
                                  anchor.features = features,
                                  normalization.method = "SCT")

int_data <- IntegrateData(anchors, normalization.method = "SCT")

saveRDS(int_data, "rpca_integrated_tissue_finalised_2_sct.rds")


# load in and process integrated data
df <- readRDS("rpca_integrated_tissue_finalised_2_sct.rds")

# normalise rna assay
df <- NormalizeData(df, assay = "RNA")

df$tis_stim <- paste(df$tissue, df$stimulation)

# run dimred and visulalisation
df <- RunPCA(df, verbose = F)
ElbowPlot(df, ndims = 50)
df <- RunUMAP(df, dims = 1:15)
df <- FindNeighbors(df, dims = 1:15)
df <- FindClusters(df, resolution = 0.5)

# find all markers
DefaultAssay(df) <- "RNA"

mks <- FindAllMarkers(df, only.pos = T, logfc.threshold = 0.4)

write.csv(mks, file = "cluster_markers.csv", row.names = F)

# define populations
df$broad <- case_when(df$seurat_clusters %in% c(3, 4, 9, 10, 12) ~ "cd8",
                      df$seurat_clusters %in% c(0, 1, 2, 6, 7, 8, 11, 13, 14) ~ "cd4",
                      df$seurat_clusters %in% 5 ~ "treg")

df$fine <- case_when(df$seurat_clusters == 0 ~ "cd4_tcm",
                     df$seurat_clusters == 1 ~ "cd4_tem_1",
                     df$seurat_clusters == 2 ~ "cd4_tem_2",
                     df$seurat_clusters == 3 ~ "cd8_temra_1",
                     df$seurat_clusters == 4 ~ "cd8_tem",
                     df$seurat_clusters == 5 ~ "treg_1",
                     df$seurat_clusters == 6 ~ "cd4_resting",
                     df$seurat_clusters == 7 ~ "cd4_th1",
                     df$seurat_clusters == 8 ~ "cd4_tem_3",
                     df$seurat_clusters == 9 ~ "cd8_temra_2",
                     df$seurat_clusters == 10 ~ "cd8_tc1",
                     df$seurat_clusters == 11 ~ "cd4_ifn",
                     df$seurat_clusters == 12 ~ "cd8_resting",
                     df$seurat_clusters == 13 ~ "cd4_activated",
                     df$seurat_clusters == 14 ~ "treg_2")

saveRDS(df, file = "tissue_final.rds")

# make names nicer 
df$fine %>% unique
fine_names <- df$fine
fine_neat <- gsub(fine_names, pattern = "_", replacement = " ") %>%
  str_to_title() %>%
  gsub(pattern = "Cd4", replacement = "CD4") %>%
  gsub(pattern = "Cd8", replacement = "CD8") %>%
  gsub(pattern = "_", replacement = " ") %>%
  gsub(pattern = "Temra", replacement = "TEMRA") %>%
  gsub(pattern = "Ifn", replacement = "IFN")

df$fine_neat <- fine_neat



# 2.  plot umaps ----------------------------------------------------------
DimPlot(df, group.by = "stimulation", split.by = "tissue", ncol = 2) +
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_blank(),
        strip.text = element_blank()) +
  NoLegend()



# 3.  plot heatmap of markers ---------------------------------------------
mks <- read.csv("cluster_markers.csv")

markers <- data.frame("0" = c("CCR7", "MAL", "EIF3E"),
                      "1" = c("ENO1", "PKM", "LGALS1"),
                      "2" = c("MDM4", "CDC42SE1", "MALAT1"),
                      "3" = c("GZMB", "CCL5", "LAG3"),
                      "4" = c("IL7R", "IL7R", "S100A4"),
                      "5" = c("FOXP3", "TIGIT", "IL2RA"),
                      "6" = c("IL7R", "SELL", "KLF2"),
                      "7" = c("IFNG", "LTA", "CCL20"),
                      "8" = c("KLRB1", "IL7R", "LTB"),
                      "9" = c("GZMK", "NKG7", "EOMES"),
                      "10" = c("XCL1", "CCL3", "IFNG"),
                      "11" = c("IFI6", "IFIT3", "OAS1"),
                      "12" = c("REG4", "CCR7", "GAS5"),
                      "13" = c("CD69", "MYC", "IRF4"),
                      "14" = c("TNFRSF4", "FOXP3", "CTLA4")) %>%
  set_colnames(gsub(x = colnames(.), pattern = "X", replacement = "")) %>%
  pivot_longer(cols = current_vars(), values_to = "gene", names_to = "cluster") %>%
  mutate(cluster = as.numeric(cluster)) %>%
  arrange(cluster) %>%
  data.frame() %>%
  mutate(cell = rep(c("cd4_tcm", "cd4_tem_1", "cd4_tem_2", "cd8_temra_1", "cd8_tem",
                      "treg_1", "cd4_resting", "cd4_th1", "cd4_tem_3", "cd8_temra_2",
                      "cd8_tc1", "cd4_ifn", "cd8_resting", "cd4_activated", "treg_2"), each = 3))

# function for pulling gene expression
get_exprs <- function(cell_type) {
  
  df_sub <- subset(df, fine == cell_type)
  exprs <- rowMeans(df_sub@assays$RNA@data) %>%
    data.frame() %>%
    set_colnames(cell_type) %>%
    rownames_to_column("gene") %>%
    filter(gene %in% markers$gene)
  
  return(exprs)
}

# pull expression data
cell_types <- df$fine %>% unique

exprs_vals <- list()
for (i in cell_types) {
  print(i)
  exprs_vals[[i]] <- get_exprs(i)
}


s_data <- reduce(exprs_vals, full_join, "gene") %>%
  inner_join(., markers, "gene") %>%
  arrange(cluster) %>%
  mutate(gene = make.names(gene, unique = T)) %>%
  column_to_rownames("gene") %>%
  select(-cluster, -cell) %>%
  t() %>%
  scale() %>%
  data.frame()

# tidy up names
cell_names <- markers$cell %>%
  sub(pattern = "_", replacement = " ") %>%
  str_to_title() %>%
  sub(pattern = "Cd", replacement = "CD") %>%
  gsub(pattern = "Temra", replacement = "TEMRA") %>%
  gsub(pattern = "Ifn", replacement = "IFN") %>%
  sub(pattern = "_", replacement = " ") %>%
  unique()

# plot heatmap
s_data[markers$cell %>% unique(), ] %>%
  Heatmap(name = "Z-score",
          cluster_rows = F,
          cluster_columns = F, 
          border = T, 
          rect_gp = gpar(lwd = 0.2, col = 'white'),
          row_labels = cell_names,
          row_names_gp = gpar(fontsize = 9),
          column_names_gp = gpar(fontsize = 9))


# 4. comparisons across tissues ----------------------------------------------
# this analysis was not done locally -- done on nih cluster

# read in data
df <- readRDS("tissue_final.rds")

# pathways

pws <- c("kegg", "reactome", "biocarta", "wiki", "pid")
pathways <- msigdbr("Homo sapiens") %>%
  filter(grepl(paste(pws, collapse = "|"), gs_subcat, ignore.case = T) |
           grepl("HALLMARK", x = gs_name, ignore.case = T)) %>%
  format_pathways()

# across tissue resting 

cell_types <- unique(df$fine)

split_tissue <- SplitObject(df, split.by = "tissue")

bl_bm <- list(); bl_ln <- list(); bl_lung <- list()
for (i in cell_types) {
  
  blood <- seurat_extract(split_tissue$bl, 
                          meta1 = "fine", value_meta1 = i,
                          meta2 = "stimulation", value_meta2 = "none")
  
  bm <- seurat_extract(split_tissue$bm, 
                       meta1 = "fine", value_meta1 = i,
                       meta2 = "stimulation", value_meta2 = "none")
  
  ln <- seurat_extract(split_tissue$ln, 
                       meta1 = "fine", value_meta1 = i,
                       meta2 = "stimulation", value_meta2 = "none")
  
  lung <- seurat_extract(split_tissue$lung, 
                         meta1 = "fine", value_meta1 = i,
                         meta2 = "stimulation", value_meta2 = "none")
  
  print(paste("comparing", i))
  bl_bm[[i]] <- compare_pathways(list(blood, bm), pathways)
  bl_ln[[i]] <- compare_pathways(list(blood, ln), pathways)
  bl_lung[[i]] <- compare_pathways(list(blood, lung), pathways)
  
}

get_qvals <- function(scpa_out, name) {
  
  df <- list()
  for (i in names(scpa_out)) {
    df[[i]] <- scpa_out[[i]] %>%
      select(Pathway, qval)
  }
  
  col_names <- names(df)
  for (i in 1:length(df)) {
    df[[i]] <- set_colnames(df[[i]], c("pathway", paste(name, col_names[[i]], sep = "_")))
  }
  
  return(df)
  
}

pway <- Reduce(full_join, c(get_qvals(bl_bm, "bm"),
                            get_qvals(bl_ln, "ln"),
                            get_qvals(bl_lung, "lung")))

save(bl_bm, bl_lung, bl_ln, pway, file = "between_tissue_resting_scpa_output.RData")


# activated

df <- readRDS("tissue_final.rds")
print(DefaultAssay(df))

# pathways 

pws <- c("kegg", "reactome", "biocarta", "wiki", "pid")
pathways <- msigdbr("Homo sapiens") %>%
  filter(grepl(paste(pws, collapse = "|"), gs_subcat, ignore.case = T) |
           grepl("HALLMARK", x = gs_name, ignore.case = T)) %>%
  format_pathways()


# accross tissue stim

cell_types <- unique(df$fine)

split_tissue <- SplitObject(df, split.by = "tissue")

bl_bm <- list(); bl_ln <- list(); bl_lung <- list()
for (i in cell_types) {
  
  blood <- seurat_extract(split_tissue$bl, 
                          meta1 = "fine", value_meta1 = i,
                          meta2 = "stimulation", value_meta2 = "stim")
  
  bm <- seurat_extract(split_tissue$bm, 
                       meta1 = "fine", value_meta1 = i,
                       meta2 = "stimulation", value_meta2 = "stim")
  
  ln <- seurat_extract(split_tissue$ln, 
                       meta1 = "fine", value_meta1 = i,
                       meta2 = "stimulation", value_meta2 = "stim")
  
  lung <- seurat_extract(split_tissue$lung, 
                         meta1 = "fine", value_meta1 = i,
                         meta2 = "stimulation", value_meta2 = "stim")
  
  print(paste("comparing", i))
  bl_bm[[i]] <- compare_pathways(list(blood, bm), pathways)
  bl_ln[[i]] <- compare_pathways(list(blood, ln), pathways)
  bl_lung[[i]] <- compare_pathways(list(blood, lung), pathways)
  
}

get_qvals <- function(scpa_out, name) {
  
  df <- list()
  for (i in names(scpa_out)) {
    df[[i]] <- scpa_out[[i]] %>%
      select(Pathway, qval)
  }
  
  col_names <- names(df)
  for (i in 1:length(df)) {
    df[[i]] <- set_colnames(df[[i]], c("pathway", paste(name, col_names[[i]], sep = "_")))
  }
  
  return(df)
  
}


pway <- Reduce(full_join, c(get_qvals(bl_bm, "bm"),
                            get_qvals(bl_ln, "ln"),
                            get_qvals(bl_lung, "lung")))

save(bl_bm, bl_lung, bl_ln, pway, file = "between_tissue_activated_scpa_out.RData")



# 5.  plot heatmap of all comparisons -------------------------------------
load(file = "between_tissue_resting_scpa_output.RData")
pway_rest <- pway
load("between_tissue_activated_scpa_out.RData")
pway_act <- pway

# tidy up output
pway_rest <- pway_rest %>%
  column_to_rownames("pathway") %>%
  set_colnames(paste("rest", colnames(.), sep = "_")) %>%
  rownames_to_column("pathway")

pway_act <- pway_act %>%
  column_to_rownames("pathway") %>%
  set_colnames(paste("stim", colnames(.), sep = "_")) %>%
  rownames_to_column("pathway")

# combine data
all_data <- full_join(pway_rest, pway_act, "pathway") %>%
  column_to_rownames("pathway")

# define 
stim <- all_data %>%
  colnames() %>%
  substr(1, 4) %>%
  str_to_sentence()

tissue <- colnames(all_data) %>%
  substr(6, 9) %>%
  sub(pattern = "_[a-z]", replacement = "") %>%
  gsub(pattern = "bm", replacement = "BM") %>%
  gsub(pattern = "ln", replacement = "LN") %>%
  gsub(pattern = "lung", replacement = "Lung")

lineage <- colnames(all_data) %>%
  str_extract(pattern = "cd[48]") %>%
  str_to_upper() %>%
  replace_na("CD4")

col_an <- HeatmapAnnotation(Stimulation = stim,
                            Tissue = tissue,
                            Lineage = lineage,
                            col = list(Stimulation = c("Rest" = "gray70", "Stim" = "orangered2"),
                                       Tissue = c("BM" = "#cccccc", "LN" = "#be83e6", "Lung" = "#84c476"),
                                       Lineage = c("CD4" = "#4589ff", "CD8" = "#ff6363")),
                            gp = gpar(col = "white", lwd = 0.05),
                            annotation_name_gp = gpar(fontsize = 9),
                            simple_anno_size = unit(3, "mm"))

hm <- all_data %>%
  Heatmap(name = "Qval",
          show_row_names = F, 
          top_annotation = col_an,
          border = T,
          show_row_dend = F,
          show_column_dend = F,
          show_column_names = F)

hm <- draw(hm)



# pca of qvals ------------------------------------------------------------
pca_calc <- all_data %>%
  t() %>%
  prcomp()

pca_plot <- data.frame(pc1 = pca_calc$x[, 1],
                       pc2 = pca_calc$x[ ,2]) %>%
  rownames_to_column("cell") %>%
  mutate(stim = substr(cell, 1, 4)) %>%
  mutate(tissue = str_extract(cell, "bm|ln|lung")) %>%
  mutate(lineage = str_extract(cell, "cd[48]") %>% 
           str_to_upper() %>%
           replace_na("CD4"))


pcaVar <- pca_calc$sdev^2
pcaVarPercentage <- round(pcaVar/sum(pcaVar)*100, 1)

ggplot(pca_plot, aes(pc1, pc2, fill = tissue, shape = stim)) +
  geom_point(cex = 2.8) +
  scale_shape_manual(values = c(22, 21)) +
  scale_fill_discrete(labels = c("bm" = "Bone marrow", "ln" = "Lymph node", "lung" = "Lung")) +
  scale_fill_manual(values = c("bm" = "#cccccc", "ln" = "#be83e6", "lung" = "#84c476")) +
  labs(fill = "Tissue", shape = "Stimulation", 
       x = paste("PC1:", pcaVarPercentage[1], "%"), 
       y = paste("PC2:", pcaVarPercentage[2], "%")) +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.key = element_blank())



# qvals with highest variance ---------------------------------------------

apply(all_data, 1, var) %>%
  data.frame() %>% 
  set_colnames("variation") %>%
  arrange(desc(variation)) %>% 
  rownames_to_column("pathway") %>%
  ggplot(aes(reorder(pathway, variation), variation)) +
  geom_point(shape = 21, cex = 3, fill = "royalblue2", color = 'black', stroke = 0.2) +
  scale_x_discrete(expand = c(0.04, 0.04)) +
  labs(x = "Pathway", y = "Variance") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA))



# plot defa data ----------------------------------------------------------
plots <- VlnPlot(df, features = c("DEFA1", "DEFA3"), pt.size = 0, group.by = "tissue", combine = F)

p1 <- VlnPlot(df, "DEFA1", pt.size = 0, group.by = "tissue") +
  theme(plot.title = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none",
        axis.title = element_blank())

p2 <- VlnPlot(df, "DEFA3", pt.size = 0, group.by = "tissue") +
  scale_x_discrete(labels = c("bl" = "Blood", "bm" = "BM", "ln" = "LN", "lung" = "Lung")) +
  theme(plot.title = element_blank(),
        legend.position = "none",
        axis.title = element_blank())

patchwork::wrap_plots(p1, p2, ncol = 1)



# defa microarray data ----------------------------------------------------
df <- read.table("GSE50677_series_matrix.txt", sep = "\t", header = T)

sample_location <- str_extract(string = colnames(df), pattern = "bone|peripheral")
sample_location <- sample_location[!is.na(sample_location)]

df <- set_colnames(df, make.names(c("PROBEID", sample_location), unique = T))

gene_names <- select(hgu133plus2.db, keys = df$PROBEID, keytype = "PROBEID", columns = "SYMBOL")

df <- full_join(df, gene_names, "PROBEID") %>%
  dplyr::select(-PROBEID) %>% 
  dplyr::select(SYMBOL, current_vars())

df <- clean_names(df)

defensin <- df[df$symbol %in% grep(x = df$symbol, pattern = "DEFA1$", value = T),]

defensin %>%
  remove_rownames() %>%
  column_to_rownames("symbol") %>%
  set_colnames(sample_location) %>%
  pivot_longer(cols = current_vars(), values_to = "expression", names_to = "location") %>%
  mutate(location = factor(location, levels = c("peripheral", "bone"))) %>%
  mutate(expression = expression + 1) %>%
  ggplot(aes(location, log2(expression))) +
  scale_x_discrete(labels = c("peripheral" = "Blood", "bone" = "BM")) +
  scale_y_continuous(breaks = c(0, 5, 10, 15), limits = c(0, 17)) +
  geom_boxplot(fill = c("gray80", "#B5B6FA"), alpha = 0.3) +
  scale_fill_manual(values = c("gray80", "#B5B6FA")) +
  geom_jitter(aes(fill = location), shape = 21, cex = 2.5, width = 0.2) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 10)) +
  labs(y = "DEFA1 expression (log2)")


stats_test <- defensin %>%
  remove_rownames() %>%
  column_to_rownames("symbol") %>%
  set_colnames(sample_location) %>%
  pivot_longer(cols = current_vars(), values_to = "expression", names_to = "location")


t.test(stats_test %>% filter(location == "peripheral") %>% pull(expression),
       stats_test %>% filter(location == "bone") %>% pull(expression))




