
# packages ----------------------------------------------------------------
library(SCPA)
library(Seurat)
library(msigdbr)
library(tidyverse)
library(magrittr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(patchwork)
library(drc)


# All samples -------------------------------------------------------------
load("int_data_N4.RData")
load("int_data_N8.RData")
load("int_data_M4.RData")
load("int_data_M8.RData")

# 1.  comparison over three time points -----------------------------------
# extract hours -----------------------------------------------------------
int_data_N4 <- SplitObject(int_data_N4, "Hour") %>%
  lapply(function(x) seurat_extract(x))

int_data_M4 <- SplitObject(int_data_M4, "Hour") %>%
  lapply(function(x) seurat_extract(x))

int_data_N8 <- SplitObject(int_data_N8, "Hour") %>%
  lapply(function(x) seurat_extract(x))

int_data_M8 <- SplitObject(int_data_M8, "Hour") %>%
  lapply(function(x) seurat_extract(x))


# compare across time -----------------------------------------------------
pathways <- "hm_react_kegg_combined.csv"

n4 <- compare_pathways(int_data_N4, pathways)
n8 <- compare_pathways(int_data_M4, pathways)
m4 <- compare_pathways(int_data_N8, pathways)
m8 <- compare_pathways(int_data_M8, pathways)

hm_mat <- n4 %>%
  select(Pathway, qval) %>%
  rename(qval_n4 = qval) %>%
  full_join(., m4 %>% select(Pathway, qval) %>%
              rename(qval_m4 = qval) %>%
              full_join(., n8 %>%
                          select(Pathway, qval) %>%
                          rename(qval_n8 = qval) %>%
                          full_join(., m8 %>%
                                      select(Pathway, qval) %>%
                                      rename(qval_m8 = qval)))) %>%
  column_to_rownames("Pathway")
save(hm_mat, file = "hm_mat.RData")

col_hm <- colorRamp2(colors = c("white", "gray85", "royalblue2"), breaks = c(0, 5, 10))
hm <- Heatmap(hm_mat,
              show_row_names = F,
              row_km = 4,
              border = T,
              col = col_hm,
              column_labels = c("Naive CD4", "Memory CD4", "Naive CD8", "Memory CD8"),
              column_names_gp = gpar(fontsize = 9))
ht <- draw(hm)

# 2. boxplots of qval distributions ------------------------------------------
load("hm_mat.RData")

hm_mat %>% 
  rownames_to_column("pathway") %>%
  set_colnames(c("pathway", "N4", "M4", "N8", "M8")) %>%
  pivot_longer(cols = N4:M8, names_to = "cell", values_to = "qval") %>%
  mutate(cell = factor(cell, levels = c("N4", "M4", "N8", "M8"))) %>%
  
  ggplot(aes(cell, qval)) +
  geom_boxplot(notch = T, outlier.shape = 21, outlier.fill = 'gray60',
               color = 'black',
               fill = c("#eda193", "#eb7059", "#c4d7f2", "#84b0f0")) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10)) +
  ylab("Qval")

# 3. heatmap of cytokine response --------------------------------------------
load("hm_mat.RData")

pways <- read.csv("pathway_quant.csv") %>%
  filter(role == "cytokine_response") %>%
  pull(pathway)

cytokine_names <- c("IFN", "General", "IFN", "IFN", "General", "TNF", "General", "IFN",
                    "IL4/13", "IFN", "IL2", "IL12", "General", "IL6", "IL1", "IL1",
                    "General", "IFN", "IL12")

cell_names <- c("nCD4", "mCD4", "nCD8", "mCD8")

hm_col <- colorRamp2(colors = c("cornflowerblue", "white", "orangered2"), breaks = c(50, 100, 200))

hm_mat %>%
  rownames_to_column("pathway") %>%
  filter(pathway %in% pways) %>%
  column_to_rownames("pathway") %>%
  .^2 %>%
  Heatmap(name = "-log10(FDR)",
          show_row_names = T,
          border = T,
          rect_gp = gpar(lwd = 0.2, col = "white"),
          cluster_columns = F,
          col = hm_col,
          row_labels = cytokine_names,
          column_labels = cell_names)



# 4. highlight ifn pathways in naive cd4 -------------------------------------
pws <- c("INTERFE", "IFN")
ifn_paths <- hm_mat %>%
  rownames_to_column("pathway") %>%
  filter(grepl(paste(pws, collapse = "|"), pathway)) %>%
  .[1:6, ] %>%
  pull(pathway)

ifn_rank <- hm_mat %>%
  rownames_to_column("pathway") %>%
  arrange(desc(qval_n4)) %>%
  mutate(path_rank = percent_rank(qval_n4))

ggplot(ifn_rank, aes(qval_n4, path_rank*100)) +
  geom_hline(yintercept = c(0, 25, 50, 75, 100), linetype = 'dotted', lwd = 0.3, color = 'gray40') +
  geom_point(shape = 21, cex = 2, color = 'black', fill = 'royalblue2', stroke = 0.05) +
  geom_point(data = subset(ifn_rank, pathway %in% ifn_paths), shape = 21, cex = 3, color = 'black', fill = 'orangered2') +
  xlab("Qval") +
  ylab("Pathway rank") +
  ggtitle("Naive CD4 T cells") +
  scale_y_continuous(expand = c(0.03, 0.03), breaks = c(0, 25, 50, 75, 100)) +
  scale_x_continuous(expand = c(0.2, 0.2)) +
  theme(panel.border = element_rect(fill = NA),
        panel.background = element_blank(),
        title = element_text(size = 9),
        axis.title = element_text(size = 10))



# 5. heatmap of ifn genes across populations and stim ------------------------
# function for extracting gene average ------------------------------------
path_avg <- function(obj, name) {
  df1 <- obj[[1]] %>%
    data.frame() %>%
    .[rownames(.) %in% ifn_genes, ] %>%
    rowMeans() %>%
    data.frame %>%
    set_colnames(paste(name, "0", sep = "_")) %>%
    rownames_to_column("gene")
  df2 <- obj[[2]] %>%
    data.frame() %>%
    .[rownames(.) %in% ifn_genes, ] %>%
    rowMeans() %>%
    data.frame %>%
    set_colnames(paste(name, "24", sep = "_")) %>%
    rownames_to_column("gene")
  
  df <- full_join(df1, df2, "gene")
  
  return(df)
}

ifn_hm <- Reduce(full_join, list(
  path_avg(list(rest, act), "naive"),
  path_avg(list(rest_treg, act_treg), "naive_treg"),
  path_avg(th1, "th1"),
  path_avg(th2, "th2"),
  path_avg(treg, "treg"),
  path_avg(tcm, "tcm")))

ifn_hm <- ifn_hm %>%
  replace(is.na(.), 0.001) %>%
  column_to_rownames("gene") %>%
  t() %>%
  scale() %>%
  t() %>%
  data.frame() %>%
  drop_na()

# row annotation
highlight_genes <- c("GBP1", "GBP2", "IFI35", "IRF1", "STAT1", "GBP4", "JAK2",
                     "NUP35", "NUP54", "EIF4E2", "NUP85", "PTPN11", "IRF4", "IRF8",
                     "IFI30", "ADAR", "STAT2", "HLA-DRB1", "ISG20", "SOCS3", "IRF9", "ISG15", "IRF7",
                     "JAK1", "IFNAR2", "IRF3", "TRIM2", "EIF4A2",
                     "IRF6", "UBC", "TRIM45")
position <- which(rownames(ifn_hm) %in% highlight_genes)
row_an <- rowAnnotation(Genes = anno_mark(at = position,
                                          labels = rownames(ifn_hm)[position],
                                          labels_gp = gpar(fontsize = 6),
                                          link_width = unit(2.5, "mm"),
                                          padding = unit(1, "mm"),
                                          link_gp = gpar(lwd = 0.5)))

# column annotation
col_cols <- rep(c("0", "24"), times = 6)
col_an <- HeatmapAnnotation(Hour = col_cols,
                            col = list(Hour = c("0" = "gray60", "24" = "purple")),
                            gp = gpar(col = "white", lwd = 0.2),
                            simple_anno_size = unit(3, "mm"), 
                            annotation_name_gp = gpar(fontsize = 8))
col_names <- rep(c("Naive", "Treg", "Th1", "Th2", "mTreg", "Tcm"), each = 2)

png(filename = "test.png", res = 600, width = 6, height = 5.5, units = "in")
Heatmap(ifn_hm,
        column_names_gp = gpar(fontsize = 7.5),
        border = T,
        rect_gp = gpar(col = 'white', lwd = 0.1),
        row_km = 5,
        column_km = 2,
        top_annotation = col_an,
        right_annotation = row_an,
        show_row_names = F,
        column_names_side = "top",
        column_labels = col_names,
        column_dend_height = unit(3, "mm"))
dev.off()



# 6. plot qpcr data for ifng genes -------------------------------------------
df <- read.csv("ifn_analysis.csv") %>%
  mutate(hour = as.character(hour)) %>%
  mutate(target = toupper(target)) %>%
  mutate(target = factor(as.character(target), levels = c("IFNA1/13/21", "IFNA4", "IFNA7", "IFNA17")))


point_col <- rep(c("royalblue2", "tomato", "mediumpurple3", "gray60"), each = 12)
box_col <- rep(c("royalblue2", "tomato", "mediumpurple3", "gray60"), each = 3)

ggplot(df, aes(hour, final_val)) +
  geom_boxplot(outlier.shape = NA, fill = box_col, alpha = 0.2, lwd = 0.45) +
  geom_jitter(shape = 21, fill = point_col, width = 0.2, cex = 2, stroke = 0.3, height = 0) +
  facet_wrap(~target, ncol = 2, scales = "free_y") +
  scale_y_continuous(limits = function(x){c(0, max(0.001, x))}) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.title.x = element_blank(),
        strip.text = element_text(face = "italic")) +
  ylab("Expression relative to RPL7")

ggsave("ifn_expression.png",
       height = 4, width = 3.8, units = "in", device = "png", dpi = 600)

# and get ct values for each
df %>%
  group_by(target) %>%
  dplyr::summarize(sd = sd(ct_mean, na.rm=TRUE))



# 7. plot ifn elisa ----------------------------------------------------------
df <- read.csv("values.csv", header = T) %>%
  data.frame() %>%
  rename(Standard = std) %>%
  data.frame() %>%
  mutate(Conc = c(125, 62.5, 31.25, 15.63, 7.81, 3.91, 1.95, 0))

plot(df$Standard, df$Conc)

r_sq <- cor(df$Standard, df$Conc)^2

lm(df$Standard ~ df$Conc)

model_1 <- drm(data = df, 
               Standard~Conc, 
               fct = LL.4())

dat <- df %>%
  dplyr::select(d1:d3) %>%
  mutate(hour = c(0, 1, 3, 6, 8, 10, 12, 24)) %>%
  pivot_longer(cols = d1:d3, names_to = "donor", values_to = "od") %>%
  arrange(donor)

dat <- dat %>%
  mutate(conc = ED(model_1, dat$od, type = "absolute", display = T)[, "Estimate"]) %>%
  mutate(hour = factor(as.character(hour), levels = c("0", "1", "3", "6", "8", "10", "12", "24")))

dat$conc[is.nan(dat$conc)] <- 0

dat %>%
  filter(hour %in% c("0", "3", "10", "12", "24")) %>%
  ggplot(aes(hour, conc)) +
  geom_hline(yintercept = 1.95, col = "red", linetype = "dashed") +
  geom_boxplot(lwd = 0.4, fill = "royalblue2", alpha = 0.2) +
  geom_jitter(cex = 2.5, shape = 21, fill = "royalblue2", width = 0.15, height = 0.1, stroke = 0.3) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "none") +
  ylab("Concentration (pg/ml)") +
  xlab("Hour")

ggsave("ifn_secretion.png",
       width = 2.5, height = 4, units = "in", device = "png", dpi = 600)



# 8. stat1 pseudotime --------------------------------------------------------
# get stat1 gene sets from msigdb
stat <- msigdbr("Homo sapiens", "C3") %>%
  filter(grepl(pattern = "STAT1_", gs_name)) %>%
  pull(gene_symbol)

# load in pseudotime calculations
load("int_data_N4.RData")
load(file = "NaiveCD4_ti_slingshot_model_1000var_feat.RData")

int_data_N4 <- subset(int_data_N4, idents = "Treg", invert = T)
hm_dat <- seurat_extract(int_data_N4)

hm_dat <- hm_dat[rownames(hm_dat) %in% stat, ]

hm_dat <- hm_dat[rowSums(hm_dat) > 1000, ]

# scale data
hm_scale <- t(scale(t(as.matrix(hm_dat))))

# column annotations
cd4_pseudo <- calculate_pseudotime(model_n4)

colors_cd4_pseudo <- colorRamp2(colors = viridis(100),
                                breaks = seq(0, max(cd4_pseudo), max(cd4_pseudo)/99))

col_an <- HeatmapAnnotation(Pseudotime = cd4_pseudo[order(cd4_pseudo)],
                            col = list(Pseudotime = colors_cd4_pseudo))


# 8. t7 quantification ----------------------------------------------------
col_point <- rep(c("gray60", "royalblue2"), times = 4)

data.frame(mock = c(2.5, 0, 4.5, 0),
           ifn = c(9.2, 37.8, 38.4, 89.6)) %>%
  pivot_longer(cols = mock:ifn, values_to = "mod", names_to = "treatment") %>%
  mutate(treatment = factor(treatment, levels = c("mock", "ifn"))) %>%
  ggplot(aes(treatment, mod)) +
  geom_boxplot(outlier.colour = "transparent") +
  geom_jitter(cex = 2.3, width = 0.2, shape = 21, fill = col_point) +
  scale_x_discrete(labels = c("mock" = "Mock", "ifn" = expression(italic("IFNA")))) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.title.x = element_blank()) +
  labs(y = "% Modification")

ggsave("t7_quant.png", device = "png", dpi = 600,
       height = 2, width = 1.7, units = "in")


# 9.  ifna crispr viability quantification --------------------------------
data.frame(scr = c(19, 6, 8, 19, 18),
           ifn = c(32, 41, 20, 38, 51)) %>%
  mutate(scr = 100 - scr) %>%
  mutate(ifn = 100 - ifn) %>%
  pivot_longer(cols = scr:ifn, names_to = "condition", values_to = "dead_cells") %>%
  mutate(condition = factor(condition, levels = c("scr", "ifn"))) %>%
  mutate(donor = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)) %>%
  ggplot(aes(condition, dead_cells)) +
  geom_boxplot() +
  geom_jitter(shape = 21, cex = 3, fill = rep(c("tomato", "royalblue2"), times = 5),
              width = 0.2) +
  scale_x_discrete(labels = c("scr" = "Mock", "ifn" = expression(italic("IFNA")))) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 100)) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10)) +
  ylab("% Viability")

ggsave(filename = "ifn_viability_human.png",
       height = 3.5, width = 2.2, device = "png", dpi = 600)


# 10.  ifnar1-/- t cell viability quantification --------------------------
# summary data
df <- read.csv("summary_stats.csv") %>%
  mutate(treatment = factor(treatment, levels = c("Unstim", "18hr", "24hr", "48hr", "72hr"))) %>%
  mutate(genotype = factor(genotype, levels = c("wt", "ko")))

point_col <- rep(c("gray60", "tomato"), each = 5)

# plot
ggplot(df, aes(treatment, live_mean, group = genotype, color = genotype)) +
  geom_line() +
  geom_errorbar(aes(ymin = live_mean - live_sd, ymax = live_mean + live_sd, width = 0.12)) +
  geom_point(shape = 21, color = "black", fill = point_col, cex = 3, stroke = 0.3) +
  scale_color_manual(name = "Genotype", values = c("wt" = "gray60", "ko" = "tomato"),
                     label = c("wt" = "WT", "ko" = bquote("IFNAR1"^"-/-"))) +
  ylim(0, NA) +
  labs(y = "% Viability") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.key = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

ggsave("ifn_viability_mouse.png",
       height = 3.5, width = 2.5, units = "in", device = "png", dpi = 600)











