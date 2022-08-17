library(Seurat)
library(tidyverse)
library(dyno)
library(magrittr)
library(SCPA)
library(ComplexHeatmap)
library(msigdbr)
library(SCPA)


# 1. trajectory analysis --------------------------------------------------
load("int_data_N4.RData")
int_data_N4 <- subset(int_data_N4, idents = "Treg", invert = T)
df <- as.matrix(int_data_N4[["RNA"]]@data)
var_genes <- names(sort(apply(df, 1, var), decreasing = TRUE))[1:1000]

counts <- Matrix::t(as(as.matrix(int_data_N4@assays$RNA@counts[var_genes,]), 'sparseMatrix'))
expression <- Matrix::t(as(as.matrix(int_data_N4@assays$RNA@data[var_genes,]), 'sparseMatrix'))

dataset_n4 <- wrap_expression(expression = expression,
                              counts = counts)

dataset_n4 <- add_grouping(dataset_n4,
                           grouping = int_data_N4$Cell_Type)

dataset_n4 <- add_dimred(dataset_n4,
                         Embeddings(int_data_N4, "umap"))

model_n4 <- infer_trajectory(dataset_n4, method = ti_slingshot(), verbose = T)

plot_dimred(model_n4, grouping = group_onto_nearest_milestones(model_n4), hex_cells = F,
            plot_trajectory = T, size_cells = 1, alpha_cells = 0.8) +
  theme(aspect.ratio = 1)

ggsave(file = "Milestone_CD4.png", device = "png", width = 4.5, height = 4.5, units = "in")



# 2.  test pathways over nodes --------------------------------------------
# get groups
mile_group <- data.frame(group_onto_nearest_milestones(model_n4)) %>%
  set_colnames("Milestone") %>%
  rownames_to_column("Cell")

n4_seu <- subset(int_data_N4, idents = "Treg", invert = T)
n4_seu$Milestone <- mile_group$Milestone

mstones <- list()
for (i in 1:max(as.numeric(mile_group$Milestone))) {
  mstones[[i]] <- seurat_extract(n4_seu, meta1 = "Milestone", value_meta1 = i)
}



# run scpa
pathways <- "Combined_metab_paths_noGO.csv"

mstone_out <- compare_pathways(mstones, pathways, downsample = 250)

mstone_out <- mstone_out %>%
  data.frame() %>%
  select(Pathway, qval) %>%
  remove_rownames() %>%
  column_to_rownames("Pathway")

col_hm <- colorRamp2(colors = c("white", "red"), breaks = c(0, max(mstone_out)))

png("Heatmap_pseudotime.png", width = 7, height = 1.22, units = "in", res = 600)
Heatmap(t(mstone_out),
        name = "Qvalue",
        col = col_hm,
        border = T,
        rect_gp = gpar(col = "white", lwd = 0.1),
        show_column_dend = F,
        show_row_names = F,
        show_column_names = F)
dev.off()




# 3. highlight glycolysis and lin acid pathway expression -----------------
gly_genes <- msigdbr("Homo sapiens", category = "H") %>%
  filter(gs_name == "HALLMARK_GLYCOLYSIS") %>%
  pull(gene_symbol)
lino_genes <- msigdbr("Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
  filter(gs_name == "KEGG_LINOLEIC_ACID_METABOLISM") %>%
  pull(gene_symbol)

df <- int_data_N4@assays$RNA@data
gly_dat <- df[rownames(df) %in% gly_genes, ] %>%
  as.matrix() %>%
  colMeans(.) %>%
  data.frame() %>%
  set_colnames("Gly") %>%
  rownames_to_column("Cell")

lino_dat <- df[rownames(df) %in% lino_genes, ] %>%
  as.matrix() %>%
  colMeans(.) %>%
  data.frame() %>%
  set_colnames("Lino") %>%
  rownames_to_column("Cell")

gly_lino <- full_join(gly_dat, lino_dat, "Cell")

pseudo <- calculate_pseudotime(model_n4) %>%
  data.frame() %>%
  set_colnames("Pseudotime") %>%
  rownames_to_column("Cell")

gly_lino_pseudo <- left_join(gly_lino, pseudo, "Cell")

ggplot(gly_lino_pseudo, aes(x = Pseudotime)) +
  geom_point(aes(y = Gly), cex = 0.8, shape = 21, fill = 'red', color = 'black', alpha = 0.6, stroke = 0.2) +
  geom_smooth(aes(y = Gly), col = 'orangered2', lwd = 0.6) +
  geom_point(aes(y = Lino), cex = 0.8, shape = 21, fill = 'seagreen3', color = 'black', alpha = 0.6, stroke = 0.2) +
  geom_smooth(aes(y = Lino), col = 'seagreen3', lwd = 0.6) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA)) +
  labs(y = "Pathway Expression")

ggsave("Pathway_pseudotime.png", device = "png", width = 3.8, height = 2.3, units = "in", dpi = 600)
dev.off()



# 4. compare metablic pathways over real time -----------------------------
rest <- seurat_extract(int_data_N4,
                       meta1 = "Cell_Type", value_meta1 = "Resting",
                       meta2 = "Hour", value_meta2 = 0)
act <- seurat_extract(int_data_N4,
                      meta1 = "Cell_Type", value_meta1 = "Activated",
                      meta2 = "Hour", value_meta2 = 24)

rest_act <- compare_pathways(list(rest, act), pathways = "Combined_metab_paths_noGO.csv")

rest_act$FC <- -rest_act$FC
rest_act <- rest_act %>%
  mutate(color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
                           FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
                           FC < -5 & adjPval < 0.01 ~ 'seagreen2',
                           FC < 5 & FC > -5 & adjPval > 0.01 ~ 'black'))

aa_path <- rest_act %>% 
  filter(grepl(pattern = "ome_arachi", ignore.case = T, x = Pathway))

ggplot(rest_act, aes(FC, qval)) +
  geom_vline(xintercept = c(-5, 5), linetype = "dashed", col = 'black', lwd = 0.3) +
  geom_point(cex = 2.6, shape = 21, fill = rest_act$color, stroke = 0.3) +
  geom_point(data = aa_path, shape = 21, cex = 2.8, fill = "orangered2", color = "black", stroke = 0.3) +
  xlim(-20, 80) +
  ylim(0, 11) +
  xlab("Enrichment") +
  ylab("Qval") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1)

ggsave(filename = "V_plot_rest_act.png", device = "png",
       width = 2.3, height = 2.3, units = "in", dpi = 600)


# 5. highlight arachidonic acid enzymes ----------------------------------
aa <- msigdbr::msigdbr("Homo sapiens") %>%
  filter(grepl("arachidonic", gs_name, ignore.case = T))

aa_genes <- aa %>% pull(gene_symbol) %>% unique()

levels(int_data_N4) <- c("Resting", "Intermediate", "Activated", "Treg") %>%
  rev()

subset(int_data_N4, idents = "Treg", invert = T) %>%
  DotPlot(features = c("PTGES3", "LTA4H", "PTGES2", "CBR1", "GPX4", "ALOX5")) + 
  RotatedAxis() +
  NoLegend() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(face = "italic"))


# 6.  gsea aa -------------------------------------------------------------


# 7.  aa elisa ------------------------------------------------------------

df <- read.csv("aa_elisa.csv") %>%
  mutate(conc = conc*1000)

ggplot(df, aes(factor(hour), conc, fill = stim)) +
  geom_boxplot(lwd = 0.4) +
  geom_point(position = position_jitterdodge(), cex = 2.3, shape = 21, stroke = 0.4) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "none") +
  ylim(0, 2500) +
  scale_fill_manual(values = c("no" = "gray60", "yes" = "#9dc1fa")) +
  labs(x = "Hour", y = "Concentration pg/ml")

ggsave("arachidonic_acid/plot.png",
       device = "png", dpi = 600,
       width = 2.5, height = 2.5, units = "in")


# 8. cd69 pla2i -----------------------------------------------------------
df <- read.csv("flow_summary.csv")

boxplot_fill <- c("gray80", "gray80", RColorBrewer::brewer.pal(n = 6, "Blues"))
point_fill <- rep(c("gray80", "gray80", rev(RColorBrewer::brewer.pal(n = 6, "Blues"))), times = 3)

df %>%
  filter(treatment == "cay") %>%
  mutate(sample = as.character(sample)) %>%
  mutate(sample = factor(sample, levels = c("Unstim", "0", "0.1", "1", "10", "25", "50", "100"))) %>%
  ggplot(aes(sample, normalised)) +
  geom_boxplot(fill = boxplot_fill, alpha = 0.5, lwd = 0.4) +
  geom_jitter(width = 0.1, shape = 21, cex = 2.2, fill = point_fill, color = "black", stroke = 0.3) +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1.2, 0.2)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA)) +
  labs(y = "Normalised CD69 gMFI", x = "uM PLA2i")
ggsave("pla2_plot.png",
       width = 3, height = 3, device = "png", units = "in", dpi = 600)


# 9.  cba and viability ---------------------------------------------------
# cytokine
df <- read.csv("33070921 cba pla2i/data.csv") %>%
  mutate(sample = factor(as.character(sample), levels = c("na", "0", "1", "10"))) %>%
  mutate(live = c(0, 83, 93, 0, 0, 95, 83, 0, 0, 91, 95, 0))

blue_colors <- RColorBrewer::brewer.pal(8, "Blues")
boxplot_fill <- c("gray80", blue_colors[6])
point_fill <- rep(boxplot_fill, times = 3)
facet_labels <- c("IFNy", "IL-10", "TNFa", "% live")
names(facet_labels) <- c("ifng", "il10", "tnf", "live")

df %>%
  filter(sample %in% c("0", "10")) %>%
  pivot_longer(cols = ifng:live, values_to = "pg_ml", names_to = "cytokine") %>%
  filter(cytokine %in% c("ifng", "il10", "tnf")) %>%
  ggplot(aes(sample, pg_ml)) +
  geom_boxplot(fill = rep(boxplot_fill, times = 3), alpha = 0.4) +
  geom_jitter(width = 0.2, shape = 21, cex = 2, fill = rep(point_fill, each = 3)) +
  scale_x_discrete(labels = c("na" = "Unstim", "0" = "Ctrl", "10" = "PLA2i")) +
  facet_wrap(~cytokine, scales = "free_y", labeller = labeller(cytokine = facet_labels)) +
  ylim(0, NA) +
  labs(y = "Cytokine (pg/ml)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.background = element_rect(fill = "gray92", linetype = "solid", color = "black"))

ggsave("cytokine.png", device = "png",
       width = 4.5, height = 3, units = "in", dpi = 600)

# viability
df %>%
  filter(sample %in% c("0", "10")) %>%
  pivot_longer(cols = ifng:live, values_to = "pg_ml", names_to = "cytokine") %>%
  filter(cytokine %in% c("live")) %>%
  ggplot(aes(sample, pg_ml)) +
  geom_boxplot(fill = boxplot_fill, alpha = 0.4) +
  geom_jitter(width = 0.2, shape = 21, cex = 2, fill = point_fill) +
  scale_x_discrete(labels = c("na" = "Unstim", "0" = "Ctrl", "10" = "PLA2i")) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25)) +
  labs(y = "% Live cells") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.background = element_rect(fill = "gray92", linetype = "solid", color = "black"))

ggsave("viability.png", device = "png",
       width = 1.7, height = 3, units = "in", dpi = 600)


# 10. compare resting naive and memory ----------------------------------------
load("int_data_N4.RData")
load("int_data_M4.RData")

# get expression data
n4 <- seurat_extract(int_data_N4,
                     meta1 = "Cell_Type", value_meta1 = "Resting",
                     meta2 = "Hour", value_meta2 = 0)
m4 <- seurat_extract(int_data_M4,
                     meta1 = "Cell_Type", value_meta1 = "Tcm",
                     meta2 = "Hour", value_meta2 = 0)
th1 <- seurat_extract(int_data_M4,
                      meta1 = "Cell_Type", value_meta1 = "Th1",
                      meta2 = "Hour", value_meta2 = 0)
th2 <- seurat_extract(int_data_M4,
                      meta1 = "Cell_Type", value_meta1 = "Th2",
                      meta2 = "Hour", value_meta2 = 0)
treg <- seurat_extract(int_data_M4,
                       meta1 = "Cell_Type", value_meta1 = "Treg",
                       meta2 = "Hour", value_meta2 = 0)

# define pathways
pathways <- "Combined_metab_paths_noGO.csv"

# compare pathways
rest_tcm <- compare_pathways(samples = list(n4, m4), pathways = pathways, min_genes = 5)
rest_th1 <- compare_pathways(samples = list(n4, th1), pathways = pathways, min_genes = 5)
rest_th2 <- compare_pathways(samples = list(n4, th2), pathways = pathways, min_genes = 5)
rest_treg <- compare_pathways(samples = list(n4, treg), pathways = pathways, min_genes = 5)

# extract qvals
qvals <- function(object, name) {
  df <- object %>%
    select(Pathway, qval) %>%
    set_colnames(c("Pathway", name))
}

tcm <- qvals(rest_tcm, "Tcm")
th1 <- qvals(rest_th1, "Th1")
th2 <- qvals(rest_th2, "Th2")
treg <- qvals(rest_treg, "Treg")

qval_cell <- Reduce(full_join, list(tcm, th1, th2, treg))

# heatmap for plotting
metab_names <- c("Arachidonic acid", "Glycogen", "Glycogen", "Glycogen", "Glycogen",
                 "Amino acid", "Cholesterol", "Vitamins/cofactors", "Nucleotides",
                 "Purine", "Glycolysis", "Glycolysis", "Glycolysis", "Glycolysis",
                 "Fatty acid", "Metabolic genes", "Ployamines", "Glycolysis",
                 "ETC", "Glycolysis", "OXPHOS", "OXPHOS")

qval_cell %>%
  filter(Pathway != "REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES") %>%
  column_to_rownames("Pathway") %>%
  filter_all(any_vars(abs(.) > 5)) %>%
  Heatmap(name = "Qval",
          show_row_names = T,
          row_names_gp = gpar(fontsize = 9),
          border = T,
          row_labels = metab_names,
          rect_gp = gpar(col = 'white', lwd = 0.2),
          column_dend_height = unit(4, "mm"),
          row_dend_width = unit(4, "mm"))


# 11. highlight amino acid metab ----------------------------------------------
metab_path <- msigdbr("Homo sapiens") %>%
  filter(gs_name == "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES") %>%
  pull(gene_symbol)

metab_path <- metab_path[metab_path %in% rownames(int_data_M4)]

rpl_psm <- grep(x = rownames(int_data_M4), pattern = "^rpl|^psm", value = T, ignore.case = T) %>%
  sort()

test <- FindMarkers(int_data_M4, ident.1 = "Tcm", ident.2 = "Th1", features = metab_path %>% unique(), logfc.threshold = 0, min.pct = 0) %>%
  mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-300, p_val_adj)) %>%
  rownames_to_column("gene")

dot_col <- case_when(startsWith(x = test$gene, prefix = "PSM") == T ~ "red",
                     grepl(x = test$gene, pattern = "RPS|RPL") == T ~ "blue",
                     !grepl(x = test$gene, pattern = "PSM|RPL|RPS") == T ~ "gray")

ggplot(test, aes(-avg_log2FC, -log10(p_val_adj))) +
  geom_jitter(shape = 21, cex = 2.3, fill = dot_col, alpha = 0.6, height = 4) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA)) +
  labs(x = "Log2 fold cahnge", y = "-log10(adj pval)")

# plot specific genes
plots <- VlnPlot(subset(int_data_M4, Cell_Type %in% c("Tcm", "Th1")), 
                 features = c("ASNS", "OAZ1", "ODC1", "SLC3A2", "SLC7A5", "SRM"), 
                 cols = c("gray80", "tomato"),
                 pt.size = 0,
                 combine = F)

plots <- lapply(plots, function(x) {
  x + theme(plot.title = element_text(size = 10),
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 9),
            axis.title.y = element_blank()) +
    NoLegend()
})

patchwork::wrap_plots(plots, ncol = 2)











