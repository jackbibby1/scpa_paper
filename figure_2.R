
# packages ----------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(Seurat)

# each broad T cell population processed independently
# 1. naive cd4 ---------------------------------------------------------------
# read files --------------------------------------------------------------

N40 <- Read10X("C4N/0H/filtered_feature_bc_matrix/")
N412 <- Read10X("C4N/12H/filtered_feature_bc_matrix/")
N424 <- Read10X("C4N/24H/filtered_feature_bc_matrix/")


# create object and filter ------------------------------------------------
N40 <- CreateSeuratObject(N40, min.cells = 3, min.features = 800) %>% 
  PercentageFeatureSet(pattern = "MT-", col.name = "percent.mt") %>%
  subset(subset = percent.mt < 10 & nFeature_RNA < median(N40$nFeature_RNA)*2) %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
## 5219 cells filtered to 4428

N412 <- CreateSeuratObject(N412, min.cells = 3, min.features = 800) %>% 
  PercentageFeatureSet(pattern = "MT-", col.name = "percent.mt") %>%
  subset(subset = percent.mt < 10 & nFeature_RNA < median(N412$nFeature_RNA)*2) %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
## 5062 cells filtered to 4547

N424 <- CreateSeuratObject(N424, min.cells = 3, min.features = 700) %>% 
  PercentageFeatureSet(pattern = "MT-", col.name = "percent.mt") %>%
  subset(subset = percent.mt < 10 & nFeature_RNA < median(N424$nFeature_RNA)*2) %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
## 6637 cells filtered to 5919

load("CC_genes.RData")

`%notin%` <- Negate(`%in%`)
VariableFeatures(N40) <- VariableFeatures(N40)[VariableFeatures(N40) %notin% cc_genes_manual]
VariableFeatures(N412) <- VariableFeatures(N412)[VariableFeatures(N412) %notin% cc_genes_manual]
VariableFeatures(N424) <- VariableFeatures(N424)[VariableFeatures(N424) %notin% cc_genes_manual]


# add cell metadata -------------------------------------------------------
N40$Hour <- 0
N412$Hour <- 12
N424$Hour <- 24


# integrate data ----------------------------------------------------------
cell_list <- list(N40, N412, N424)
anchors <- FindIntegrationAnchors(cell_list, dims = 1:30, normalization.method = "LogNormalize")
int_data_N4 <- IntegrateData(anchors, dims = 1:30)


# scale and regress -------------------------------------------------------
s_genes <- cc.genes.updated.2019$s.genes
g2m_genes <- cc.genes.updated.2019$g2m.genes
int_data_N4 <- CellCycleScoring(int_data_N4, 
                                s.features = s_genes, g2m.features = g2m_genes, set.ident = T)
int_data_N4 <- ScaleData(int_data_N4, 
                         vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
int_data_N4 <- RunPCA(int_data_N4, verbose = F)
int_data_N4 <- RunUMAP(int_data_N4, dims = 1:30)
int_data_N4 <- FindNeighbors(int_data_N4, dims = 1:30)
int_data_N4 <- FindClusters(int_data_N4, resolution = seq(0.1, 1, 0.1))


# find cluster markers ----------------------------------------------------
DefaultAssay(int_data_N4) <- "RNA"
cluster_markers <- FindAllMarkers(int_data_N4, min.pct = 0.1, logfc.threshold = 0.25,
                                  only.pos = T)
write.csv(cluster_markers, file = "NaiveCD4_Cluster_Markers.csv", row.names = F)


# define cell types -------------------------------------------------------
cluster_ids <- c("Resting", "Intermediate", "Activated", "Treg")
names(cluster_ids) <- levels(int_data_N4)
int_data_N4 <- RenameIdents(int_data_N4, cluster_ids)
int_data_N4 <- AddMetaData(int_data_N4, int_data_N4@active.ident, "Cell_Type")

save(int_data_N4, file = "int_data_N4.RData")


# 2. memory cd4 --------------------------------------------------------------
# read files --------------------------------------------------------------
M40 <- Read10X("C4M/0H/filtered_feature_bc_matrix/")
M412 <- Read10X("C4M/12H/filtered_feature_bc_matrix/")
M424 <- Read10X("C4M/24H/filtered_feature_bc_matrix/")


# create object and filter ------------------------------------------------
M40 <- CreateSeuratObject(M40, min.cells = 3, min.features = 1200) %>% 
  PercentageFeatureSet(pattern = "MT-", col.name = "percent.mt") %>%
  subset(subset = percent.mt < 10 & nFeature_RNA < median(M40$nFeature_RNA)*2) %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
## 4304 cells filtered to 2953

M412 <- CreateSeuratObject(M412, min.cells = 3, min.features = 1000) %>% 
  PercentageFeatureSet(pattern = "MT-", col.name = "percent.mt") %>%
  subset(subset = percent.mt < 10 & nFeature_RNA < median(M412$nFeature_RNA)*2) %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
## 4974 cells filtered to 3543

M424 <- CreateSeuratObject(M424, min.cells = 3, min.features = 800) %>% 
  PercentageFeatureSet(pattern = "MT-", col.name = "percent.mt") %>%
  subset(subset = percent.mt < 10 & nFeature_RNA < median(M424$nFeature_RNA)*2) %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
## 5754 cells filtered to 5283

`%notin%` <- Negate(`%in%`)
VariableFeatures(M40) <- VariableFeatures(M40)[VariableFeatures(M40) %notin% words] # Removed 9 CC genes
VariableFeatures(M412) <- VariableFeatures(M412)[VariableFeatures(M412) %notin% words] # Removed 7 CC genes
VariableFeatures(M424) <- VariableFeatures(M424)[VariableFeatures(M424) %notin% words] # Removed 19 CC genes


# add cell metadata -------------------------------------------------------
M40$Hour <- 0
M412$Hour <- 12
M424$Hour <- 24


# integrate data ----------------------------------------------------------
cell_list <- list(M40, M412, M424)
anchors <- FindIntegrationAnchors(cell_list, dims = 1:30, normalization.method = "LogNormalize")
int_data_M4 <- IntegrateData(anchors, dims = 1:30)
rm(M40, M412, M424, anchors, cell_list)


# scale and regress -------------------------------------------------------
s_genes <- cc.genes.updated.2019$s.genes
g2m_genes <- cc.genes.updated.2019$g2m.genes
int_data_M4 <- CellCycleScoring(int_data_M4, 
                                s.features = s_genes, g2m.features = g2m_genes, set.ident = F)
int_data_M4 <- ScaleData(int_data_M4, 
                         vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
int_data_M4 <- RunPCA(int_data_M4, verbose = F)
int_data_M4 <- RunUMAP(int_data_M4, dims = 1:30)
int_data_M4 <- FindNeighbors(int_data_M4, dims = 1:30)
int_data_M4 <- FindClusters(int_data_M4, resolution = seq(0.1, 1, 0.1))


# find cluster markers ----------------------------------------------------
DefaultAssay(int_data_M4) <- "RNA"
cluster_markers <- FindAllMarkers(int_data_M4, min.pct = 0.1, logfc.threshold = 0.25,
                                  only.pos = T)
write.csv(cluster_markers, file = "MemoryCD4_Cluster_Markers.csv", row.names = F)


# define cell types -------------------------------------------------------
cluster_ids <- c("Tcm", "Teff", "Th1", "Th2", "Treg", "Proliferating")
names(cluster_ids) <- levels(int_data_M4)
int_data_M4 <- RenameIdents(int_data_M4, cluster_ids)
int_data_M4 <- AddMetaData(int_data_M4, int_data_M4@active.ident, "Cell_Type")

save(int_data_M4, file = "int_data_M4.RData")


# 3. naive cd8 ---------------------------------------------------------------
# read files --------------------------------------------------------------
N80 <- Read10X("C8N/0H/filtered_feature_bc_matrix/")
N812 <- Read10X("C8N/12H/filtered_feature_bc_matrix/")
N824 <- Read10X("C8N/24H/filtered_feature_bc_matrix/")


# create object and filter ------------------------------------------------
N80 <- CreateSeuratObject(N80, min.cells = 3, min.features = 1300) %>% 
  PercentageFeatureSet(pattern = "MT-", col.name = "percent.mt") %>%
  subset(subset = percent.mt < 10 & nFeature_RNA < median(N80$nFeature_RNA)*2) %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
## 2794 cells filtered to 1048

N812 <- CreateSeuratObject(N812, min.cells = 3, min.features = 1300) %>% 
  PercentageFeatureSet(pattern = "MT-", col.name = "percent.mt") %>%
  subset(subset = percent.mt < 10 & nFeature_RNA < median(N812$nFeature_RNA)*2) %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
## 3601 cells filtered to 2066

N824 <- CreateSeuratObject(N824, min.cells = 3, min.features = 950) %>% 
  PercentageFeatureSet(pattern = "MT-", col.name = "percent.mt") %>%
  subset(subset = percent.mt < 10 & nFeature_RNA < median(N824$nFeature_RNA)*2) %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
## 5932 cells filtered to 3927

load("CC_genes.RData")

`%notin%` <- Negate(`%in%`)
VariableFeatures(N80) <- VariableFeatures(N80)[VariableFeatures(N80) %notin% words] # Removed 7 CC genes
VariableFeatures(N812) <- VariableFeatures(N812)[VariableFeatures(N812) %notin% words] # Removed 11 CC genes
VariableFeatures(N824) <- VariableFeatures(N824)[VariableFeatures(N824) %notin% words] # Removed 20 CC genes


# add cell metadata -------------------------------------------------------
N80$Hour <- 0
N812$Hour <- 12
N824$Hour <- 24


# integrate data ----------------------------------------------------------
cell_list <- list(N80, N812, N824)
anchors <- FindIntegrationAnchors(cell_list, dims = 1:30, normalization.method = "LogNormalize")
int_data_N8 <- IntegrateData(anchors, dims = 1:30)
rm(N80, N812, N824, anchors, cell_list)


# scale and regress -------------------------------------------------------
s_genes <- cc.genes.updated.2019$s.genes
g2m_genes <- cc.genes.updated.2019$g2m.genes
int_data_N8 <- CellCycleScoring(int_data_N8, 
                                s.features = s_genes, g2m.features = g2m_genes, set.ident = T)
int_data_N8 <- ScaleData(int_data_N8, 
                         vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
int_data_N8 <- RunPCA(int_data_N8, verbose = F)
int_data_N8 <- RunUMAP(int_data_N8, dims = 1:30)
int_data_N8 <- FindNeighbors(int_data_N8, dims = 1:30)
int_data_N8 <- FindClusters(int_data_N8, resolution = seq(0.1, 1, 0.1))


# find cluster markers ----------------------------------------------------
DefaultAssay(int_data_N8) <- "RNA"
cluster_markers <- FindAllMarkers(int_data_N8, min.pct = 0.1, logfc.threshold = 0.25,
                                  only.pos = T)
write.csv(cluster_markers, file = "NaiveCD8_Cluster_Markers.csv", row.names = F)


# defince cell types ------------------------------------------------------
cluster_ids <- c("Resting", "Intermediate", "Activated", "TEMRA", "Undefined1", "Undefined2")
names(cluster_ids) <- levels(int_data_N8)
int_data_N8 <- RenameIdents(int_data_N8, cluster_ids)
int_data_N8 <- AddMetaData(int_data_N8, int_data_N8@active.ident, "Cell_Type")

save(int_data_N8, file = "int_data_N8.RData")


# 4. memory cd8 --------------------------------------------------------------
# read files --------------------------------------------------------------
M80 <- Read10X("C8M/0H/filtered_feature_bc_matrix/")
M812 <- Read10X("C8M/12H/filtered_feature_bc_matrix/")
M824 <- Read10X("C8M/24H/filtered_feature_bc_matrix/")

## Create object and filter
M80 <- CreateSeuratObject(M80, min.cells = 3, min.features = 1300) %>% 
  PercentageFeatureSet(pattern = "MT-", col.name = "percent.mt") %>%
  subset(subset = percent.mt < 10 & nFeature_RNA < median(M80$nFeature_RNA)*2) %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
## 2317 cells filtered to 737

M812 <- CreateSeuratObject(M812, min.cells = 3, min.features = 1300) %>% 
  PercentageFeatureSet(pattern = "MT-", col.name = "percent.mt") %>%
  subset(M812, subset = percent.mt < 10 & nFeature_RNA < median(M812$nFeature_RNA)*2) %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
## 3246 cells filtered to 2189

M824 <- CreateSeuratObject(M824, min.cells = 3, min.features = 800) %>% 
  PercentageFeatureSet(pattern = "MT-", col.name = "percent.mt") %>%
  subset(M824, subset = percent.mt < 10 & nFeature_RNA < median(M824$nFeature_RNA)*2) %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
## 5820 cells filtered to 4846

`%notin%` <- Negate(`%in%`)
VariableFeatures(M80) <- VariableFeatures(M80)[VariableFeatures(M80) %notin% words] # Removed 7 CC genes
VariableFeatures(M812) <- VariableFeatures(M812)[VariableFeatures(M812) %notin% words] # Removed 15 CC genes
VariableFeatures(M824) <- VariableFeatures(M824)[VariableFeatures(M824) %notin% words] # Removed 20 CC genes


# add cell metadata -------------------------------------------------------
M80$Hour <- 0
M812$Hour <- 12
M824$Hour <- 24


# integrate ---------------------------------------------------------------
cell_list <- list(M80, M812, M824)
anchors <- FindIntegrationAnchors(cell_list, dims = 1:30, normalization.method = "LogNormalize")
int_data_M8 <- IntegrateData(anchors, dims = 1:30)
rm(M80, M812, M824, anchors, cell_list)


# scale and regress -------------------------------------------------------
DefaultAssay(int_data_M8) <- "integrated"
int_data_M8 <- CellCycleScoring(int_data_M8, 
                                s.features = s_genes, g2m.features = g2m_genes, set.ident = T)
int_data_M8 <- ScaleData(int_data_M8, 
                         vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
int_data_M8 <- RunPCA(int_data_M8, verbose = F)
int_data_M8 <- RunUMAP(int_data_M8, dims = 1:30)
int_data_M8 <- FindNeighbors(int_data_M8, dims = 1:30)
int_data_M8 <- FindClusters(int_data_M8, resolution = seq(0.1, 1, 0.1))


# find cluster markers ----------------------------------------------------
DefaultAssay(int_data_M8) <- "RNA"
cluster_markers <- FindAllMarkers(int_data_M8, min.pct = 0.1, logfc.threshold = 0.25,
                                  only.pos = T)
write.csv(cluster_markers, file = "MemoryCD8_Cluster_Markers.csv", row.names = F)


# define cell types -------------------------------------------------------
cluster_ids <- c("Tcm", "Tc1", "Tc2", "Tem1", "Tem2", "GZM+", "Proliferating")
names(cluster_ids) <- levels(int_data_M8)
int_data_M8 <- RenameIdents(int_data_M8, cluster_ids)

save(int_data_M8, file = "int_data_M8.RData")


# 5. umaps -------------------------------------------------------------------
# umap extraction and theme -----------------------------------------------
extract_umap <- function(seu_obj) {
  
  embed <- Embeddings(seu_obj, "umap") %>%
    data.frame() %>%
    rownames_to_column("cell")
  cell_type <- seu_obj$Cell_Type %>%
    data.frame() %>% 
    set_colnames("cell_type") %>%
    rownames_to_column("cell")
  hour <- seu_obj$Hour %>%
    data.frame() %>%
    set_colnames("hour") %>%
    rownames_to_column("cell")
  
  umap_embeds <- Reduce(full_join, list(embed, cell_type, hour))
  return(umap_embeds)
  
}

umap_theme <- theme(panel.background = element_blank(),
                    panel.border = element_blank(),
                    strip.background = element_blank(),
                    strip.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title = element_blank(),
                    axis.text = element_blank(),
                    aspect.ratio = 1)

# naive cd4
n4_umap <- extract_umap(int_data_N4)

p1 <- ggplot(n4_umap, aes(UMAP_1, UMAP_2, fill = cell_type)) +
  geom_point(shape = 21, cex = 0.65, color = 'black', stroke = 0.02, alpha = 0.8) +
  facet_wrap(~hour, ncol = 3) +
  scale_fill_manual(values = c("Resting" = "cornflowerblue",
                               "Intermediate" = "gray80",
                               "Activated" = "tomato",
                               "Treg" = "#1e9c2f")) +
  umap_theme

# memory cd4
m4_umap <- extract_umap(int_data_M4)

p2 <- ggplot(m4_umap, aes(UMAP_1, UMAP_2, fill = cell_type)) +
  geom_point(shape = 21, cex = 0.65, color = 'black', stroke = 0.02, alpha = 0.8) +
  facet_wrap(~hour, ncol = 3) +
  scale_fill_manual(values = c("Tcm" = "mediumpurple1",
                               "Teff" = "gray70",
                               "Th1" = "tomato",
                               "Th2" = "orange",
                               "Treg" = "#1e9c2f",
                               "Proliferating" = "cornflowerblue")) +
  umap_theme

# naive cd8
n8_umap <- extract_umap(int_data_N8)

p3 <- ggplot(n8_umap, aes(UMAP_1, UMAP_2, fill = cell_type)) +
  geom_point(shape = 21, cex = 0.65, color = 'black', stroke = 0.02, alpha = 0.8) +
  facet_wrap(~hour, ncol = 3) +
  scale_fill_manual(values = c("Resting" = "cornflowerblue",
                               "Intermediate" = "gray80",
                               "Activated" = "tomato",
                               "TEMRA" = "mediumorchid2",
                               "Undefined1" = "#1e9c2f",
                               "Undefined2" = "goldenrod3")) +
  umap_theme

# memory cd8
m8_umap <- extract_umap(int_data_M8)

p4 <- ggplot(m8_umap, aes(UMAP_1, UMAP_2, fill = cell_type)) +
  geom_point(shape = 21, cex = 0.65, color = 'black', stroke = 0.02, alpha = 0.8) +
  facet_wrap(~hour, ncol = 3) +
  scale_fill_manual(values = c("Tcm" = "#477afc",
                               "Tc1" = "#ff7f1f",
                               "Tc2" = "mediumorchid2",
                               "Tem1" = "#c2c2c2",
                               "Tem2" = "steelblue2",
                               "GZM+" = "#1e9c2f",
                               "Proliferating" = "goldenrod3")) +
  umap_theme


plots <- list(p1, p2, p3, p4)
plots <- lapply(plots, function(x){
  x + theme(axis.line = element_blank(), 
            axis.ticks = element_blank(),
            axis.text = element_blank(), 
            axis.title.x = element_blank(),
            axis.title.y = element_blank(), 
            strip.text = element_blank(),
            legend.text = element_text(size = 11),
            legend.title = element_blank(),
            legend.key = element_blank(),
            legend.key.height = unit(0.4, "cm"), 
            plot.margin = margin(1, 0.75, 0.45, 0, "cm")) +
    guides(fill = guide_legend(override.aes = list(size = 3)))
})

wrap_plots(plots, ncol = 2)
ggsave("test_umap.png", device = "png", 
       width = 13, height = 4, units = "in")
dev.off()

# 6. cell proportions --------------------------------------------------------
# proportion extraction function ------------------------------------------

proportions <- function(seu_obj, colors) {
  df <- as.data.frame(table(subset(seu_obj, Hour == 0)$Cell_Type)) %>%
    set_colnames(c("Cell", "f0")) %>%
    mutate(f12 = as.data.frame(table(subset(seu_obj, Hour == 12)$Cell_Type))[,2]) %>%
    mutate(f24 = as.data.frame(table(subset(seu_obj, Hour == 24)$Cell_Type))[,2]) %>%
    mutate(`0` = prop.table(f0)*100) %>%
    mutate(`12` = prop.table(f12)*100) %>%
    mutate(`24` = prop.table(f24)*100) %>%
    select(Cell, `0`, `12`, `24`) %>%
    mutate(color = colors) %>%
    pivot_longer(., cols = `0`:`24`, names_to = "Time", values_to = "Proportion")
}


# get proportions ---------------------------------------------------------
n4_cell <- proportions(int_data_N4, colors = c("cornflowerblue", "gray80", "tomato", "#1e9c2f"))
m4_cell <- proportions(int_data_M4, colors = c("mediumpurple1", "gray70", "tomato", "orange",
                                               "#1e9c2f", "cornflowerblue"))
n8_cell <- proportions(int_data_N8, colors = c("cornflowerblue", "gray80", "tomato", "mediumorchid2",
                                               "#1e9c2f", "goldenrod3"))
m8_cell <- proportions(int_data_N8, colors = c("#477afc", "#ff7f1f", "mediumorchid2", "#c2c2c2",
                                               "steelblue2", "#1e9c2f", "goldenrod3"))


p1 <- ggplot(n4_cell, aes(x = Time, y = Proportion)) +
  geom_col(fill = n4_cell$color, col = 'white', size = 0.2) +
  labs(title = "Naive CD4") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 11),
        plot.title = element_text(size = 12, hjust = 0.5))

p2 <- ggplot(m4_cell, aes(x = Time, y = Proportion)) +
  geom_col(fill = m4_cell$color, col = 'white', size = 0.2) +
  labs(title = "Memory CD4") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 11),
        plot.title = element_text(size = 12, hjust = 0.5))

p3 <- ggplot(n8_cell, aes(x = Time, y = Proportion)) +
  geom_col(fill = n8_cell$color, col = 'white', size = 0.2) +
  labs(title = "Naive CD8") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 11),
        plot.title = element_text(size = 12, hjust = 0.5))

p4 <- ggplot(m8_cell, aes(x = Time, y = Proportion)) +
  geom_col(fill = m8_cell$color, col = 'white', size = 0.2) +
  labs(title = "Memory CD8") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 11),
        plot.title = element_text(size = 12, hjust = 0.5))

## Plot together
plot <- wrap_plots(p1, p2, p3, p4, ncol = 2)
ggsave(plot = plot, filename = "Cluster_proportions.png", 
       width = 3.6, height = 6.5, units = "in", dpi = 600)

# 7. dot plots ---------------------------------------------------------------
# function for top markers ------------------------------------------------

top_markers <- function(filepath, cell_type) {
  df <- read.csv(filepath) %>%
    filter(cluster == cell_type) %>% 
    top_n(7, wt = avg_logFC) %>% 
    pull(gene) %>% 
    as.data.frame() %>%
    mutate(Cluster = cell_type) %>% 
    set_colnames(c("Gene", "Cluster"))
}

# cd4 markers
memory_cd4_markers <- list()
for (i in levels(int_data_M4)) {
  memory_cd4_markers[[i]] <- top_markers("MemoryCD4_Cluster_Markers.csv", i)
}
memory_cd4_markers <- bind_rows(memory_cd4_markers)

naive_cd4_markers <- list()
for (i in levels(int_data_N4)) {
  naive_cd4_markers[[i]] <- top_markers("NaiveCD4_Cluster_Markers.csv", i)
}
naive_cd4_markers <- bind_rows(naive_cd4_markers)

mks_cd4 <- c(naive_cd4_markers$Gene, memory_cd4_markers$Gene) %>% unique()

int_data_N4 <- RenameIdents(int_data_N4, "Treg" = "nTreg")
int_data_M4 <- RenameIdents(int_data_M4, "Treg" = "mTreg")

cd4 <- merge(int_data_N4, int_data_M4, merge.data = T)

levels(cd4) <- c("Resting", "Intermediate", "Activated", "nTreg", "mTreg", "Tcm", "Teff", "Th1", "Th2", "Proliferating") %>%
  rev()

p1 <- DotPlot(cd4, features = mks_cd4, dot.scale = 4) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 12))

# cd8 markers
memory_cd8_markers <- list()
for (i in levels(int_data_M8)) {
  memory_cd8_markers[[i]] <- top_markers("MemoryCD8_Cluster_Markers.csv", i)
}
memory_cd8_markers <- bind_rows(memory_cd8_markers)

naive_cd8_markers <- list()
for (i in levels(int_data_N8)) {
  naive_cd8_markers[[i]] <- top_markers("NaiveCD8_Cluster_Markers.csv", i)
}
naive_cd8_markers <- bind_rows(naive_cd8_markers)

mks_cd8 <- c(naive_cd8_markers$Gene, memory_cd8_markers$Gene) %>% unique()

cd8 <- merge(int_data_N8, int_data_M8, merge.data = T)

levels(cd8) <- c("Resting", "Intermediate", "Activated", "TEMRA", "Undefined1", "Undefined2",
                 "Tcm", "Tc1", "Tc2", "Tem1", "Tem2", "GZM+", "Proliferating") %>%
  rev()

p2 <- DotPlot(cd8, features = mks_cd8, dot.scale = 4) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 12))


patchwork::wrap_plots(p1, p2, ncol = 1)
ggsave(file = "dot_plots.png",
       device = "png", units = "cm", dpi = 600,
       height = 19, width = 44)
dev.off()










