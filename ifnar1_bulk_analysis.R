library(tidyverse)
library(magrittr)
library(edgeR)

# read in and clean data --------------------------------------------------
df <- read.table("counts.txt", header = T, row.names = 1) 

# process data ------------------------------------------------------------
groups <- factor(substr(colnames(df), 1, 2), levels = c("wt", "ko"))

y <- DGEList(counts = df, group = groups)
keep_genes <- filterByExpr(y)
y <- y[keep_genes, , keep.lib.sizes = F]
y <- calcNormFactors(y)

# normalised data ---------------------------------------------------------
norm_data <- data.frame(cpm(y))

par(mfrow = c(1, 2))
boxplot(log2(y$counts), outline = F, main = "Pre normalisation",
        col = "#84abf0", las = 2)
boxplot(norm_data, outline = F, main = "Post normalisation",
        col = "#84abf0", las = 2)


# sequencing depth --------------------------------------------------------
par(mfrow = c(1, 1))

seq_depth <- y$samples$lib.size/1e6
barplot(seq_depth, col = "#84abf0", las = 2, main = "Sequencing depth (millions) \n WT WT WT KO KO KO",
        ylab = "Read number / 1e6")

# stats testing -----------------------------------------------------------
design <- model.matrix(~0 + groups)

y <- estimateDisp(y, design)

de <- exactTest(y, pair = c("wt", "ko"))

degs <- topTags(de, n = nrow(df)) %>%
  data.frame()

# pca ---------------------------------------------------------------------
pca_df <- norm_data %>%
  + 1 %>%
  log2() %>%
  t() %>%
  prcomp()

pca_gg <- data.frame(pc1 = pca_df$x[, "PC1"],
                     pc2 = pca_df$x[, "PC2"],
                     genotype = groups)

stdev_val <- pca_df$sdev^2
stdev_val <- round(stdev_val/sum(stdev_val)*100, 1)

pca_gg %>%
  ggplot(aes(pc1, pc2, fill = genotype)) +
  ggrepel::geom_text_repel(label = rownames(pca_gg), force = 50) +
  geom_point(shape = 21, size = 3.5) +
  scale_fill_manual(values = c("gray80", "tomato")) +
  labs(x = paste0("PC1: ", stdev_val[1], "%"),
       y = paste0("PC2: ", stdev_val[2], "%")) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1)

# volcano -----------------------------------------------------------------
up_wt <- degs %>%
  filter(FDR < 0.05 & logFC < -0.6)

up_ko <- degs %>%
  filter(FDR < 0.05 & logFC > 0.6)

label_ko <- degs[rownames(degs) %in% c("Mx1", "Mx2", "Oas1a"), ]

degs %>%
  ggplot(aes(logFC, -log10(FDR))) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "black", linetype = "dashed", lwd = 0.3) +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed", lwd = 0.3) +
  geom_point(shape = 21, fill = "gray80", stroke = 0.2, alpha = 0.6) +
  geom_point(data = up_wt, shape = 21, fill = "cornflowerblue", stroke = 0.2, size = 2.3) +
  geom_point(data = up_ko, shape = 21, fill = "tomato", stroke = 0.2, size = 2.3) +
  ggrepel::geom_label_repel(data = label_ko, label = rownames(label_ko), 
                            fontface = "italic",
                            nudge_x = -1.5, label.padding = unit(0.7, "mm"),
                            label.r = 0.1, point.padding = 0.3) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1)

# highlight specific genes ------------------------------------------------
cyt_genes <- c("Il1a", "Il1b", "Ifng", "Il27",
               "Oas1a", "Mx1", "Mx2")

norm_data[rownames(norm_data) %in% cyt_genes, ] %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = wt:ko_2, names_to = "condition", values_to = "expression") %>%
  mutate(condition = substr(condition, 1, 2)) %>%
  mutate(condition = factor(toupper(condition), levels = c("WT", "KO"))) %>%
  ggplot(aes(condition, expression, fill = condition)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(shape = 21, width = 0.2, alpha = 0.8) +
  scale_fill_manual(values = c("WT" = "gray80", "KO" = "tomato")) +
  scale_x_discrete(label = c("WT", expression(italic("Ifnar1-/-")))) +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(~gene, scales = "free_y", ncol = 4) +
  labs(y = "TMM normalised expression") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "none",
        axis.title.x = element_blank(),
        strip.text = element_text(face = "italic"))















