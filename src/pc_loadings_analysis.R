library(reshape2)
library(edgeR)
library(limma)
library(gtools)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(magrittr)
library(openxlsx)

out_dir <- "results/"

dir.create(out_dir, showWarnings = F)

# pca plot with both cells

meta <- readRDS("data/fbx_oxp_count_metadata.RDS")
counts <- readRDS("data/fbx_oxp_gene_counts.RDS")

counts_dge <- DGEList(counts, group = meta$treatment)

keep <- filterByExpr(counts_dge, min.count=40)

counts_dge <- counts_dge[keep,keep.lib.sizes=FALSE]
counts_dge <- calcNormFactors(counts_dge)

counts_y <- voom(counts_dge, plot=F)

counts_pca <- prcomp(t(counts_y$E), scale. = T)


# make the plot
counts_pca_x <- counts_pca$x
counts_pca_x <- merge(counts_pca_x,
                      meta, by="row.names")

ggplot(counts_pca_x, aes(x=PC1,
                              y=PC2,
                              color=treatment,
                              shape=cell_type)) + 
  geom_point(size=4) +
  scale_color_manual(values=c("green4","purple","black")) +
  theme_bw() +
  theme(legend.text = element_text(size=15),
        axis.text = element_text(size=15),
        axis.title = element_text(size=15),
        legend.direction = "vertical",
        legend.position = "top") +
  labs(x="Principal Component 1 (PC1)", 
       y="Principal Component 2 (PC2)",  title=NULL, color=NULL, shape=NULL)


ggsave(paste0(out_dir, "fbx_obx.both_cells.pca.png"),
       width=6, height=5)

# heatmaps

for (cell in unique(meta$cell_type)) {
  
  subset_meta <- meta[meta$cell_type == cell,]
  subset_counts <- counts[,rownames(subset_meta)]
  
  counts_dge <- DGEList(subset_counts, group = subset_meta$treatment)
  
  keep <- filterByExpr(counts_dge, min.count=40)
  
  counts_dge <- counts_dge[keep,keep.lib.sizes=FALSE]
  counts_dge <- calcNormFactors(counts_dge)
  
  counts_y <- voom(counts_dge, plot=F)
  
  counts_pca <- prcomp(t(counts_y$E), scale. = T)
  
  # pc loadings analysis
  counts_rot <- as.data.frame(counts_pca$rotation)
  
  top_PC1 <- counts_rot[order(abs(counts_rot$PC1), decreasing=T),]$PC1[1:250]
  names(top_PC1) <- rownames(counts_rot[order(abs(counts_rot$PC1), decreasing=T),])[1:250]
  
  genes_pc1 <- counts_y$E[names(top_PC1),]
  
  # scale em and heatmap em
  genes_pc1_s <- t(scale(t(genes_pc1)))
  
  plot_colors <- c("green4","purple","black")
  names(plot_colors) <- levels(subset_meta$treatment)
  
  ca <- columnAnnotation(TREATMENT=subset_meta$treatment,
                         col=list("TREATMENT"=plot_colors))
  
  genes_pc1_s_dist <- dist(genes_pc1_s, method = "euclidean")
  genes_pc1_s_clust <- hclust(genes_pc1_s_dist)
  
  genes_pc1_s_up_down <- cutree(genes_pc1_s_clust, 2)
  
  png(paste0(out_dir, "fbx_obx.", cell,".cpm_z_score.heatmap.png"), width=8, height=6, units="in", res=300)
  print(Heatmap(genes_pc1_s,
                col=colorRamp2(c(-2,0,2), c("blue", "white", "red")),
                cluster_columns = T,
                cluster_rows = T,
                name = "log2CPM Z-Score",
                split = genes_pc1_s_up_down[rownames(genes_pc1_s)],
                top_annotation = ca,
                show_row_names = F))
  dev.off()
}


