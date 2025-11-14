library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
# here we don't use the filtered data after QC for demonstration
seurat_obj_A1 <- readRDS("./data/seurat_obj_A1.rds")
seurat_obj_B1 <- readRDS("./data/seurat_obj_B1.rds")
# sctransform respectively
seurat_obj_A1 <- SCTransform(seurat_obj_A1, assay = "Spatial", verbose = FALSE)
seurat_obj_B1 <- SCTransform(seurat_obj_B1, assay = "Spatial", verbose = FALSE)

# Integrate the two spatial datasets
spatial_list <- list(seurat_obj_A1, seurat_obj_B1)
features <- SelectIntegrationFeatures(object.list = spatial_list, nfeatures = 3000)
# Preprocessing for SCT integration
# reference: https://github.com/satijalab/seurat/issues/8216
spatial_list[[1]][["RNA"]] <-spatial_list[[1]][["Spatial"]]
spatial_list[[2]][["RNA"]] <-spatial_list[[2]][["Spatial"]]
# add sample information
spatial_list[[1]]$sample <- "A1_colon_d0"
spatial_list[[2]]$sample <- "B1_colon_d14"
spatial_list <- PrepSCTIntegration(object.list = spatial_list, anchor.features = features)

# Using anchors to integrate data
anchors <- FindIntegrationAnchors(
  object.list = spatial_list,
  normalization.method = "SCT",
  anchor.features = features
)

# create the integrated data
seurat_combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
# Run the standard workflow for visualization and clustering
seurat_combined <- RunPCA(seurat_combined, verbose = FALSE)
# Choose the best number of pc to use
pct <- seurat_combined[["pca"]]@stdev / sum(seurat_combined[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1], sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1)

ElbowPlot(seurat_combined, ndims = pc.use+5)$data %>% ggplot() +
  geom_point(aes(x = dims, y = stdev)) +
  geom_vline(xintercept = pc.use, color = "darkred") +
  theme_bw() + labs(title = "Elbow plot: quantitative approach")

seurat_combined <- RunUMAP(seurat_combined, dims = 1:pc.use)
seurat_combined <- FindNeighbors(seurat_combined, dims = 1:pc.use)
seurat_combined <- FindClusters(seurat_combined, resolution = 0.5)
DimPlot(seurat_combined, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()
p <- DimPlot(seurat_combined, group.by = "sample", label = FALSE)
ggsave("./results/umap_integrated_by_sample.tiff", p, width = 8, height = 6)
p <- DimPlot(seurat_combined, group.by = "seurat_clusters", label = TRUE)
ggsave("./results/umap_integrated_by_cluster.tiff", p, width = 8, height = 6)
DimPlot(seurat_combined, group.by = "seurat_clusters", label = TRUE)
# Save the integrated data
saveRDS(seurat_combined, file = "./data/seurat_integrated_spatial.rds")
# visualization the clustering result on the spatial location
p <- Seurat::SpatialDimPlot(seurat_combined,
               pt.size.factor = 3,
               image.scale = "hires", # important when using high resolution image
               group.by = "seurat_clusters"
)
ggsave("./results/spatial_cluster_integrated.tiff", p, width = 8, height = 6)
# Find markers for each cluster
cluster_markers <- FindAllMarkers(seurat_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster_markers, file = "./results/integrated_cluster_markers.csv")
top10 <- cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
p <- DoHeatmap(seurat_combined, features = top10$gene) + NoLegend()
ggsave("./results/integrated_cluster_markers_heatmap.tiff", p, width = 18, height = 18)

# DEG and Vocalno Plot between A1_colon_d0 and B1_colon_d14
de_genes <- FindMarkers(
  seurat_combined,
  ident.1 = "B1_colon_d14",
  ident.2 = "A1_colon_d0", #control group
  group.by = "sample",
  min.pct = 0.1,
  logfc.threshold = 0.25
)
write.csv(de_genes, file = "./results/integrated_de_genes_A1_vs_B1.csv")
de_genes$sig <- ifelse(de_genes$avg_log2FC>0, "Up", "Down")
de_genes$sig <- ifelse(de_genes$p_val_adj < 0.05 & abs(de_genes$avg_log2FC) >= 1, de_genes$sig , "no")
ggplot(de_genes, aes(x = avg_log2FC, y = -log10(p_val_adj), color=sig)) +
  geom_point(alpha = 1) +
  xlab("Average Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value") +
  ggtitle("Volcano Plot: A1_colon_d0 vs B1_colon_d14") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  scale_color_manual(values = list(
    "Up" = "red",
    "Down" = "blue",
    "no" = "grey"
  )) +
  ggprism::theme_prism()
ggsave("./results/volcano_plot_A1_vs_B1.tiff", width = 6, height = 3, dpi = 300)

#optional: visualize some gene expression
p <- SpatialFeaturePlot(seurat_combined, rownames(de_genes)[1:2], image.scale = "hires", pt.size.factor = 2)
ggsave("./results/spatial_feature_de_genes.tiff", p, width = 10, height = 6, dpi=300)


# Spatially Variable Features Detection -----------------------------
# Identify spatially variable features
# attention here to run FindSpatiallyVariableFeatures on each sample separately
seurat_obj_A1 <- FindSpatiallyVariableFeatures(
  seurat_obj_A1,
  assay = "SCT",
  features = VariableFeatures(seurat_combined)[1:1000],
  selection.method = "moransi", # faster than markvariogram, install.packages('Rfast2') for faster computation
  verbose = TRUE
)
top_sv_genes <- head(SpatiallyVariableFeatures(seurat_obj_A1, selection.method = "moransi"), 2)
#  visualization of spatially variable features
p <- SpatialFeaturePlot(seurat_obj_A1, features = top_sv_genes, ncol = 2, image.scale = "hires", pt.size.factor = 2.5)
ggsave("./results/spatially_variable_genes_A1.tiff", p, width = 10, height = 6, dpi=300)

seurat_obj_B1 <- FindSpatiallyVariableFeatures(
  seurat_obj_B1,
  assay = "SCT",
  features = VariableFeatures(seurat_combined)[1:1000],
  selection.method = "moransi", # faster than markvariogram, install.packages('Rfast2') for faster computation
  verbose = TRUE
)
top_sv_genes <- head(SpatiallyVariableFeatures(seurat_obj_B1, selection.method = "moransi"), 2)
#  visualization of spatially variable features
p <- SpatialFeaturePlot(seurat_obj_B1, features = top_sv_genes, ncol = 2, image.scale = "hires", pt.size.factor = 2)
ggsave("./results/spatially_variable_genes_B1.tiff", p, width = 10, height = 6, dpi=300)




