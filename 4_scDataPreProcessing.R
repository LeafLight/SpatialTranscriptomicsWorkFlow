library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
sc_seurat <- Read10X("./data/GSE264408_RAW/controlExample")
  sc_seurat <- CreateSeuratObject(
    counts = sc_seurat,
    project = "Control",
    min.cells = 3,
    min.features = 200
  )
sc_seurat[["percent.mt"]] <- PercentageFeatureSet(sc_seurat, pattern = "^mt-")
sc_seurat[["percent.mt"]]
qc_violin_before <- VlnPlot(sc_seurat,
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                     group.by = "orig.ident",
                    ncol = 4, pt.size = 0)

sc_seurat <- subset(sc_seurat, subset =  nFeature_RNA > 200 & nFeature_RNA < 7500 &
                    nCount_RNA > 500 & nCount_RNA < 50000 &
                    percent.mt < 25
                    )
qc_violin_after <- VlnPlot(sc_seurat,
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                     group.by = "orig.ident",
                    ncol = 4, pt.size = 0)
qc_violin <- qc_violin_before / qc_violin_after + plot_annotation(title = "scRNA-seq QC before and after filtering")
ggsave("./results/sc_qc_violin.tiff", qc_violin, width = 8, height = 8, dpi = 300)
p1 <- FeatureScatter(sc_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt") +NoLegend()+
  geom_hline(yintercept = 20, linetype = "dashed", color = "red")
p2 <- FeatureScatter(sc_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+NoLegend()
qc_scatter <- p1 + p2
ggsave("./results/sc_qc_scatter.tiff", qc_scatter, width = 8, height = 4, dpi = 300)

sc_seurat <- SCTransform(sc_seurat, vars.to.regress = "percent.mt", verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000)
sc_seurat <- RunPCA(sc_seurat)
# Choose the best number of pc to use
pct <- sc_seurat[["pca"]]@stdev / sum(sc_seurat[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1], sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1)

ElbowPlot(sc_seurat, ndims = pc.use+5)$data %>% ggplot() +
  geom_point(aes(x = dims, y = stdev)) +
  geom_vline(xintercept = pc.use, color = "darkred") +
  theme_bw() + labs(title = "Elbow plot: quantitative approach")
ggsave("./results/sc_elbow_plot.tiff", width = 6, height = 4, dpi = 300)
sc_seurat <- RunUMAP(sc_seurat, dims = 1:pc.use)
sc_seurat <- FindNeighbors(sc_seurat, dims = 1:pc.use)
sc_seurat <- FindClusters(sc_seurat, resolution = 0.5)
DimPlot(sc_seurat, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()
p <- DimPlot(sc_seurat, label = TRUE, pt.size = 1)
ggsave("./results/sc_umap_by_cluster.tiff", p, width = 8, height = 6)
p <- DimPlot(sc_seurat, group.by = "orig.ident", label = FALSE, pt.size = 1)
ggsave("./results/sc_umap_by_sample.tiff", p, width = 8, height = 6)
# Save the processed scRNA-seq data
saveRDS(sc_seurat, file = "./data/sc_seurat_processed.rds")


intestinal_markers <- list(
  "Epithelial" = c("Epcam", "Cdh1"),
  "Enterocyte" = c("Fabp1", "Fabp2", "Apoa1", "Apoa4", "Si", "Alpi"),
  "Goblet" = c("Muc2", "Tff3", "Spink4", "Fcgpb", "Clca1"),
  "Stem" = c("Lgr5", "Olfm4", "Ascl2", "Smoc2", "Axin2"),
  "Stromal" = c("Vim", "Col1a1", "Col3a1", "Pdgfra", "Pdgfrb", "Acta2"),
  "Endothelial" = c("Pecam1", "Cdh5", "Vwf", "Cldn5", "Erg"),
  "Immune" = c("Ptprc", "Cd3e", "Cd4", "Cd8a", "Cd79a", "Ms4a1", "Lyz2", "Cd68", "Itgax"),
  "Fibroblast" = c("Pdgfra", "Pdgfrb", "Col1a1", "Col3a1", "Des", "Acta2"),
  "Myofibroblast" = c("Acta2", "Tagln", "Myh11", "Des"),
  "Pericyte" = c("Cspg4", "Rgs5", "Abcc9", "Kcnj8")
)



CalculateCellTypeScores <- function(seurat_obj, marker_list) {
  for (cell_type in names(marker_list)) {
    markers <- marker_list[[cell_type]]
    available_markers <- markers[markers %in% rownames(seurat_obj)]

    if (length(available_markers) >= 2) {
      cat("Calculating score for", cell_type, "using", length(available_markers), "markers\n")
      seurat_obj <- AddModuleScore(seurat_obj,
                                  features = list(available_markers),
                                  name = cell_type,
                                  ctrl = min(50, length(available_markers)))
    }
  }
  return(seurat_obj)
}


sc_seurat <- CalculateCellTypeScores(sc_seurat, intestinal_markers)



score_plots <- list()
for (cell_type in names(intestinal_markers)) {
  score_col <- paste0(cell_type, "1")
  if (score_col %in% colnames(sc_seurat@meta.data)) {
    p <- FeaturePlot(sc_seurat, features = score_col, order = TRUE) +
      scale_color_viridis_c() +
      ggtitle(paste(cell_type, "Score"))
    score_plots[[cell_type]] <- p
  }
}


if (length(score_plots) > 0) {
  score_grid <- wrap_plots(score_plots, ncol = 3)
  ggsave("./results/scRNA_celltype_scores.tiff", score_grid, width = 15, height = 15, dpi = 300)
}

names(intestinal_markers)


score_columns <- grep("1$", colnames(sc_seurat@meta.data), value = TRUE)
DotPlot(sc_seurat, features = score_columns, group.by = "seurat_clusters") +
  RotatedAxis() +
  scale_x_discrete(labels = names(intestinal_markers)) +
  ggtitle("Cell Type Scores DotPlot")

# 8 Enterocytes,
# 9 Intestinal Stem Cells,
# 15 Goblet Cells,
# 10, 7 Endothelial Cells,
# 1, 3, 4, 6, 11, 12, 13 Stromal Cells,
# 0, 2, 5, 16, 17 Immune Cells,

sc_meta <- read.csv("./data/GSE264408_RAW/GSE264408_metadata.csv", row.names = 1)
sc_seurat$celltype_major <- sc_meta[paste0("HC1-", colnames(sc_seurat)), 'celltype_major']
DimPlot(sc_seurat, group.by = "celltype_major", label = TRUE, pt.size = 1) + NoLegend()
ggsave("./results/sc_annotated_by_major_celltype.tiff", width = 8, height = 6, dpi = 300)
saveRDS(sc_seurat,"./data/sc_seurat_annotated.rds")
sc_seurat <- readRDS("./data/sc_seurat_annotated.rds")
celltype_prop <- table(sc_seurat$celltype_major) %>%
  as.data.frame() %>%
  arrange(desc(Freq))

p_bar <- ggplot(celltype_prop, aes(x = reorder(Var1, -Freq), y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_label(aes(x = 2, y = 55, label = "Downsampling threshold (n=50)"),
          color = "red", fill = "white", label.size = 0.8, label.padding = unit(0.4, "lines"),
          hjust = 0, vjust = 0)+
  theme_minimal() +
  labs(x = "Cell Type", y = "Number of Cells", fill = "Cell Type",
       title = "scRNA-seq Reference Cell Type Distribution") +
  ggprism::theme_prism() +
  theme(axis.text.x = element_text(angle = 45),
        legend.position = "none")
ggsave("./results/sc_celltype_distribution_barplot.tiff", p_bar, width = 10, height = 8, dpi = 300)