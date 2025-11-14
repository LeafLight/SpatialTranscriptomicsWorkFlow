library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
seurat_obj_A1 <- readRDS("./data/seurat_obj_A1.rds")
seurat_obj_B1 <- readRDS("./data/seurat_obj_B1.rds")
# Process QC data
seurat_obj_A1 <- PercentageFeatureSet(seurat_obj_A1, pattern = "^mt-", col.name = "percent.mt") # Mt MT
seurat_obj_B1 <- PercentageFeatureSet(seurat_obj_B1, pattern = "^mt-", col.name = "percent.mt")

k_format <- function(x) {
  case_when(
    x >= 1e6 ~ paste0(round(x/1e6, 1), "M"),
    x >= 1e3 ~ paste0(round(x/1e3, 1), "k"),
    TRUE ~ as.character(x)
  )
}
# Visulization of QC results
p1 <- VlnPlot(seurat_obj_A1, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"), ncol = 3, pt.size = 0)
p2 <- SpatialFeaturePlot(seurat_obj_A1, features = "nCount_Spatial", pt.size.factor = 3, image.scale = "hires")+
  scale_fill_distiller(
    palette = "Spectral",
    labels = k_format  # using comma format for large numbers
    )
p3 <- SpatialFeaturePlot(seurat_obj_A1, features = "percent.mt", pt.size.factor = 3, image.scale = "hires")
ggsave("./results/QC_A1.tiff", p1/(p2|p3), width = 8, height = 12, dpi=300)


p1 <- VlnPlot(seurat_obj_B1, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"), ncol = 3, pt.size = 0)
p2 <- SpatialFeaturePlot(seurat_obj_B1, features = "nCount_Spatial", pt.size.factor = 2, image.scale = "hires")+
  scale_fill_distiller(
    palette = "Spectral",
    labels = k_format  # using comma format for large numbers
    )
p3 <- SpatialFeaturePlot(seurat_obj_B1, features = "percent.mt", pt.size.factor = 2, image.scale = "hires")
ggsave("./results/QC_B1.tiff", p1/(p2|p3), width = 8, height = 12, dpi=300)


# optional: filtering spots by nFeature and nCount
sum(seurat_obj_A1$nFeature_Spatial < 200)
sum(seurat_obj_B1$nFeature_Spatial < 200)
sum(seurat_obj_A1$nCount_Spatial < 500)
sum(seurat_obj_B1$nCount_Spatial < 500)
# Filter the data based on QC metrics
seurat_obj_A1 <- subset(seurat_obj_A1, subset = nFeature_Spatial > 200 & nCount_Spatial > 500)
seurat_obj_B1 <- subset(seurat_obj_B1, subset = nFeature_Spatial > 200 & nCount_Spatial > 500)


# Visulization of QC results
p1 <- VlnPlot(seurat_obj_A1, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"), ncol = 3, pt.size = 0)
p2 <- SpatialFeaturePlot(seurat_obj_A1, features = "nCount_Spatial", pt.size.factor = 3, image.scale = "hires")+
  scale_fill_distiller(
    palette = "Spectral",
    labels = k_format  # using comma format for large numbers
    )

p3 <- SpatialFeaturePlot(seurat_obj_A1, features = "percent.mt", pt.size.factor = 3, image.scale = "hires")
ggsave("./results/QC_A1_after_qc.tiff", p1/p2/p3, width = 12, height = 24)


p1 <- VlnPlot(seurat_obj_B1, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"), ncol = 3, pt.size = 0)
p2 <- SpatialFeaturePlot(seurat_obj_B1, features = "nCount_Spatial", pt.size.factor = 2, image.scale = "hires") +
  scale_fill_distiller(
    palette = "Spectral",
    labels = k_format  # using comma format for large numbers
    )
p3 <- SpatialFeaturePlot(seurat_obj_B1, features = "percent.mt", pt.size.factor = 2, image.scale = "hires")
ggsave("./results/QC_B1_after_qc.tiff", p1/p2/p3, width = 12, height = 24)
# Save the filtered data
saveRDS(seurat_obj_A1, file = "./data/seurat_obj_A1_after_qc.rds")
saveRDS(seurat_obj_B1, file = "./data/seurat_obj_B1_after_qc.rds")
