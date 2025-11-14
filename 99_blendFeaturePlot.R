library(Seurat)
seurat_obj_A1 <- readRDS("./data/seurat_obj_A1.rds")
position_A1 <- read.csv("./data/GSE169749_RAW/GSM5213483_V19S23-097_A1_S1_tissue_positions_list.csv.gz",header = FALSE,row.names = 1)
position_A1 <- position_A1[, 4:5]
colnames(position_A1) <-  paste0("Spatial_",1:ncol(position_A1))
# find common barcodes between seurat object and position file
common_cells <- intersect(Cells(seurat_obj_A1), rownames(position_A1))
cat("Common cells:", length(common_cells), "\n")
# only keep the common cells

# Attention: all the params below is required to be specified
position_A1 <- position_A1[common_cells, ]
seurat_obj_A1[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(position_A1),
                                                    key = "Spatial",assay = "RNA")
FeaturePlot(seurat_obj_A1, features = c("Epcam", "Vim"), blend = TRUE, blend.threshold = 0.5,
                          reduction = "spatial", pt.size = 1, slot="count")

