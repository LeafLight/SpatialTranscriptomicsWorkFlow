library(Giotto)
library(Seurat)
seurat_obj_A1 <- readRDS("./data/seurat_obj_A1_deconvoluted.rds")
seurat_obj_B1 <- readRDS("./data/seurat_obj_B1_deconvoluted.rds")
fromSeuratToGiotto <- function(obj){
expr_matrix <- GetAssayData(obj, assay = "Spatial", slot = "counts")
# 提取空间坐标信息
coords <- GetTissueCoordinates(obj)
giotto_obj <- createGiottoObject(expression = expr_matrix, spatial_locs = coords)
}
giotto_A1 <- fromSeuratToGiotto(seurat_obj_A1)
spatPlot(giotto_A1, show_image = FALSE, point_alpha = 0.7,
         save_param = list(save_name = "spatplot_with_image"))
# create spatial network
giotto_A1  <- createSpatialNetwork(giotto_A1 ,
                                   method = "Delaunay",
                                   name = "Delaunay_network")
# loading lr network data(ligand-receptor)
# reference: https://zenodo.org/api/records/15168114/files/lr_network_mouse.csv/content

prop_mat <- readRDS("./data/SPOTlight_deconvolution_matrix_A1.rds")
celltype_major <- colnames(prop_mat)[max.col(prop_mat, ties.method = "first")]
cell_metadata <- pDataDT(giotto_A1)

cell_metadata[, celltype_major := celltype_major]
# add meta to giotto object
giotto_A1 <- addCellMetadata(
  gobject = giotto_A1,
  new_metadata = cell_metadata,
  by_column = TRUE,
  column_cell_ID = "cell_ID"
)
giotto_A1  <- normalizeGiotto(giotto_A1)
lr_net <- read.csv("./data/mouse_lr_network.csv")
lr_from <- lr_net$from[lr_net$from %in% rownames(giotto_A1) & lr_net$to %in% rownames(giotto_A1)]
lr_to <- lr_net$to[lr_net$from %in% rownames(giotto_A1) & lr_net$to %in% rownames(giotto_A1)]
ccc_res <- exprCellCellcom(
  giotto_A1,
  cluster_column = "celltype_major",
  feat_set_1=lr_from,
  feat_set_2=lr_to,
  verbose = T
)
p <- ggplot(ccc_res[(p.adj < 0.05) & log2fc>2, ], aes(x = LR_cell_comb, y = LR_comb, size = log2fc, color = LR_cell_comb)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Custom Cell-Cell Communication Plot",
       x = "Cell Types",
       y = "Ligand-Receptor") +
  scale_size_continuous(range = c(1, 10))+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./results/Giotto_A1_cellcell_communication.tiff", p, width = 8, height = 24)


giotto_B1 <- fromSeuratToGiotto(seurat_obj_B1)
spatPlot(giotto_B1, show_image = FALSE, point_alpha = 0.7,
         save_param = list(save_name = "spatplot_with_image"))
# create spatial network
giotto_B1  <- createSpatialNetwork(giotto_B1 ,
                                   method = "Delaunay",
                                   name = "Delaunay_network")
# loading lr network data(ligand-receptor)
# reference: https://zenodo.org/api/records/15168114/files/lr_network_mouse.csv/content

prop_mat <- readRDS("./data/SPOTlight_deconvolution_matrix_B1.rds")
celltype_major <- colnames(prop_mat)[max.col(prop_mat, ties.method = "first")]
cell_metadata <- pDataDT(giotto_B1)

cell_metadata[, celltype_major := celltype_major]
# add meta to giotto object
giotto_B1 <- addCellMetadata(
  gobject = giotto_B1,
  new_metadata = cell_metadata,
  by_column = TRUE,
  column_cell_ID = "cell_ID"
)
giotto_B1  <- normalizeGiotto(giotto_B1)
lr_net <- read.csv("./data/mouse_lr_network.csv")
lr_from <- lr_net$from[lr_net$from %in% rownames(giotto_B1) & lr_net$to %in% rownames(giotto_B1)]
lr_to <- lr_net$to[lr_net$from %in% rownames(giotto_B1) & lr_net$to %in% rownames(giotto_B1)]
ccc_res <- exprCellCellcom(
  giotto_B1,
  cluster_column = "celltype_major",
  feat_set_1=lr_from,
  feat_set_2=lr_to,
  verbose = T
)
p <- ggplot(ccc_res[(p.adj < 0.05) & log2fc>1, ], aes(x = LR_cell_comb, y = LR_comb, size = log2fc, color = LR_cell_comb)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Custom Cell-Cell Communication Plot",
       x = "Cell Types",
       y = "Ligand-Receptor") +
  scale_size_continuous(range = c(1, 10))+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./results/Giotto_B1_cellcell_communication.tiff", p, width = 14, height = 24)
