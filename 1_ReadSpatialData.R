library(Seurat)
library(SeuratObject)

# 1. Read the data ---------------------------------------------------
# it's recommended to organize the data as follows if you have multiple samples when using Seurat:
# ./data/
# ├── sample_A1/
# │   ├── filtered_feature_bc_matrix.h5
# │   ├── spatial/
# │   │   ├── tissue_positions_list.csv
# │   │   ├── tissue_hires_image.png
# │   │   └── scalefactors_json.json
# │   └── ... (other files)
# └── sample_B1/
#     ├── filtered_feature_bc_matrix.h5
#     ├── spatial/
#     │   ├── tissue_positions_list.csv
#     │   ├── tissue_hires_image.png
#     │   └── scalefactors_json.json
#     └── ... (other files)
# # But since the there are only two samples here,
# we can alse read them directly one by one if you want to,
# or just **organize the data and skip** this region(recommended).
# # set the path of the data
# data_dir <- "./data/GSE169749_RAW"
# # read sample A1
# seurat_obj_A1 <- Seurat::Load10X_Spatial(
#   data.dir = data_dir,
#   filename = "GSM5213483_V19S23-097_A1_S1_filtered_feature_bc_matrix.h5", # you have to specify the h5 file name here in this one-by-one case
#   image = "GSM5213483_V19S23-097_A1_S1_tissue_hires_image.png.gz", # or use .tif.gz file
#   filter.matrix = TRUE,
#   assay = "Spatial",
#   slice = "slice_A1" # the name of the slice
# )
#
# # read sample B1
# seurat_obj_B1 <- Seurat::Load10X_Spatial(
#   data.dir = data_dir,
#   filename = "GSM5213484_V19S23-097_B1_S2_filtered_feature_bc_matrix.h5",
#   image = "tissue_hires_image.png",
#   filter.matrix = TRUE,
#   assay = "Spatial",
#   slice = "slice_B1"
# )
#
# # optional: add sample information
# seurat_obj_A1$sample <- "A1_colon_d0"
# seurat_obj_B1$sample <- "B1_colon_d14"
#
# # you can have a visulization of the spatial object
# # Reading the information of spot positions
# position_A1 <-  read.csv("./data/GSE169749_RAW/GSM5213483_V19S23-097_A1_S1_tissue_positions_list.csv.gz",header = FALSE,row.names = 1)
# position_B1 <-  read.csv("./data/GSE169749_RAW/GSM5213484_V19S23-097_B1_S2_tissue_positions_list.csv.gz",header = FALSE,row.names = 1)
# head(position_A1)
# position_A1 <- position_A1[, 4:5]
# position_B1 <- position_B1[, 4:5]
#
#
# # work around with no image available using DimReducObject, Or just organize the data and use the image file directly and skip here
# colnames(position_A1) = paste0("Spatial_",1:ncol(position_A1))
# # find common barcodes between seurat object and position file
# common_cells <- intersect(Cells(seurat_obj_A1), rownames(position_A1))
# cat("Common cells:", length(common_cells), "\n")
# # only keep the common cells
# position_A1 <- position_A1[common_cells, ]
# seurat_obj_A1[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(position_A1),
#                                                     key = "Spatial",assay = "RNA")
# DimPlot(seurat_obj_A1,reduction = "spatial",pt.size = 1)

# Read Organized Data ---------------------------------------------------
data_dir_A1 <- "./data/GSE169749_RAW/sample_A1"
data_dir_B1 <- "./data/GSE169749_RAW/sample_B1"

# Loading Image Manually and Load Spatial Data -----------------------------
image_A1 <- Seurat::Read10X_Image(image.dir = "./data/GSE169749_RAW/sample_A1/spatial",
                       image.name = "tissue_hires_image.png")#do it manually when using high resolution image file
seurat_obj_A1 <- Seurat::Load10X_Spatial(
  data.dir = data_dir_A1,
  image = image_A1  # use high resolution image file
)

image_B1 <- Seurat::Read10X_Image(image.dir = "./data/GSE169749_RAW/sample_B1/spatial",
                       image.name = "tissue_hires_image.png")
seurat_obj_B1 <- Seurat::Load10X_Spatial(
  data.dir = data_dir_B1,
   image = image_B1
)

Seurat::SpatialDimPlot(seurat_obj_A1,
               pt.size.factor = 3,
               image.scale = "hires" # important when using high resolution image
)

Seurat::SpatialDimPlot(seurat_obj_B1,
               pt.size.factor = 2,
               image.scale = "hires" # important when using high resolution image
)

cell_deconv <- read.csv("./data/download_page_VISDS000128_celltype_deco.csv", row.names = 1)
saveRDS(seurat_obj_A1, "./data/seurat_obj_A1.rds")
saveRDS(seurat_obj_B1, "./data/seurat_obj_B1.rds")

