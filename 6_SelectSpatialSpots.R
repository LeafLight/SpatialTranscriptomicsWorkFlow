library(Seurat)
# Create data required for SelectSpatialSpots
seurat_obj_A1 <- readRDS("./data/seurat_obj_A1.rds")
position_A1 <- GetTissueCoordinates(seurat_obj_A1)
position_A1 <- position_A1[, c('cell', 'x', 'y')]
colnames(position_A1) <- c("CELL_ID", "X", "Y")
write.csv(position_A1, "./data/sss_A1.csv")


