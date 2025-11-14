# Load required libraries
library(Seurat)
library(SPOTlight)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(scran)
# Load the data
cat("Loading spatial and single-cell data...\n")
seurat_obj_A1 <- readRDS("./data/seurat_obj_A1.rds")
sc_seurat <- readRDS("./data/sc_seurat_annotated.rds")
sc_sce <- as.SingleCellExperiment(sc_seurat, assay = "RNA")
colData(sc_sce)$celltype <- sc_seurat$celltype_major
sc_sce <- logNormCounts(sc_sce)
spatial_sce <- as.SingleCellExperiment(seurat_obj_A1, assay = "Spatial")

# Get vector indicating which genes are neither ribosomal or mitochondrial
genes <- !grepl(pattern = "^Rp[l|s]|mt", x = rownames(sc_sce))

dec <- modelGeneVar(sc_sce, subset.row = genes)
plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
curve(metadata(dec)$trend(x), col = "blue", add = TRUE)
# Get the top 3000 genes.
hvg <- getTopHVGs(dec, n = 3000)
colLabels(sc_sce) <- colData(sc_sce)$celltype

# Compute marker genes
mgs <- scoreMarkers(sc_sce, subset.row = genes)
mgs_fil <- lapply(names(mgs), function(i) {
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.8
    x <- x[x$mean.AUC > 0.8, ]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)
# split cell indices by identity
idx <- split(seq(ncol(sc_sce)), sc_sce$celltype)
# downsample to at most 50 per identity & subset

# life analysis
n_cells <- 50
cs_keep <- lapply(idx, function(i) {
    n <- length(i)
    if (n < n_cells)
        n_cells <- n
    sample(i, n_cells)
})
sc_sce <- sc_sce[, unlist(cs_keep)]
res <- SPOTlight(
    x = sc_sce,
    y = spatial_sce,
    groups = as.character(sc_sce$celltype),
    mgs = mgs_df,
    hvg = hvg,
    weight_id = "mean.AUC",
    group_id = "cluster",
    gene_id = "gene")
# Extract NMF model fit
mod <- res$NMF
plotTopicProfiles(
    x = mod,
    y = sc_sce$celltype,
    facet = FALSE,
    min_prop = 0.01,
    ncol = 1) +
    theme(aspect.ratio = 1)
plotTopicProfiles(
    x = mod,
    y = sc_sce$celltype,
    facet = TRUE,
    min_prop = 0.01,
    ncol = 7)
mat <- res$mat
plotCorrelationMatrix(mat)
plotInteractions(mat, which = "heatmap", metric = "prop")
plotInteractions(mat, which = "heatmap", metric = "jaccard")
plotInteractions(mat, which = "network")
ct <- colnames(mat)
mat[mat < 0.1] <- 0

# Define color palette
# (here we use 'paletteMartin' from the 'colorBlindness' package)
paletteMartin <- c(
    "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db",
    "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff",
    "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")

pal <- colorRampPalette(paletteMartin)(length(ct))
names(pal) <- ct

 spatial_coords <- GetTissueCoordinates(seurat_obj_A1)

library(SpatialExperiment)
# create SpatialExperiment object
# reference: https://github.com/drighelli/SpatialExperiment/issues/115
## Extract and process image data
img <- SpatialExperiment::SpatialImage(
    x = as.raster(seurat_obj_A1@images[["slice1"]]@image))

imgData <- DataFrame(
    sample_id = "A1_colon_d0",
    image_id = "slice1",
    data = I(list(img)),
    scaleFactor = seurat_obj_A1@images[["slice1"]]@scale.factors$hires)
spatial_spe <- SpatialExperiment(
  assays = list(counts = counts(spatial_sce)),
  colData = colData(spatial_sce),
  spatialCoords = as.matrix(spatial_coords[, c("x", "y")]),
  sample_id = "A1_colon_d0",
  imgData = imgData
)
p1 <- plotSpatialScatterpie(
    x = spatial_spe,
    y = mat,
    cell_types = colnames(mat),
    img = FALSE,
    scatterpie_alpha = 1,
    pie_scale = 0.4
  ) +
    scale_fill_manual(
        values = pal,
        breaks = names(pal),
    )
p2 <- plotSpatialScatterpie(
    x = spatial_spe,
    y = mat,
    cell_types = colnames(mat),
    img = TRUE,
    scatterpie_alpha = 1,
    pie_scale = 0.4
  ) +
    scale_fill_manual(
        values = pal,
        breaks = names(pal),
    )
p <- p1 | p2
ggsave("./results/SPOTlight_deconvolution_A1.tiff", p, width = 16, height = 16)

# Add deconvolution results to spatial object
seurat_obj_A1 <- AddMetaData(seurat_obj_A1, as.data.frame(mat))

# Visualize major cell types
major_celltypes <- colnames(mat) # Top 6 cell types
spatial_plots <- list()

for (ct in major_celltypes) {
  p <- SpatialFeaturePlot(
    seurat_obj_A1,
    features = ct,
    pt.size.factor = 2,
    image.scale = "hires"
  ) +
    scale_fill_viridis_c() +
    ggtitle(paste("Cell Type:", ct))

  spatial_plots[[ct]] <- p
}

wrap_plots(spatial_plots, ncol = 3) %>%
  ggsave("./results/celltype_spatial_distribution_A1.tiff", ., width = 15, height = 10)

# Save the deconvoluted Seurat object
saveRDS(seurat_obj_A1, file = "./data/seurat_obj_A1_deconvoluted.rds")
cat("Deconvolution completed and results saved.\n")
# Save the mat
saveRDS(mat, file = "./data/SPOTlight_deconvolution_matrix_A1.rds")
#### all the same in B1 sample ----
# Load the data
cat("Loading spatial and single-cell data...\n")
seurat_obj_B1 <- readRDS("./data/seurat_obj_B1.rds")
sc_seurat <- readRDS("./data/sc_seurat_annotated.rds")
sc_sce <- as.SingleCellExperiment(sc_seurat, assay = "RNA")
colData(sc_sce)$celltype <- sc_seurat$celltype_major
sc_sce <- logNormCounts(sc_sce)
spatial_sce <- as.SingleCellExperiment(seurat_obj_B1, assay = "Spatial")

# Get vector indicating which genes are neither ribosomal or mitochondrial
genes <- !grepl(pattern = "^Rp[l|s]|mt", x = rownames(sc_sce))

dec <- modelGeneVar(sc_sce, subset.row = genes)
plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
curve(metadata(dec)$trend(x), col = "blue", add = TRUE)
# Get the top 3000 genes.
hvg <- getTopHVGs(dec, n = 3000)
colLabels(sc_sce) <- colData(sc_sce)$celltype

# Compute marker genes
mgs <- scoreMarkers(sc_sce, subset.row = genes)
mgs_fil <- lapply(names(mgs), function(i) {
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.8
    x <- x[x$mean.AUC > 0.8, ]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)
# split cell indices by identity
idx <- split(seq(ncol(sc_sce)), sc_sce$celltype)
# downsample to at most 20 per identity & subset
# We are using 5 here to speed up the process but set to 75-100 for your real
# life analysis
n_cells <- 50
cs_keep <- lapply(idx, function(i) {
    n <- length(i)
    if (n < n_cells)
        n_cells <- n
    sample(i, n_cells)
})
sc_sce <- sc_sce[, unlist(cs_keep)]
res <- SPOTlight(
    x = sc_sce,
    y = spatial_sce,
    groups = as.character(sc_sce$celltype),
    mgs = mgs_df,
    hvg = hvg,
    weight_id = "mean.AUC",
    group_id = "cluster",
    gene_id = "gene")
# Extract NMF model fit
mod <- res$NMF
plotTopicProfiles(
    x = mod,
    y = sc_sce$celltype,
    facet = FALSE,
    min_prop = 0.01,
    ncol = 1) +
    theme(aspect.ratio = 1)
plotTopicProfiles(
    x = mod,
    y = sc_sce$celltype,
    facet = TRUE,
    min_prop = 0.01,
    ncol = 7)
mat <- res$mat
plotCorrelationMatrix(mat)
plotInteractions(mat, which = "heatmap", metric = "prop")
plotInteractions(mat, which = "heatmap", metric = "jaccard")
plotInteractions(mat, which = "network")
ct <- colnames(mat)
mat[mat < 0.1] <- 0

# Define color palette
# (here we use 'paletteMartin' from the 'colorBlindness' package)
paletteMartin <- c(
    "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db",
    "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff",
    "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")

pal <- colorRampPalette(paletteMartin)(length(ct))
names(pal) <- ct

 spatial_coords <- GetTissueCoordinates(seurat_obj_B1)

library(SpatialExperiment)
# create SpatialExperiment object
# reference: https://github.com/drighelli/SpatialExperiment/issues/115
## Extract and process image data
img <- SpatialExperiment::SpatialImage(
    x = as.raster(seurat_obj_B1@images[["slice1"]]@image))

imgData <- DataFrame(
    sample_id = "B1_colon_d0",
    image_id = "slice1",
    data = I(list(img)),
    scaleFactor = seurat_obj_B1@images[["slice1"]]@scale.factors$hires)
spatial_spe <- SpatialExperiment(
  assays = list(counts = counts(spatial_sce)),
  colData = colData(spatial_sce),
  spatialCoords = as.matrix(spatial_coords[, c("x", "y")]),
  sample_id = "B1_colon_d0",
  imgData = imgData
)
p1 <- plotSpatialScatterpie(
    x = spatial_spe,
    y = mat,
    cell_types = colnames(mat),
    img = FALSE,
    scatterpie_alpha = 1,
    pie_scale = 0.4
  ) +
    scale_fill_manual(
        values = pal,
        breaks = names(pal),
    )
p2 <- plotSpatialScatterpie(
    x = spatial_spe,
    y = mat,
    cell_types = colnames(mat),
    img = TRUE,
    scatterpie_alpha = 1,
    pie_scale = 0.4
  ) +
    scale_fill_manual(
        values = pal,
        breaks = names(pal),
    )
p <- p1 | p2
ggsave("./results/SPOTlight_deconvolution_B1.tiff", p, width = 16, height = 16)

# Add deconvolution results to spatial object
seurat_obj_B1 <- AddMetaData(seurat_obj_B1, as.data.frame(mat))

# Visualize major cell types
major_celltypes <- colnames(mat) # Top 6 cell types
spatial_plots <- list()

for (ct in major_celltypes) {
  p <- SpatialFeaturePlot(
    seurat_obj_B1,
    features = ct,
    pt.size.factor = 2,
    image.scale = "hires"
  ) +
    scale_fill_viridis_c() +
    ggtitle(paste("Cell Type:", ct))

  spatial_plots[[ct]] <- p
}

wrap_plots(spatial_plots, ncol = 3) %>%
  ggsave("./results/celltype_spatial_distribution_B1.tiff", ., width = 15, height = 10)

# Save the deconvoluted Seurat object
saveRDS(seurat_obj_B1, file = "./data/seurat_obj_B1_deconvoluted.rds")
cat("Deconvolution completed and results saved.\n")
# Save the mat
saveRDS(mat, file = "./data/SPOTlight_deconvolution_matrix_B1.rds")
FeaturePlot(seurat_obj_A1, features = c("Fgf1", "Lgr5"), blend = TRUE, reduction = "Spatial", ) +
  scale_fill_viridis_c()


