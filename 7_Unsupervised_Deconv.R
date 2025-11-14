library(STdeconvolve)
library(Seurat)
# reference: https://github.com/JEFworks-Lab/STdeconvolve
seurat_obj_A1 <- readRDS("./data/seurat_obj_A1.rds")
# extract count matrix
count_matrix_A1 <- as.matrix(GetAssayData(seurat_obj_A1, assay = "Spatial", slot = "counts"))
# run STdeconvolve
counts <- cleanCounts(count_matrix_A1 , min.lib.size = 100)
# `pos` must have exactly 2 columns named `x` and `y`.
pos <- GetTissueCoordinates(seurat_obj_A1)
pos <- pos[, c('x', 'y')]
## feature select for genes
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
## choose optimal number of cell-types
ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(2, 9, by = 1))
## get best model results
optLDA <- optimalModel(models = ldas, opt = "min")
## extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta
saveRDS(deconProp, file = "./data/STdeconvolve_A1_celltype_proportions.rds")
annot <- read.csv("./data/sss_result_A1.csv", row.names = 1)
annot <- annot[rownames(deconProp),]
group <- annot$Custom.Name
group <- factor(group)
## visualize deconvolved cell-type proportions
tiff("./results/STdeconvolve_A1_celltype_proportions.tiff", width = 8, height = 6, units = 'in', res = 300)
vizAllTopics(deconProp, pos,
             groups = group,
             group_cols = rainbow(length(levels(group))),
             r=40,
             lwd = 0.2
             )
dev.off()



