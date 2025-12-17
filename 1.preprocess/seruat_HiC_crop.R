set.seed(123)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
# library(hdf5r)
library(Matrix)
library(sctransform)
library(plyr)
library(gridExtra)
library(magrittr)
library(tidyr)
# library(raster)
library(OpenImageR)
library(ggpubr)
library(grid)
library(wesanderson)
library(png)
library(pdftools)
library(cowplot)
library(magick)
# library(DBITtools)
library(SPtoolbox)
library(purrr)

sample <- "E13-18-HiC"  # "E1305-hic-24"
pixels <- 2500 # 2500
res_col <- 4812  # res=1Mb

# 1 Load data from higashi
setwd(paste0("/media/maoni/data/CZP/spatial_hic/", sample, "/"))
stat_clean <- read.table("5_bedpe/stat_clean.csv", sep =",", header = TRUE, dec =".", stringsAsFactors = F)
iA_iB <- paste(stat_clean$iA, stat_clean$iB, sep="x")   ## for sample E13-18-HiC (正常顺序)
# iA_iB <- paste(51-stat_clean$iA, stat_clean$iB, sep="x")  ## for sample E1305-hic-24
res_col_names <- read.table("/media/maoni/data/CZP/spatial_hic/Res_500kb_col_names_4812.txt", sep="\t", header=FALSE)

count.m <- matrix(0, pixels, res_col)
rownames(count.m) <- iA_iB
colnames(count.m) <- res_col_names[,1]
dim(count.m) # 9216 2413

# 2 Load image information
assay = "Spatial"
slice = sample
Seurat_object = CreateSeuratObject(counts = t(count.m), project = sample, assay = assay) # min.cells min.features

image.dir = "/media/maoni/data/CZP/spatial_hic/seurat/core/50*50_blank_img" # "./Img"
image.nam = paste0(sample,"_fix.png") # "grey_pixel_1080p.png"
coord.nam = "combine_barcode.round2round1_index1_index2.Seurat.txt"
image <- readPNG(source = file.path(image.dir, image.nam))[,,1:3]
scale.factors <- c("tissue_hires_scalef"=1, "fiducial_diameter_fullres"=1, "tissue_lowres_scalef"=1)
tissue.positions <- read.table(file = file.path(image.dir,coord.nam), col.names = c("barcodes", "tissue", "row", "col", "imagerow", "imagecol"), header = FALSE, as.is = TRUE, row.names = 1)
spot.radius <- 0.015 # estiamte:(0.13)*50/410/2
image <- new(Class = "VisiumV1", image = image, scale.factors = scalefactors(spot = scale.factors[1], fiducial = scale.factors[2], hires = scale.factors[1], lowres = scale.factors[3]), coordinates = tissue.positions, spot.radius = spot.radius)
image <- image[Cells(Seurat_object)]
DefaultAssay(object = image) <- assay
Seurat_object[[slice]] <- image


# 3. crop
Seurat_object$crop <- "original"

Seurat_object_cropped = shiny_st(Seurat_object, assay = "Spatial", slot = "data", image = Images(Seurat_object), python_env = "/home/maoni/miniconda3/envs/r4.3/bin/python", script = "/media/maoni/data/Software/DBIT_tools/filter_pixel_AI.py")

table(Seurat_object_cropped$crop)

save(Seurat_object_cropped, file = paste0("6_crop/", sample, ".crop.seurat.RData"))

stat_clean <- read.table("5_bedpe/stat_clean.csv", sep =",", header = TRUE, dec =".", stringsAsFactors = F)
iA_iB = paste(stat_clean$iA, stat_clean$iB, sep="x")  ####
stat_cropped <- stat_clean[match(colnames(Seurat_object_cropped)[which(Seurat_object_cropped$crop=="original")],iA_iB), ]
write.csv(stat_cropped, "6_crop/stat_cropped.csv", sep = ",", col.names = TRUE, row.names = FALSE)
dim(stat_cropped)

stat_cropped_barcode <- data.frame(stat_cropped[,5])
write.table(stat_cropped_barcode,  "6_crop/stat_cropped_barcode", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# 4. check crop
load(paste0("6_crop/", sample, ".crop.seurat.RData"))

Seurat_object_cropped <- subset(Seurat_object_cropped, crop == "original")

pdf("6_crop/0_nCounts.Seurat.pdf", width = 14, height = 7)
plot1 <- VlnPlot(Seurat_object_cropped, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(Seurat_object_cropped, features = "nCount_Spatial", crop = F, max.cutoff = 2000, pt.size.factor = 0.8, stroke = NA) + theme(legend.position = "right")
plot2$layers[[1]]$aes_params=c(plot2$layers[[1]]$aes_params, shape=22)
plot_grid(plot1, plot2)
dev.off()

pdf("6_crop/0_nCounts.Seurat_filter.pdf", width = 14, height = 7)
plot1 <- SpatialFeaturePlot(Seurat_object_cropped, features = "nCount_Spatial",crop = F,max.cutoff = 2000, pt.size.factor = 0.8, stroke = NA) + theme(legend.position = "right")
plot1$layers[[1]]$aes_params=c(plot1$layers[[1]]$aes_params, shape=22)
plot2 <- SpatialFeaturePlot(Seurat_object_cropped, features = "nFeature_Spatial",crop = F,max.cutoff = 1500, pt.size.factor = 0.8, stroke = NA) + theme(legend.position = "right")
plot2$layers[[1]]$aes_params=c(plot2$layers[[1]]$aes_params, shape=22)
plot_grid(plot1, plot2)
dev.off()


