library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(sctransform)
library(plyr)
library(gridExtra)
library(magrittr)
library(tidyr)
library(OpenImageR)
library(ggpubr)
library(grid)
library(wesanderson)
library(png)
library(cowplot)
library(viridis)
library(qpdf)
library(RColorBrewer)
# crop packages
# conda install -c conda-forge imagemagick
# BiocManager::install("magick")
# devtools::install_github("EddieLv/STvis/STvis", force = TRUE)
library(magick) # 
# library(STvis) # 

# 0 file preparing
source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_hic/E1305-hic-24"
sample <- "E1305-hic-24"
higashi_out <- "8_higashi_out_higashi_500kb"
work_path_seurat <- "10_seurat_bigSpotShape_higashi_500kb"
embed_type <- "scAB_scale"

# "scAB_scale" "embedding", embed_all","embed_raw","embed_l2_norm","embed_correct_coverage_fh","embed_l2_norm_correct_coverage_fh"
setwd(work_path)
dir.create(work_path_seurat)
setwd(paste0(work_path, "/", work_path_seurat))
dir.create(embed_type)
setwd(paste0(work_path,  "/", work_path_seurat, "/", embed_type))
dir.create("Img")
dir.create("Matrix")
dir.create("Out")
dir.create("Plot")
file.copy(from=paste0(work_path, "/5_bedpe/stat_clean.csv"), 
          to=paste0(work_path,  "/", work_path_seurat, "/", embed_type, "/Matrix/"), overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)

file.copy(from=paste0(work_path, "/6_crop/stat_cropped.csv"), 
          to=paste0(work_path,  "/", work_path_seurat, "/", embed_type, "/Matrix/"), overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)

file.copy(from=paste0(work_path, "/", higashi_out, "/", sample, "_scAB_scale.csv"), 
          to=paste0(work_path,  "/", work_path_seurat, "/", embed_type, "/Matrix/"), overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)

file.copy(from=paste0(work_path, "/", higashi_out, "/", sample, "_", embed_type, ".csv"), 
          to=paste0(work_path,  "/", work_path_seurat, "/", embed_type, "/Matrix/"), overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)


#################################################################### E13-24 ##########################################################
# 1 Load data from st pipe
setwd(paste0(work_path, "/", work_path_seurat, "/", embed_type))
stat_cropped <- read.table("Matrix/stat_cropped.csv", sep =",", header = TRUE, dec =".", stringsAsFactors = F)
iA_iB <- paste(51-stat_cropped$iA, stat_cropped$iB, sep="x")
count_matrix <- read.table(paste0("Matrix/", sample, "_scAB_scale.csv"), sep =",", header = F, dec =".", stringsAsFactors = F)
rownames(count_matrix) <- iA_iB
count_matrix <- as.matrix(count_matrix)
dim(count_matrix) # 1999 2413/4813

bin_spatial <- read.table(paste0(work_path ,"/", higashi_out, "/", sample, "_scAB_bin.csv"), sep =",", header = FALSE, dec =".", stringsAsFactors = F)
bin_spatial[,1] <-  gsub("b","",bin_spatial$V1)
bin_spatial_merge <- paste(bin_spatial$V1, paste(bin_spatial$V3, bin_spatial$V2, sep="-"), sep=":")
colnames(count_matrix) <- bin_spatial_merge

emb.m <- as.matrix(read.table(paste0("Matrix/", sample, "_", embed_type, ".csv"), sep=",", header = F))
rownames(emb.m) <- rownames(count_matrix)
dim(emb.m)  # 1999  64/256

assay = "Spatial"
slice = sample
image.dir = "/media/maoni/data/CZP/spatial_hic/seurat_yuhao/core/50*50_blank_img" # "./Img"
image.nam = paste0(sample, "_fix.png") # "white_background_1080p.png" # "grey_pixel_1080p.png"
coord.nam = "combine_barcode.round2round1_index1_index2.Seurat.txt"
Seurat_object = CreateSeuratObject(counts = t(count_matrix), project = sample, assay = assay) # min.cells min.features
image <- readPNG(source = file.path(image.dir, image.nam))[,,1:3]
scale.factors <- c("tissue_hires_scalef"=1, "fiducial_diameter_fullres"=1, "tissue_lowres_scalef"=1)
tissue.positions <- read.table(file = file.path(image.dir,coord.nam), col.names = c("barcodes", "tissue", "row", "col", "imagerow", "imagecol"), header = FALSE, as.is = TRUE, row.names = 1)
spot.radius <- 0.015 # estiamte:(0.13)*50/410/2
image <- new(Class = "VisiumV1", image = image, scale.factors = scalefactors(spot = scale.factors[1], fiducial = scale.factors[2], hires = scale.factors[1], lowres = scale.factors[3]), coordinates = tissue.positions, spot.radius = spot.radius)
image <- image[Cells(Seurat_object)]
DefaultAssay(object = image) <- assay
Seurat_object[[slice]] <- image

# 3 Data quality
Seurat_object$valid_contact <- log10(stat_cropped$valid_contact)
Seurat_object$cis_more_10kb <- log10(stat_cropped$cis_more_10kb)

pdf("Plot/0_nCounts.Seurat_valid_contact_cis_more_10kb.pdf", width = 14, height = 7)
plot1 <- SpatialFeaturePlot(Seurat_object, features = "valid_contact",crop = F, pt.size.factor = 1.6, stroke = NA, min.cutoff = "q3") + theme(aspect.ratio = 1)
plot1$layers[[1]]$aes_params=c(plot1$layers[[1]]$aes_params, shape=22)
plot2 <- SpatialFeaturePlot(Seurat_object, features = "cis_more_10kb",crop = F, pt.size.factor = 1.6, stroke = NA, min.cutoff = "q3") + theme(aspect.ratio = 1) 
plot2$layers[[1]]$aes_params=c(plot2$layers[[1]]$aes_params, shape=22)
plot_grid(plot1, plot2)
dev.off()

pdf("Plot/0_nCounts.Seurat.pdf", width = 14, height = 7)
plot1 <- VlnPlot(Seurat_object, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(Seurat_object, features = "nCount_Spatial", max.cutoff = 2000, pt.size.factor = 1.6, stroke = NA) + theme(legend.position = "right")
plot2$layers[[1]]$aes_params=c(plot2$layers[[1]]$aes_params, shape=22)
plot_grid(plot1, plot2)
dev.off()

pdf("Plot/0_nCounts.Seurat_filter.pdf", width = 14, height = 7)
plot1 <- SpatialFeaturePlot(Seurat_object, features = "nCount_Spatial",crop = F,max.cutoff = 2000, pt.size.factor = 1.6, stroke = NA) + theme(legend.position = "right")
plot1$layers[[1]]$aes_params=c(plot1$layers[[1]]$aes_params, shape=22)
plot2 <- SpatialFeaturePlot(Seurat_object, features = "nFeature_Spatial",crop = F,max.cutoff = 1500, pt.size.factor = 1.6, stroke = NA) + theme(legend.position = "right")
plot2$layers[[1]]$aes_params=c(plot2$layers[[1]]$aes_params, shape=22)
plot_grid(plot1, plot2)
dev.off()

# 4 Run PCA
Seurat_object <- SCTransform(Seurat_object, assay = "Spatial", verbose = T, variable.features.n = 2000)
Seurat_object <- RunPCA(Seurat_object, assay = "SCT", verbose = FALSE)
ElbowPlot(Seurat_object, ndims = 50)
# Seurat_object[["pca"]] = CreateDimReducObject(embeddings = emb.m, key = "PC_", assay = 'SCT')

 
# 5 Cluster
dim_p <- list()
spatial_p <- list()
pcNO <- 36
for(res in c(0.1,0.3,0.5,0.7,0.9,1,1.2,1.5,2)){
  Seurat_object <- FindNeighbors(Seurat_object, reduction = "pca", dims = 1:pcNO)
  Seurat_object <- FindClusters(Seurat_object, verbose = FALSE, resolution = res)
  Seurat_object <- RunUMAP(Seurat_object, reduction = "pca", dims = 1:pcNO)
  dim_p[[paste0("res",res)]] <- DimPlot(Seurat_object, reduction = "umap", label = TRUE) + ggtitle(paste0("res",res)) # the UMAP plot
  spatial_p[[paste0("res",res)]] <- SpatialDimPlot(Seurat_object, label = TRUE, label.size = 3, pt.size.factor = 1.6, stroke = NA) + ggtitle(paste0("res",res)) # the spatial plot
  spatial_p[[paste0("res",res)]]$layers[[1]]$aes_params=c(spatial_p[[paste0("res",res)]]$layers[[1]]$aes_params, shape=22)
}

pdf(file = paste("Plot/2_UMAP_clusters_PC",pcNO,".pdf",sep =""), width = 35, height = ceiling(length(dim_p)/5)*7)
ggarrange(plotlist = dim_p, ncol = 5, nrow = ceiling(length(dim_p)/5))
dev.off()

pdf(file = paste("Plot/2_Spatial_clusters_PC",pcNO,".pdf",sep =""), width = 35, height = ceiling(length(spatial_p)/5)*7)
ggarrange(plotlist = spatial_p, ncol = 5, nrow = ceiling(length(spatial_p)/5))
dev.off()


dim_p <- list()
spatial_p <- list()
res <- 1
for(pcNO in seq(10,50,by=5)){
  Seurat_object <- FindNeighbors(Seurat_object, reduction = "pca", dims = 1:pcNO)
  Seurat_object <- FindClusters(Seurat_object, verbose = FALSE, resolution = res)
  Seurat_object <- RunUMAP(Seurat_object, reduction = "pca", dims = 1:pcNO)
  dim_p[[paste0("PC",pcNO)]] <- DimPlot(Seurat_object, reduction = "umap", label = TRUE) + ggtitle(paste0("PC",pcNO)) # the UMAP plot
  spatial_p[[paste0("PC",pcNO)]] <- SpatialDimPlot(Seurat_object, label = TRUE, label.size = 3, pt.size.factor = 1.6, stroke = NA) + ggtitle(paste0("PC",pcNO)) # the spatial plot
  spatial_p[[paste0("PC",pcNO)]]$layers[[1]]$aes_params=c(spatial_p[[paste0("PC",pcNO)]]$layers[[1]]$aes_params, shape=22)
}

pdf(file = paste("Plot/2_UMAP_clusters_res",res,".pdf",sep =""), width = 28, height = ceiling(length(dim_p)/4)*7)
ggarrange(plotlist = dim_p, ncol = 4, nrow = ceiling(length(dim_p)/4))
dev.off()

pdf(file = paste("Plot/2_Spatial_clusters_res",res,".pdf",sep =""), width = 28, height = ceiling(length(spatial_p)/4)*7)
ggarrange(plotlist = spatial_p, ncol = 4, nrow = ceiling(length(spatial_p)/4))
dev.off()


# Select PC and res
PC=36
res=1
Seurat_object <- FindNeighbors(Seurat_object, reduction = "pca", dims = 1:PC)
Seurat_object <- FindClusters(Seurat_object, verbose = FALSE, resolution = res)
Seurat_object <- RunUMAP(Seurat_object, reduction = "pca", dims = 1:PC)
umap_Spa <- list()
umap_Spa[["p1"]] <- DimPlot(Seurat_object, reduction = "umap", label = TRUE)
umap_Spa[["p2"]] <- SpatialDimPlot(Seurat_object, label = TRUE, label.size = 3, pt.size.factor = 1.6, stroke = NA)
umap_Spa[["p2"]]$layers[[1]]$aes_params=c(umap_Spa[["p2"]]$layers[[1]]$aes_params, shape=22)
umap_Spa[["p3"]] <- DimPlot(Seurat_object, reduction = "umap", label = TRUE)
umap_Spa[["p4"]] <- SpatialDimPlot(Seurat_object, label = TRUE, label.size = 3, pt.size.factor = 1.6, stroke = NA)
umap_Spa[["p4"]]$layers[[1]]$aes_params=c(umap_Spa[["p4"]]$layers[[1]]$aes_params, shape=22)

pdf(file = paste0("Plot/2_UMAP_clusters_PC",PC,"_res",res,".pdf"), width = 28, height = ceiling(length(umap_Spa)/4)*7)
p <- ggarrange(plotlist = umap_Spa, ncol = 4, nrow = ceiling(length(umap_Spa)/4))
print(p)
dev.off()

spatial_p_cluster <- list()
for(i in sort(unique(Idents(Seurat_object)))){
  spatial_p_cluster[[paste0("idents_",i)]] <- SpatialDimPlot(Seurat_object, cells.highlight = CellsByIdentities(object = Seurat_object, idents = i), cols.highlight = c("#DE2D26", "grey90"), facet.highlight = TRUE, ncol = 1, pt.size.factor = 1.6, stroke = NA, label.size = 12) + ggtitle(paste0("idents",i)) # the spatial plot
  spatial_p_cluster[[paste0("idents_",i)]]$layers[[1]]$aes_params=c(spatial_p_cluster[[paste0("idents_",i)]]$layers[[1]]$aes_params, shape=22)
}

pdf(file = paste0("Plot/2_Spatial_clusters_PC",PC,"_res",res,".pdf"), width = 28, height = ceiling(length(spatial_p_cluster)/4)*7)
p <- ggarrange(plotlist = spatial_p_cluster, ncol = 4, nrow = ceiling(length(spatial_p_cluster)/4))
print(p)
dev.off()

pdf_combine(input = c(paste0("Plot/2_UMAP_clusters_PC",PC,"_res",res,".pdf"), paste0("Plot/2_Spatial_clusters_PC",PC,"_res",res,".pdf")), output = paste0("Plot/2_UMAP_Spatial_clusters_PC",PC,"_res",res,".pdf"))


# 6.DE
de_markers_24 = FindAllMarkers(Seurat_object, test.use = "wilcox", only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
de_markers_24 %>% group_by(cluster) %>% slice_max(n=5, order_by = avg_log2FC) -> top5

DefaultAssay(Seurat_object) <- "Spatial"

pdf("Plot/3_findmarker_2D_impose.data_pc36_res1.pdf",height = length(table(top5$cluster))*7, width = 5*7)
plot <- SpatialFeaturePlot(object = Seurat_object, features = top5$gene, slot = "counts", alpha = c(0.1, 1), ncol = 5, pt.size.factor = 1.6)
plot$layers[[1]]$aes_params=c(plot$layers[[1]]$aes_params, shape=22)
plot
dev.off()

# scale_fill_gradientn(colours = viridis::inferno(100))
# scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")))
# scale_fill_gradientn(values=c(0,0.1,0.2,0.4,0.6,0.8,1),colors = c("grey86","grey86","#fed976","#fc4e2a","#e31a1c","#bd0026","#800026"), na.value = "white")
# scale_color_gradientn(values=c(0,0.0001,0.2,0.4,0.5,0.7,0.8,1),colors = c("#EEE595","#E6E18F","#E6D27D","#E7B553","#D1742E","#84235D","#452458","#100818"))
# scale_fill_gradientn(values=c(0,0.1,0.2,0.4,0.6,0.8,1),colors = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), na.value = "white")


pdf("Plot/3_findmarker_UMAP.data_pc36_res1.pdf",height = length(table(top5$cluster))*7, width = 5*7)
FeaturePlot(Seurat_object, features = top5$gene,ncol = 5,pt.size = 0.1,order = T)
dev.off()

pdf(file = paste("Plot/3_top5_ABcom_heat.pdf"), width = 20, height = 15)
DoHeatmap(Seurat_object, features = top5$gene, slot = "data", group.by = "seurat_clusters", size = 10) + scale_fill_gradientn(colors = c("white", "grey", "firebrick3")) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank())
dev.off()

pdf(file = paste("Plot/3_allDE_ABcom_heat.pdf"), width = 20, height = 15)
DoHeatmap(Seurat_object, features = de_markers_24$gene, slot = "data", group.by = "seurat_clusters", size = 10) + scale_fill_gradientn(colors = c("white", "grey", "firebrick3")) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank())
dev.off()


write.csv(de_markers_24,row.names = F,file = "Out/1_de_ABcom.csv")
write.csv(Seurat_object@meta.data,row.names = T,file = "Out/1_cell_metadata.csv")
save(de_markers_24, top5, Seurat_object, file=paste0("Out/", sample, "_Seurat_object.RData"))


# 7. annotaton
setwd(paste0(work_path,  "/", work_path_seurat, "/", embed_type))
load(paste0("Out/", sample, "_Seurat_object.RData"))
PC=36
res=1

anno <- rep(NA,ncol(Seurat_object))
anno[Seurat_object$seurat_clusters %in% c(0,2,4,8)] <- "Other"
anno[Seurat_object$seurat_clusters %in% c(1,3)] <- "CNS"
anno[Seurat_object$seurat_clusters %in% c(5)] <- "Radial glia"
anno[Seurat_object$seurat_clusters %in% c(6)] <- "Liver"
anno[Seurat_object$seurat_clusters %in% c(7)] <- "Heart"

Seurat_object$cell_type <- factor(anno,levels = c("CNS", "Radial glia", "Heart", "Liver", "Other"))

my_color_palette <- brewer.pal(9, "Set1")[c(4,6,8,7,9)]
names(my_color_palette) <- c("CNS", "Radial glia", "Heart", "Liver", "Other")

idents_cols_raw <- brewer.pal(9, "Set1")
names(idents_cols_raw) <- levels(Seurat_object$seurat_clusters)


umap_Spa <- list()
umap_Spa[["p1"]] <- DimPlot(Seurat_object, reduction = "umap", group.by = "seurat_clusters", label = FALSE, pt.size = 0.5, cols = idents_cols_raw) + theme(aspect.ratio = 1)
umap_Spa[["p2"]] <- SpatialDimPlot(Seurat_object, group.by = "seurat_clusters", label = FALSE, label.size = 3, pt.size.factor = 1.6, stroke = NA, cols = idents_cols_raw) + theme(aspect.ratio = 1)
umap_Spa[["p2"]]$layers[[1]]$aes_params=c(umap_Spa[["p2"]]$layers[[1]]$aes_params, shape=22)
umap_Spa[["p3"]] <- DimPlot(Seurat_object, reduction = "umap", group.by = "seurat_clusters", label = FALSE, pt.size = 0.5, cols = idents_cols_raw)+ theme(aspect.ratio = 1)
umap_Spa[["p4"]] <- SpatialDimPlot(Seurat_object, group.by = "seurat_clusters", label = FALSE, label.size = 3, pt.size.factor = 1.6, stroke = NA, cols = idents_cols_raw)+ theme(aspect.ratio = 1)
umap_Spa[["p4"]]$layers[[1]]$aes_params=c(umap_Spa[["p4"]]$layers[[1]]$aes_params, shape=22)

pdf(file = paste0("Plot/2_UMAP_clusters_PC",PC,"_res",res,"_final.pdf"), width = 28, height = ceiling(length(umap_Spa)/4)*7)
p <- ggarrange(plotlist = umap_Spa, ncol = 4, nrow = ceiling(length(umap_Spa)/4))
print(p)
dev.off()

spatial_p_cluster <- list()
for(i in names(table(Seurat_object$seurat_clusters))){
  Seurat_object_subset <- Seurat_object[, Seurat_object$seurat_clusters==i]
  spatial_p_cluster[[i]] <- SpatialDimPlot(Seurat_object_subset, group.by = "seurat_clusters", crop = FALSE, label = FALSE, label.size = 3, pt.size.factor = 1.6, stroke = NA,  label.box = FALSE, cols = idents_cols_raw) + theme(aspect.ratio = 1) + ggtitle(i) + theme(legend.position="none") # the spatial plot
  spatial_p_cluster[[i]]$layers[[1]]$aes_params=c(spatial_p_cluster[[i]]$layers[[1]]$aes_params, shape=22)
}

pdf(file = paste0("Plot/2_Spatial_clusters_PC",PC,"_res",res,"_final.pdf"), width = 28, height = ceiling(length(spatial_p_cluster)/4)*7)
p <- ggarrange(plotlist = spatial_p_cluster, ncol = 4, nrow = ceiling(length(spatial_p_cluster)/4))
print(p)
dev.off()

pdf_combine(input = c(paste0("Plot/2_UMAP_clusters_PC",PC,"_res",res,"_final.pdf"), paste0("Plot/2_Spatial_clusters_PC",PC,"_res",res,"_final.pdf")), output = paste0("Plot/2_UMAP_Spatial_clusters_PC",PC,"_res",res,"_final.pdf"))

save(Seurat_object, file=paste0("Out/", sample, "_Seurat_object_anno.RData"))

# 8. hightlight
source("/media/maoni/data/R_functions_ST_E13.R")
Seurat_object_magic <- magic.Seurat(Seurat_object)
DefaultAssay(Seurat_object_magic) <- "MAGIC_Spatial"

save(Seurat_object_magic, file=paste0("Out/", sample, "_Seurat_object_magic.RData"))

load(paste0("Out/", sample, "_Seurat_object_magic.RData"))

pdf("Plot/3_findmarker_2D_impose.data_pc36_res1_magic.pdf",height = length(table(top5$cluster))*7, width = 5*7)
plot <- SpatialFeaturePlot(object = Seurat_object_magic, features = top5$gene, alpha = c(1, 1), ncol = 5, pt.size.factor = 1.6) + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")))
plot$layers[[1]]$aes_params=c(plot$layers[[1]]$aes_params, shape=22)
plot
dev.off()

de_spe_select <- top5$gene[c(6,27,34,37,7,28,35,38)]
# [1] "chr1:68500000-69000000"   "chr11:23000000-23500000"  "chr13:12000000-12500000"  "chr5:109000000-109500000" "chr1:68000000-68500000"   "chr2:150500000-151000000" "chr7:18500000-19000000"  
# [8] "chr4:143500000-144000000"
de_spe_select <- c("chr19:48000000-48500000", "chr19:48500000-49000000", "chr2:162000000-162500000", "chr15:41000000-41500000", "chr1:158000000-158500000", "chr1:158500000-159000000", "chr9:28500000-29000000", "chr9:29000000-29500000", "chr9:29500000-30000000")

de_spe_select <- c("chr19:47500000-48000000", "chr15:40500000-41000000", "chr16:5500000-6000000")

de_spe_select <- c("chr15:40500000-41000000")

spatial_p_feature <- list()
for(i in 1:length(de_spe_select)){
  spatial_p_feature[[paste0("feature",i)]] <- SpatialFeaturePlot(object = Seurat_object_magic, features = de_spe_select[i], alpha = c(1, 1), ncol = 1, pt.size.factor = 1.6, stroke = NA) + scale_fill_gradientn(values=c(0,0.3,0.5,0.6,0.7,0.8,1),colors = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407")[c(1:4,7)], na.value = "white") # + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) #  + theme(legend.position="right")
  spatial_p_feature[[paste0("feature",i)]]$layers[[1]]$aes_params=c(spatial_p_feature[[paste0("feature",i)]]$layers[[1]]$aes_params, shape=22)
}

pdf(file = paste("Plot/4_findmarker_2D_impose.data_pc36_res1_de_spe_select_highlightexample_4.pdf"), width = 4*5, height = ceiling(length(spatial_p_feature)/4)*7)
ggarrange(plotlist = spatial_p_feature, ncol = 4, nrow = ceiling(length(spatial_p_feature)/4))
dev.off()


########################################################################## Downstream analysis #################################################################
###########################  prepare cropped_stat fro each cluster
source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_hic/E1305-hic-24"
sample <- "E1305-hic-24"
higashi_out <- "8_higashi_out_higashi_500kb"
work_path_seurat <- "10_seurat_bigSpotShape_higashi_500kb"
embed_type <- "scAB_scale"

# hicup: prepare stat_cropped_barcode files for clusters
setwd(paste0(work_path,  "/", work_path_seurat, "/", embed_type))
load(paste0("Out/", sample, "_Seurat_object_anno.RData"))

setwd(work_path)
stat_cropped <- read.table("6_crop/stat_cropped.csv", sep =",", header = TRUE, dec =".", stringsAsFactors = F)
iA_iB <- paste(51-stat_cropped$iA, stat_cropped$iB, sep="x")


for(i in as.numeric(names(table(Seurat_object$seurat_clusters)))){
  stat_cropped_barcode <- stat_cropped[match(colnames(Seurat_object)[Seurat_object$seurat_clusters==i], iA_iB), ]
  stat_cropped_barcode <- data.frame(stat_cropped_barcode$pixel, paste0("E13all.R3.C", i+1, ".hicup"))
  write.table(stat_cropped_barcode, paste0("6_crop/stat_cropped_barcode_C", i+1), sep="\t", row.names=F, col.names=F, quote=F)
}


# hicexplorer: prepare stat_cropped_barcode files for clusters

setwd(work_path)
stat_cropped <- read.table("6_crop/stat_cropped.csv", sep =",", header = TRUE, dec =".", stringsAsFactors = F)
iA_iB <- paste(51-stat_cropped$iA, stat_cropped$iB, sep="x")

for(i in as.numeric(names(table(Seurat_object$seurat_clusters)))){
  stat_cropped_barcode <- stat_cropped[match(colnames(Seurat_object)[Seurat_object$seurat_clusters==i], iA_iB), ]
  
  stat_cropped_barcode_R1 <- data.frame(stat_cropped_barcode$pixel, paste0(sample, "_R1_C", i+1))
  write.table(stat_cropped_barcode_R1, paste0("6_crop/stat_cropped_barcode_R1_C", i+1), sep="\t", row.names=F, col.names=F, quote=F)
  
  stat_cropped_barcode_R2 <- data.frame(stat_cropped_barcode$pixel, paste0(sample, "_R2_C", i+1))
  write.table(stat_cropped_barcode_R2, paste0("6_crop/stat_cropped_barcode_R2_C", i+1), sep="\t", row.names=F, col.names=F, quote=F)
}




########################################################################## DE_ABcom heatmap and function analysis ##########################################################
########################### DE_ABcom heatmap 
# mean of each cell type
dim(Seurat_object@assays$Spatial@counts)
dim(de_markers_24)
count_matrix <- Seurat_object@assays$Spatial@counts

cell_type <- rep(NA, length(Seurat_object$seurat_clusters))
Seurat_object$cell_type <- cell_type
Seurat_object$cell_type <- Seurat_object$seurat_clusters

count_matrix_de_markers <- count_matrix[match(de_markers_24$gene, rownames(count_matrix)),]
zscore_by_row <- function(mat) {
  apply(mat, 1, function(row) {
    (row - mean(row)) / sd(row)
  })
}
count_matrix_de_markers_norm <- zscore_by_row(count_matrix_de_markers)
count_matrix_de_markers_norm <- as.data.frame(count_matrix_de_markers_norm)

count_matrix_de_markers_norm$cell_type <- Seurat_object$cell_type
count_matrix_de_markers_norm_mean <- aggregate(. ~ cell_type, data = count_matrix_de_markers_norm, FUN = mean)
rownames(count_matrix_de_markers_norm_mean) <- count_matrix_de_markers_norm_mean$cell_type
count_matrix_de_markers_norm_mean <- count_matrix_de_markers_norm_mean[,-1]
count_matrix_de_markers_norm_mean_trans <- t(count_matrix_de_markers_norm_mean)

count_matrix_de_markers_norm_mean$cell_type <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9")
count_matrix_de_markers_norm_mean_mean <- aggregate(. ~ cell_type, data = count_matrix_de_markers_norm_mean, FUN = mean)
rownames(count_matrix_de_markers_norm_mean_mean) <- count_matrix_de_markers_norm_mean_mean$cell_type
count_matrix_de_markers_norm_mean_mean <- count_matrix_de_markers_norm_mean_mean[,-1]
count_matrix_de_markers_norm_mean_mean_trans <- t(count_matrix_de_markers_norm_mean_mean)

temp <- as.data.frame(count_matrix_de_markers_norm_mean_mean_trans)
temp <- temp[, c("C2", "C4", "C6", "C8", "C7")]

library(pheatmap)

annotation_col = data.frame(
  cluster = factor(c("C2", "C4", "C6", "C8", "C7"))
)
rownames(annotation_col) = c("C2", "C4", "C6", "C8", "C7")

ann_colors = list(
  cluster = c(C2 = "#377EB8", C4 = "#984EA3", C6 = "#FFFF33", C8 = "#F781BF", C7 = "#A65628")
)

################################ Acom_spe
C2_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==1), which(apply(temp, 1, max) > 0.75))]
C4_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==2), which(apply(temp, 1, max) > 0.75))]
C6_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==3), which(apply(temp, 1, max) > 0.25))]
C8_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==4), which(apply(temp, 1, max) > 0.25))]
C7_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==5), which(apply(temp, 1, max) > 1))]

temp_Acom_spe <- c(C2_spe, C4_spe, C6_spe, C8_spe, C7_spe)
temp_Acom_spe_matrix <- temp[match(temp_Acom_spe, rownames(temp)), ]

length(C2_spe) + length(C4_spe) + length(C6_spe) + length(C8_spe) + length(C7_spe) # 213
length(C2_spe) # 293
length(C4_spe) # 13
length(C6_spe) # 10
length(C7_spe) # 325
length(C8_spe) # 26

pdf(file = paste("Plot/3_allDE_Acom_heat_anno_average.pdf"), width = 4, height = 4)
# pheatmap(temp_Acom_spe_matrix, border = F, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, color = colorRampPalette(c("#C44C9C", "#B96397", "#AA7591", "#968686", "#86917D", "#739D6F", "#5CA65F", "#1CB24B"))(100))
pheatmap(temp_Acom_spe_matrix, border = F, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, color = colorRampPalette(c("#36359D", "#0F74DC", "#05A5C6", "#4EBD90", "#A9BE69", "#F3B947", "#F4E61A"))(100))
dev.off()


annotation_col = data.frame(
  cluster = factor(c("C2C4", "C6", "C8", "C7"))
)
rownames(annotation_col) = c("C2C4", "C6", "C8", "C7")

ann_colors = list(
  cluster = c(C2C4 = "#6766AE", C6 = "#FFFF33", C8 = "#F781BF", C7 = "#A65628")
)

temp_Acom_spe_matrix$C2C4 <- rowMeans(temp_Acom_spe_matrix[, c("C2", "C4")])
temp_Acom_spe_matrix_revision <- temp_Acom_spe_matrix[,c("C2C4", "C6", "C8", "C7")]

pdf(file = paste("Plot/3_allDE_Acom_heat_anno_average_revision.pdf"), width = 4, height = 4)
# pheatmap(temp_Acom_spe_matrix, border = F, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, color = colorRampPalette(c("#C44C9C", "#B96397", "#AA7591", "#968686", "#86917D", "#739D6F", "#5CA65F", "#1CB24B"))(100))
pheatmap(temp_Acom_spe_matrix_revision, border = F, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, color = colorRampPalette(c("#36359D", "#0F74DC", "#05A5C6", "#4EBD90", "#A9BE69", "#F3B947", "#F4E61A"))(100))
dev.off()


### AddModuleScore (Acom)
dim(Seurat_object@assays$Spatial@counts)
count_matrix <- Seurat_object@assays$Spatial@counts

zscore_by_row <- function(mat) {
  apply(mat, 1, function(row) {
    (row - mean(row)) / sd(row)
  })
}
count_matrix_norm <- t(zscore_by_row(count_matrix))
Seurat_object@assays$Spatial@counts <- count_matrix_norm
Seurat_object@assays$Spatial@data <- count_matrix_norm

Seurat_object <- AddModuleScore(Seurat_object,
                                features = list(C2_spe=C2_spe, C6_spe=C6_spe, C7_spe=C7_spe, C8_spe=C8_spe),
                                name = c("C2", "C6", "C7", "C8"),
                                assay = "Spatial")

spatial_p_feature <- list()
for(feature in c("C21", "C62", "C73", "C84")){
  spatial_p_feature[[paste0("umap_",feature)]] <- FeaturePlot(Seurat_object, features = feature, label = FALSE, repel = TRUE) + scale_color_gradientn(colors = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407")) + theme(aspect.ratio = 1)
  spatial_p_feature[[paste0("spatial_",feature)]] <- SpatialFeaturePlot(object = Seurat_object, features = feature, alpha = c(1, 1), ncol = 1, pt.size.factor = 1.6, stroke = NA) + scale_fill_gradientn(colors = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), na.value = "white") + theme(aspect.ratio = 1) #  + theme(legend.position="right")
  spatial_p_feature[[paste0("spatial_",feature)]]$layers[[1]]$aes_params=c(spatial_p_feature[[paste0("spatial_",feature)]]$layers[[1]]$aes_params, shape=22)
} 
pdf(file = paste("Plot/3_allDE_Acom_FeaturePlot_anno_average.pdf"), width = 4*7, height = 2*7)
p <- ggarrange(plotlist = spatial_p_feature, ncol = 4, nrow = 2)
print(p)
dev.off()

###########################  PCGs located in DE_Acom regions
C2_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==1), which(apply(temp, 1, max) > 0.75))]
C4_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==2), which(apply(temp, 1, max) > 0.75))]
C6_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==3), which(apply(temp, 1, max) > 0.25))]
C8_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==4), which(apply(temp, 1, max) > 0.25))]
C7_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==5), which(apply(temp, 1, max) > 1.95))]

length(C2_spe) + length(C4_spe) + length(C6_spe) + length(C8_spe) + length(C7_spe) # 213
length(C2_spe) # 293
length(C4_spe) # 13
length(C6_spe) # 10
length(C7_spe) # 325
length(C8_spe) # 26

C2_spe_matrix <- data.frame(C2_spe, "C2"); colnames(C2_spe_matrix) <- c("Acom", "cell_type")
C4_spe_matrix <- data.frame(C4_spe, "C4"); colnames(C4_spe_matrix) <- c("Acom", "cell_type")
C6_spe_matrix <- data.frame(C6_spe, "C6"); colnames(C6_spe_matrix) <- c("Acom", "cell_type")
C7_spe_matrix <- data.frame(C7_spe, "C7"); colnames(C7_spe_matrix) <- c("Acom", "cell_type")
C8_spe_matrix <- data.frame(C8_spe, "C8"); colnames(C8_spe_matrix) <- c("Acom", "cell_type")

Cluster_spe_matrix <- rbind(C2_spe_matrix, C4_spe_matrix, C6_spe_matrix, C8_spe_matrix, C7_spe_matrix)

DE_ABcom_chr <- unlist(strsplit(Cluster_spe_matrix$Acom, ":"))[seq(from = 1, to = length(unlist(strsplit(Cluster_spe_matrix$Acom, ":"))), by = 2)]
DE_ABcom_pos <- unlist(strsplit(Cluster_spe_matrix$Acom, ":"))[seq(from = 2, to = length(unlist(strsplit(Cluster_spe_matrix$Acom, ":"))), by = 2)]
DE_ABcom_start <- unlist(strsplit(DE_ABcom_pos, "-"))[seq(from = 1, to = length(unlist(strsplit(DE_ABcom_pos, "-"))), by = 2)]
DE_ABcom_end <- unlist(strsplit(DE_ABcom_pos, "-"))[seq(from = 2, to = length(unlist(strsplit(Cluster_spe_matrix$Acom, "-"))), by = 2)]
DE_ABcom_region_Acom <- data.frame(DE_ABcom_chr, DE_ABcom_start, DE_ABcom_end, Cluster_spe_matrix)

write.table(DE_ABcom_region_Acom, "Plot/DE_ABcom_region_Acom.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

### command run
Reference=/media/maoni/data/Reference/mouse/GRCm38_mm10
AB_path=/media/maoni/data/CZP/spatial_hic/E1305-hic-24/10_seurat_bigSpotShape_higashi_500kb/scAB_scale/Plot
intersectBed -a $AB_path/DE_ABcom_region_Acom.bed -b $Reference/ensembl/All_gene_Anno_PCG.bed -wa -wb > $AB_path/DE_ABcom_region_Acom_PCGs.bed

DE_ABcom_region_Acom_PCGs <- read.table("Plot/DE_ABcom_region_Acom_PCGs.bed", sep = "\t", header = FALSE)
table(DE_ABcom_region_Acom_PCGs$V5)

# C2  C4  C6  C7  C8 
# 420  19  29 203 188

DE_ABcom_region_Acom_PCGs_unique <- unique(DE_ABcom_region_Acom_PCGs$V9)
length(DE_ABcom_region_Acom_PCGs_unique)

##########################  load RNA data (PCGs ggvioplot located in DE_Acom regions)
source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_transcriptome/E13-23-RNA"
sample <- "E13-23-RNA"
setwd(paste0(work_path, "/6_seurat_bigSpotShape_new"))
PC=30
res=0.8
load(paste0("Out/", sample, "_cropped_PC",PC,"_res",res,"_anno.RData"))

C2_spe_PCGs <- unique(DE_ABcom_region_Acom_PCGs$V9[DE_ABcom_region_Acom_PCGs$V5=="C2"])
C4_spe_PCGs <- unique(DE_ABcom_region_Acom_PCGs$V9[DE_ABcom_region_Acom_PCGs$V5=="C4"])
C6_spe_PCGs <- unique(DE_ABcom_region_Acom_PCGs$V9[DE_ABcom_region_Acom_PCGs$V5=="C6"])
C8_spe_PCGs <- unique(DE_ABcom_region_Acom_PCGs$V9[DE_ABcom_region_Acom_PCGs$V5=="C8"])
C7_spe_PCGs <- unique(DE_ABcom_region_Acom_PCGs$V9[DE_ABcom_region_Acom_PCGs$V5=="C7"])

Seurat_object_cropped <- AddModuleScore(Seurat_object_cropped,
                                features = list(C2_spe_PCGs=C2_spe_PCGs, C4_spe_PCGs=C4_spe_PCGs, C6_spe_PCGs=C6_spe_PCGs, C8_spe_PCGs=C8_spe_PCGs, C7_spe_PCGs=C7_spe_PCGs),
                                name = c("C2", "C4", "C6", "C8", "C7"),
                                assay = "SCT")
head(Seurat_object_cropped@meta.data)
DE_ABcom_region_Acom_PCGs_meanExp <- Seurat_object_cropped@meta.data[grep(pattern = "Inhibitory neuron progenitors|Inhibitory interneurons|Excitatory/Inhibitory neurons|Cholinergic neurons|Radial glia|Premature oligodendrocyte|Cardiac muscle lineages|Hepatocytes & Primitive erythroid lineage" ,Seurat_object_cropped$cell_type), c("C21", "C42", "C63", "C84", "C75", "cell_type")]
head(DE_ABcom_region_Acom_PCGs_meanExp)
DE_ABcom_region_Acom_PCGs_meanExp <- DE_ABcom_region_Acom_PCGs_meanExp[,c(1,3,4,5,6)]

# ggviolin plots with p-values
library(rstatix)
source("/media/maoni/data/Zhanglin_Junjie/ColorPalettes.R")
my_color_palette <- ArchRPalettes$bear
names(my_color_palette) <- c("Inhibitory neuron progenitors", "Inhibitory interneurons", "Excitatory/Inhibitory neurons", "Cholinergic neurons", "Radial glia", "Premature oligodendrocyte",
                             "Cardiac muscle lineages", "Hepatocytes & Primitive erythroid lineage", 
                             "Connective tissue progenitors", "Ependymal cell", "Chondrocytes & osteoblasts", "Epithelial cells", "Myocytes", "Osteoblasts", "Endothelial cells", "Other")


DE_ABcom_region_Acom_PCGs_meanExp_2 <- DE_ABcom_region_Acom_PCGs_meanExp
DE_ABcom_region_Acom_PCGs_meanExp_2$cell_type <- as.character(DE_ABcom_region_Acom_PCGs_meanExp_2$cell_type)
DE_ABcom_region_Acom_PCGs_meanExp_2[grep("Inhibitory neuron progenitors|Inhibitory interneurons|Excitatory/Inhibitory neurons|Cholinergic neurons", DE_ABcom_region_Acom_PCGs_meanExp_2$cell_type),5] <- "CNS"
DE_ABcom_region_Acom_PCGs_meanExp_2[grep("Radial glia|Premature oligodendrocyte", DE_ABcom_region_Acom_PCGs_meanExp_2$cell_type),5] <- "Gila"
DE_ABcom_region_Acom_PCGs_meanExp_2[grep("Cardiac muscle lineages", DE_ABcom_region_Acom_PCGs_meanExp_2$cell_type),5] <- "Heart"
DE_ABcom_region_Acom_PCGs_meanExp_2[grep("Hepatocytes & Primitive erythroid lineage", DE_ABcom_region_Acom_PCGs_meanExp_2$cell_type),5] <- "Liver"
DE_ABcom_region_Acom_PCGs_meanExp_2$cell_type <- factor(DE_ABcom_region_Acom_PCGs_meanExp_2$cell_type, levels = c("CNS", "Gila", "Heart", "Liver"))

cluster <- c("C24", "C6", "C8", "C7")
cell_type <- c("CNS", "Gila", "Heart", "Liver")
p_all <-list()
for(i in 1:4){
  # i <- 4
  Cluster_spe_PCGs_meanExp <- data.frame(cluster[i], DE_ABcom_region_Acom_PCGs_meanExp_2[,i], DE_ABcom_region_Acom_PCGs_meanExp_2$cell_type)
  colnames(Cluster_spe_PCGs_meanExp) <- c("group", "Acom_PCGs_meanExp", "cell_type")
  
  stat.test <- Cluster_spe_PCGs_meanExp %>%
    group_by(group) %>%
    t_test(Acom_PCGs_meanExp ~ cell_type, ref.group = cell_type[i])
  stat.test
  
  stat.test <- stat.test %>% add_y_position()
  
  p_all[[i]] <- ggviolin(Cluster_spe_PCGs_meanExp, x = "cell_type", y = "Acom_PCGs_meanExp", fill = "cell_type", facet.by = "group", 
                            ncol = 1, font.label = list(size = 3, color = "black"),
                            add = "boxplot",
                            add.params = list(fill = "white", width = 0.05,linetype = 1)) +
    theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
    # scale_fill_manual(values=my_color_palette)
}

pdf("Plot/3_DE_ABcom_region_Acom_PCGs_vlnplot_2.pdf", width = 12, height= 5)
ggarrange(p_all[[1]], p_all[[2]], p_all[[3]], p_all[[4]], ncol = 4, nrow = 1)
dev.off()


################# DE_Acom function analysis
library(clusterProfiler)
library(org.Mm.eg.db)
library(pathview)
library(enrichplot)
library(ggplot2)
library(ggpubr)
library(stringi)

source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_hic/E1305-hic-24"
sample <- "E1305-hic-24"
higashi_out <- "8_higashi_out_higashi_500kb"
work_path_seurat <- "10_seurat_bigSpotShape_higashi_500kb"
embed_type <- "scAB_scale"

setwd(paste0(work_path,  "/", work_path_seurat, "/", embed_type))
DE_ABcom_region_Acom_PCGs <- read.table("Plot/DE_ABcom_region_Acom_PCGs.bed", sep = "\t", header = FALSE)
table(DE_ABcom_region_Acom_PCGs$V5)

PCG_list <- list()
for (organ in c("C2|C4","C6","C7", "C8")){
  # organ <- "C2|C4"
  gene_list_PCG <- unique(DE_ABcom_region_Acom_PCGs$V9[grep(organ, DE_ABcom_region_Acom_PCGs$V5)])
  PCG_list[[organ]] <- gene_list_PCG
  
  gene_list <- sort(gene_list_PCG, decreasing = TRUE)
  gene_id_list <- bitr(gene_list,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Mm.eg.db)
  gene_id_list <- dplyr::distinct(gene_id_list,SYMBOL,.keep_all=TRUE) # remove duplication
  gene_id_list_sort <- sort(gene_id_list$ENTREZID, decreasing = TRUE)
  
  go_enrich_BP <- enrichGO(gene = gene_id_list_sort,
                           OrgDb = org.Mm.eg.db,
                           keyType = 'ENTREZID',
                           readable = T,
                           ont = "BP",
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH")
  go_enrich_BP <- as.data.frame(go_enrich_BP)
  write.table(go_enrich_BP,paste0(organ, "_DE_ABcom_region_Acom_PCGs_enrichGO_GOBP.txt"),sep="\t",col.names=T,row.names=F,quote=F)
  
  go_enrich_KEGG <- enrichKEGG(gene = gene_id_list_sort,
                               organism = "hsa",
                               keyType = "kegg",
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH")
  go_enrich_KEGG <- as.data.frame(go_enrich_KEGG)
  write.table(go_enrich_KEGG,paste0(organ, "_DE_ABcom_region_Acom_PCGs_enrichKEGG_KEGG.txt"),sep="\t",col.names=T,row.names=F,quote=F)
}


### barplot for select GO terms
library(ggplot2)
library(forcats)

setwd(paste0(work_path, "/6_seurat_bigSpotShape_new"))

my_color <- c("#6766AD", brewer.pal(9, "Set1")[c(6,8,7)])
names(my_color) <- c("C2|C4", "C6", "C8", "C7")

GO_barplot_all <- c()
GO_barplot_list <- list()
for (organ in c("CNS","Radial_glia","Heart", "Liver")){
  # organ <- "Radial_glia"
  go_enrich_BP <- read.csv(paste0(organ, "_DE_ABcom_allDEG_enrichGO_GOBP.txt"), sep="\t", header=T)
  GO_barplot <- data.frame(go_enrich_BP[which(go_enrich_BP$highlight==1), ][, c("Description", "pvalue", "Count")], organ)[1:5,]
  GO_barplot$pvalue <- -log10(GO_barplot$pvalue)
  colnames(GO_barplot) <- c("GO_term", "pvalue", "Count", "organ")
  GO_barplot$GO_term <- factor(GO_barplot$GO_term, levels = rev(GO_barplot$GO_term))
  GO_barplot_all <- rbind(GO_barplot_all, GO_barplot)
}

GO_barplot_all$organ <- factor(GO_barplot_all$organ, levels = c("CNS","Radial_glia","Heart", "Liver"))
p <- ggplot(GO_barplot_all, aes(x=GO_term,y=pvalue, fill=organ)) + geom_bar(stat="identity", width=0.8) + coord_flip() +  xlab("GOBP") + ylab("-log10(pvalue)") +
  facet_wrap(~organ, scales = "free_y") +
  theme_bw() + 
  theme(legend.position="none") + 
  theme(panel.grid = element_blank()) + 
  scale_fill_manual(values=my_color)

pdf("Plot/DE_ABcom_region_Acom_PCGs_enrichGO_GOBP_barplot_all.pdf", width = 9, height = 4)
print(p)
dev.off()








###################################################  Identify DEGs located in DE ABcom regions ################################################### 
###########################  load Hi-C data
source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_hic/E1305-hic-24"
sample <- "E1305-hic-24"
higashi_out <- "8_higashi_out_higashi_500kb"
work_path_seurat <- "10_seurat_bigSpotShape_higashi_500kb"
embed_type <- "scAB_scale"

setwd(paste0(work_path,  "/", work_path_seurat, "/", embed_type))
load(paste0("Out/", sample, "_Seurat_object.RData"))
# top10_hic <- de_markers_24 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)  %>% arrange(cluster, desc(avg_log2FC))
DE_ABcom <- de_markers_24$gene[which(de_markers_24$cluster==1 | de_markers_24$cluster==3 | de_markers_24$cluster==5 | de_markers_24$cluster==6 | de_markers_24$cluster==7)]

DE_ABcom_chr <- unlist(strsplit(DE_ABcom, ":"))[seq(from = 1, to = length(unlist(strsplit(DE_ABcom, ":"))), by = 2)]
DE_ABcom_pos <- unlist(strsplit(DE_ABcom, ":"))[seq(from = 2, to = length(unlist(strsplit(DE_ABcom, ":"))), by = 2)]
DE_ABcom_start <- unlist(strsplit(DE_ABcom_pos, "-"))[seq(from = 1, to = length(unlist(strsplit(DE_ABcom_pos, "-"))), by = 2)]
DE_ABcom_end <- unlist(strsplit(DE_ABcom_pos, "-"))[seq(from = 2, to = length(unlist(strsplit(DE_ABcom_pos, "-"))), by = 2)]
DE_ABcom_region <- data.frame(DE_ABcom_chr, DE_ABcom_start, DE_ABcom_end, DE_ABcom)
write.table(DE_ABcom_region, "Plot/DE_ABcom_region.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

### command run
Reference=/media/maoni/data/Reference/mouse/GRCm38_mm10
AB_path=/media/maoni/data/CZP/spatial_hic/E1305-hic-24/10_seurat_bigSpotShape_higashi_500kb/scAB_scale/Plot
intersectBed -a $AB_path/DE_ABcom_region.bed -b $Reference/ensembl/All_gene_Anno_PCG.bed -wa -wb > $AB_path/DE_ABcom_region_PCGs.bed

DE_ABcom_region_PCGs <- read.table("Plot/DE_ABcom_region_PCGs.bed", sep = "\t", header = FALSE)
DE_ABcom_region_PCGs_unique <- unique(DE_ABcom_region_PCGs$V8)


################################# load RNA data
source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_transcriptome/E13-23-RNA"
sample <- "E13-23-RNA"
setwd(paste0(work_path, "/6_seurat_bigSpotShape_new"))
PC=30
res=0.8
load(paste0("Out/", sample, "_cropped_PC",PC,"_res",res,"_anno.RData"))

top50_anno <- de_markers %>% group_by(cell_type) %>% top_n(n = 50, wt = avg_log2FC)  %>% arrange(cell_type, desc(avg_log2FC))
all_anno <- de_markers

### DEG DE_ABcom overlap
DE_ABcom_allDEGs <- data.frame(all_anno[which(all_anno$gene %in% DE_ABcom_region_PCGs_unique),])
DE_ABcom_allDEGs_final <- c()
for(i in 1:length(DE_ABcom_allDEGs$gene)){
  ABcom <- paste(DE_ABcom_region_PCGs$V4[which(DE_ABcom_region_PCGs$V8 %in% DE_ABcom_allDEGs$gene[i])], collapse=";")
  tmp <- cbind(DE_ABcom_allDEGs[i,], ABcom)
  DE_ABcom_allDEGs_final <- rbind(DE_ABcom_allDEGs_final, tmp)
}
write.csv(DE_ABcom_allDEGs_final, "Plot/DE_ABcom_allDEGs_final.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

pdf(paste0("Plot/3_findmarker_dotplot_data_PC",PC,"_res",res,"_DE_ABcom_allDEGs_final.pdf"), height = 40, width = 16)
p <- DotPlot(Seurat_object_cropped, features = unique(DE_ABcom_allDEGs_final$gene), cols = c("grey", "red"), group.by = 'cell_type') + 
  coord_flip() + 
  theme(axis.text.x = element_text(hjust = 1, angle = 45), axis.title = element_blank()) +
  ggtitle(label = 'DE_ABcom_allDEGs_final')
print(p)
dev.off()


################################# DEG DE_ABcom highlight
# hic magic
source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_hic/E1305-hic-24"
sample <- "E1305-hic-24"
higashi_out <- "8_higashi_out_higashi_500kb"
work_path_seurat <- "10_seurat_bigSpotShape_higashi_500kb"
embed_type <- "scAB_scale"
setwd(paste0(work_path,  "/", work_path_seurat, "/", embed_type))
load(paste0("Out/", sample, "_Seurat_object_magic.RData"))


# RNA magic
source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_transcriptome/E13-23-RNA"
sample <- "E13-23-RNA"
setwd(paste0(work_path, "/6_seurat_bigSpotShape_new"))
Seurat_object_cropped_magic <- readRDS(paste0("Out/", sample, "_cropped_magic", ".rds"))

# highlight

for(n in 1:ceiling(length(DE_ABcom_allDEGs_final$cell_type)/100)){
  
  # i <- match("Slc17a6", DE_ABcom_allDEGs_final$gene)
  
  DE_ABcom_allDEGs_final_subset <- DE_ABcom_allDEGs_final[i,]
  
  
  spatial_p_feature <- list()
  for(i in 1:length(DE_ABcom_allDEGs_final_subset$gene)){
    DEG_spe_select <- DE_ABcom_allDEGs_final_subset$gene[i]
    DEAB_spe_select <- unlist(strsplit(DE_ABcom_allDEGs_final_subset$ABcom[i], ";"))
    
    spatial_p_feature[[DEG_spe_select]] <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = DEG_spe_select, alpha = c(1, 1), ncol = 1, pt.size.factor = 1.6, stroke = NA) + theme(aspect.ratio = 1) + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) + theme(aspect.ratio = 1) #  + theme(legend.position="right")
    spatial_p_feature[[DEG_spe_select]]$layers[[1]]$aes_params=c(spatial_p_feature[[DEG_spe_select]]$layers[[1]]$aes_params, shape=22)
    
    for(j in 1:length(DEAB_spe_select)){
      spatial_p_feature[[DEAB_spe_select[j]]] <- SpatialFeaturePlot(object = Seurat_object_magic, features = DEAB_spe_select[j], alpha = c(1, 1), ncol = 1, pt.size.factor = 1.6, stroke = NA) + scale_fill_gradientn(values=c(0,0.3,0.5,0.6,0.7,0.8,1),colors = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), na.value = "white") + theme(aspect.ratio = 1) # + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) #  + theme(legend.position="right")
      spatial_p_feature[[DEAB_spe_select[j]]]$layers[[1]]$aes_params=c(spatial_p_feature[[DEAB_spe_select[j]]]$layers[[1]]$aes_params, shape=22)
    }
  }
  
  png(file = paste0("Plot/4_findmarker_2D_impose.data_pc36_res1_de_spe_select_RNA_HT", ".png"), width = 4*480, height = ceiling(length(spatial_p_feature)/4)*480)
  p <- ggarrange(plotlist = spatial_p_feature, ncol = 4, nrow = ceiling(length(spatial_p_feature)/4))
  print(p)
  dev.off()
}


DE_ABcom_allDEGs_final_subset <- DE_ABcom_allDEGs_final[i,]

DEG_highlight_example <- c("Sox2", "Nrxn1", "Myh6", "Afp")
spatial_p_feature <- list()
for(i in 1:length(DEG_highlight_example)){
  DEG_spe_select <- DEG_highlight_example[i]
  spatial_p_feature[[DEG_spe_select]] <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = DEG_spe_select, alpha = c(1, 1), ncol = 1, pt.size.factor = 1.6, stroke = NA) + theme(aspect.ratio = 1) + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) + theme(aspect.ratio = 1) #  + theme(legend.position="right")
  spatial_p_feature[[DEG_spe_select]]$layers[[1]]$aes_params=c(spatial_p_feature[[DEG_spe_select]]$layers[[1]]$aes_params, shape=22)
}
  
pdf(file = paste0("Plot/5_Spatial_PC30_res0.8_final_gene_highlight_figS2", ".pdf"), width = 4*7, height = ceiling(length(spatial_p_feature)/4)*7)
p <- ggarrange(plotlist = spatial_p_feature, ncol = 4, nrow = ceiling(length(spatial_p_feature)/4))
print(p)
dev.off()


#### revision (highlight synapse assembly genes)
# cd /media/maoni/data/Reference/GOA
# wget http://current.geneontology.org/annotations/mgi.gaf.gz
# gunzip mgi.gaf.gz
# grep "GO:0007416" mgi.gaf | cut -f2 | sort | uniq > genes_GO_0007416_mouse.txt


synapse_assembly_genes_all <- read.table("/media/maoni/data/Reference/GOA/genes_GO_0007416_mouse.txt")[,1]
synapse_assembly_genes_enriched <- c("Lrrtm1", "Vstm5", "Magi2", "Syndig1", "Negr1", "Lrrc4b", "Clstn3", "Lrrtm3", "Lzts1", "Adgrb3", "Reln", "Ntrk3", 
                                     "Nrxn3", "Nrg3", "Lrrn3", "Lrrn1", "Large1", "Gria1", "Gap43", "Erbb4")
synapse_assembly_genes_example <- c("Nrxn3", "Nrg3")

Seurat_object_cropped_magic <- AddModuleScore(Seurat_object_cropped_magic,
                                features = list(All=synapse_assembly_genes_all, Enriched=synapse_assembly_genes_enriched, Nrxn3=synapse_assembly_genes_example[1], Nrg3=synapse_assembly_genes_example[2]),
                                name = c("All", "Enriched", "Nrxn3", "Nrg3"),
                                assay = "MAGIC_SCT")

featurelist <- c("All1", "Enriched2", "Nrxn33", "Nrg34")
spatial_p_feature <- list()
for(i in 1:4){
  spatial_p_feature[[featurelist[i]]] <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = featurelist[i], alpha = c(1, 1), ncol = 1, pt.size.factor = 1.6, stroke = NA) + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))[c(1:11)]) + theme(aspect.ratio = 1) #  + theme(legend.position="right")
  spatial_p_feature[[featurelist[i]]]$layers[[1]]$aes_params=c(spatial_p_feature[[featurelist[i]]]$layers[[1]]$aes_params, shape=22)
}

pdf("/media/maoni/data/CZP/Figures/revision/Stero-seq/synapse_assembly_genes_highlight_2.pdf", width = 4*4, height = ceiling(length(spatial_p_feature)/4)*4)
p <- ggarrange(plotlist = spatial_p_feature, ncol = 4, nrow = ceiling(length(spatial_p_feature)/4))
print(p)
dev.off()


###################################### ABcom analysis example (Figure 2)
# cluster highlight
source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_hic/E1305-hic-24"
sample <- "E1305-hic-24"
higashi_out <- "8_higashi_out_higashi_500kb"
work_path_seurat <- "10_seurat_bigSpotShape_higashi_500kb"
embed_type <- "scAB_scale"
setwd(paste0(work_path,  "/", work_path_seurat, "/", embed_type))

load(paste0("Out/", sample, "_Seurat_object.RData"))

cluster <- c("1|3", "5", "7", "6")
celltype <- c("CNS", "Radial glia", "HT", "LV")
spatial_p_cluster <- list()
for(i in 1:length(cluster)){
  Seurat_object_subset <- Seurat_object[ ,grep(cluster[i], Seurat_object$seurat_clusters)]
  spatial_p_cluster[[celltype[i]]] <- SpatialDimPlot(Seurat_object_subset, group.by = "seurat_clusters", crop = FALSE, label = FALSE, label.size = 3, pt.size.factor = 1.6, stroke = NA,  label.box = FALSE, cols = idents_cols_raw) + theme(aspect.ratio = 1) + ggtitle(celltype[i]) + theme(legend.position="none") + theme(plot.title = element_text(size = 10)) # the spatial plot
  spatial_p_cluster[[celltype[i]]]$layers[[1]]$aes_params=c(spatial_p_cluster[[celltype[i]]]$layers[[1]]$aes_params, shape=22)
}

pdf(file = paste0("Plot/5_Spatial_clusters_PC",PC,"_res",res,"_final_C2C4C6C7C8.pdf"), width = 3, height = 10)
p_spatial <- spatial_p_cluster[[celltype[1]]]/spatial_p_cluster[[celltype[2]]]/spatial_p_cluster[[celltype[3]]]/spatial_p_cluster[[celltype[4]]]
print(p_spatial)
dev.off()

# overlap with markers in paper
gene_list_paper <- c("Il1rapl2", "Meox2", "Tgfb2", "Adamts9", "Postn", "Ror1", "Runx2", "Twist2", "Prrx1", "Wt1", "Mylk", "Ednra", "Sox9", 
                     "Foxp2", "Col2a1", "Col9a1", "Col11a1", "Pax9", "Ntng1", "Car10", "Epcam", "Trp63", "Grhl2", "Pth2r", "Fabp7", "Pax3", 
                     "Fzd10", "Hes5", "Gpc5", "Smoc1", "Prmt8", "Gadd45g", "Cdkn1c", "Btg2", "Nkx6-3", "Nrn1", "Slc17a6", "Grem2", "Slc17a6", 
                     "Sox1", "Olig2", "Nkx2-1", "Fgf8", "En2", "Fgf15", "Fgf17", "Pax5", "Neb", "Myh3", "Tpm2", "Acta2", "Foxb1", "Scube2", "Prtg", 
                     "Pax2", "Slc6a5", "Il23a", "Bmpr1a", "Prtg", "Col1a1", "Camk1d", "Rbm8a", "Pax2", "Slc6a5", "Fut9", "Id4", "Pcdh19", "Cdon", 
                     "Emx1", "Ptprb", "Pecam1", "Vwf", "Klhl4", "Hbegf", "Egfl7", "ITGA11", "Atp1a2", "Lamc3", "Epha7", "Snca", 
                     "Slc4a1", "Kel", "Plp1", "Cdh19", "Syt13", "Shox2", "Ptprr", "Pcbp3", "Msx1", "Fgf10", "Wnt5a", "Lmx1b", 
                     "Hba-a1", "Hbb-y", "Hba-x", "Hbb-bh1", "Hbb-bs", "Abcb4", "Dlx1", "Dlx2", "Neurod2", "Tiam2", "Afp", "Alb", "Apoa2", "Afp29", "Pik3c2g", "Hoga1", 
                     "TFEB", "Shh", "Slit1", "Ntn1", "Apoe", "Lyz2", "Selenop", "Ptprc", "Ly86", "Ctss", "Sostdc1", "Htr2c", "Kcne2", "Ttr", "Slit2", 
                     "Slit3", "Chat", "Myl2", "Myocd", "Hcn4", "Ctnna3", "Ryr2", "Tbx20", "Pf4", "Itgb3", "Itga2b", "Ppbp", "Cd226", "Tyr", "Trpm1 ", 
                     "Pmel", "GJA8", "Cryba1", "Cryaa", "Ngp", "S100a8")

DE_ABcom_allDEGs_final <- read.csv("Plot/DE_ABcom_allDEGs_final.csv", sep=",", header=TRUE)
Found_in_paper <- DE_ABcom_allDEGs_final[which(DE_ABcom_allDEGs_final$gene %in% gene_list_paper), ]
write.csv(Found_in_paper, "Plot/Found_in_paper.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

## Dotplot example
source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_transcriptome/E13-23-RNA"
sample <- "E13-23-RNA"
setwd(paste0(work_path, "/6_seurat_bigSpotShape_new"))
PC=30
res=0.8
load(paste0("Out/", sample, "_cropped_PC",PC,"_res",res,"_anno.RData"))
DE_ABcom_allDEGs_final <- read.csv("Plot/DE_ABcom_allDEGs_final.csv", sep=",", header=TRUE)

gene_list_example <- list("CNS"=c("Slc17a6", "Nrxn3", "Syt1", "Car10","Rbfox1"),
                          "Radial glia"=c("Nkx2-1", "Ascl1","Rgs20"),
                          "Heart"=c("Ryr2", "Ctnna3", "Tbx20", "Csrp3"),
                          "Liver"=c("Afp", "Kel", "Apoe", "Afm"))

Seurat_object_cropped_subset <- Seurat_object_cropped[, grep(pattern = "Inhibitory neuron progenitors|Inhibitory interneurons|Excitatory/Inhibitory neurons|Cholinergic neurons|Radial glia|Premature oligodendrocyte|Cardiac muscle lineages|Hepatocytes & Primitive erythroid lineage" ,Seurat_object_cropped$cell_type)]

p_list <- list()
for(i in 1:length(gene_list_example)){
  p_list[[i]] <- DotPlot(Seurat_object_cropped_subset, features = rev(gene_list_example[[i]]), cols = c("grey95", "#B30000"), group.by = 'cell_type', scale = FALSE) + 
                         coord_flip() + 
                         theme(axis.text.x = element_blank(), axis.title.y = element_text(size = 10), axis.title = element_blank()) +
                         # theme(axis.text.x = element_text(hjust = 1, angle = 45), axis.title = element_blank()) +
                         theme(legend.position = "right", legend.key.size = unit(0.3, "cm")) + 
                         theme(panel.border = element_rect(color = "black", size = 0.5), legend.text = element_text(size = 10), legend.title = element_text(size = 10))
}

pdf(file = paste0("Plot/5_Spatial_DotPlot_PC",PC,"_res",res,"_final_C2C4C6C7C8.pdf"), width = 7, height = 10)
p_dotplot <- p_list[[1]]/p_list[[2]]/p_list[[3]]/p_list[[4]]
print(p_dotplot)
dev.off()

source("/media/maoni/data/Zhanglin_Junjie/ColorPalettes.R") # cols = paletteContinuous(set = "comet") (set Continuous colors)
my_color_palette <- ArchRPalettes$bear[1:8]
names(my_color_palette) <- c("Inhibitory neuron progenitors", "Inhibitory interneurons", "Excitatory/Inhibitory neurons", "Cholinergic neurons", "Radial glia", "Premature oligodendrocyte",
                             "Cardiac muscle lineages", "Hepatocytes & Primitive erythroid lineage")

pdf(file = paste0("Plot/5_Spatial_DotPlot_PC",PC,"_res",res,"_final_C2C4C6C7C8_color.pdf"), width = 10, height = 1)
scales::show_col(my_color_palette, labels = FALSE, borders = NA, cex_label = 1, ncol = length(my_color_palette))
dev.off()


# DEG DE_ABcom example hightlight
source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_transcriptome/E13-23-RNA"
sample <- "E13-23-RNA"
PC=30
res=0.8
setwd(paste0(work_path, "/6_seurat_bigSpotShape_new"))

DE_ABcom_allDEGs_final <- read.csv("Plot/DE_ABcom_allDEGs_final.csv", sep=",", header=TRUE)

gene_example_final <- c("Nrxn3", "Ascl1", "Ryr2", "Kel")

DE_ABcom_allDEGs_final_subset <- DE_ABcom_allDEGs_final[match(gene_example_final, DE_ABcom_allDEGs_final$gene), ]

spatial_p_feature <- list()
for(i in 1:length(DE_ABcom_allDEGs_final_subset$gene)){
  DEG_spe_select <- DE_ABcom_allDEGs_final_subset$gene[i]
  DEAB_spe_select <- unlist(strsplit(DE_ABcom_allDEGs_final_subset$ABcom[i], ";"))[[1]]
  
  spatial_p_feature[[DEG_spe_select]] <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = DEG_spe_select, alpha = c(1, 1), ncol = 1, pt.size.factor = 1.6, stroke = NA) + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) + theme(aspect.ratio = 1)
  spatial_p_feature[[DEG_spe_select]]$layers[[1]]$aes_params=c(spatial_p_feature[[DEG_spe_select]]$layers[[1]]$aes_params, shape=22)
  
  spatial_p_feature[[DEAB_spe_select]] <- SpatialFeaturePlot(object = Seurat_object_magic, features = DEAB_spe_select, alpha = c(1, 1), ncol = 1, pt.size.factor = 1.6, stroke = NA) + scale_fill_gradientn(values=c(0,0.3,0.5,0.6,0.7,0.8,1),colors = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), na.value = "white") + theme(aspect.ratio = 1)
  spatial_p_feature[[DEAB_spe_select]]$layers[[1]]$aes_params=c(spatial_p_feature[[DEAB_spe_select]]$layers[[1]]$aes_params, shape=22)
}

pdf(file = paste0("Plot/5_Spatial_PC",PC,"_res",res,"_final_C2C4C6C7C8_DEG_DEAB_example_final.pdf"), width = 2*7, height = 4*7)
p <- ggarrange(plotlist = spatial_p_feature, ncol = 2, nrow = 4)
print(p)
dev.off()


# homer ABcom track seperate
setwd("/media/maoni/data/CZP/spatial_hic/hicup/E13/E13all/R3/homer/")

for(class in c("C2", "C4", "C6", "C7", "C8")){
  ABcom <- read.table(paste0("R3.100k.400k.", class, ".CompartmentHomer.PC1.Final.bedGraph"), sep="\t", header = TRUE)
  Acom <- ABcom[ABcom$E1 > 0, ]
  Bcom <- ABcom[ABcom$E1 < 0, ]
  write.table(Acom, paste0("R3.100k.400k.", class, ".CompartmentHomer.PC1.Final.Acom.bedGraph"), sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(Bcom, paste0("R3.100k.400k.", class, ".CompartmentHomer.PC1.Final.Bcom.bedGraph"), sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}


# homer ABcom track seperate(Bulk)
setwd("/media/maoni/data/CZP/bulk_hic/mouseOrgan/HOMER/PC1/Comb")

for(sm in c("E13FB", "E13HB", "E13SP", "E13HT", "E13LV")){
  ABcom <- read.table(paste0(sm, ".100k.400k.Comb.CompartmentHomer.PC1.final.bedGraph"), sep="\t", header = TRUE)
  Acom <- ABcom[ABcom$E1 > 0, ]
  Bcom <- ABcom[ABcom$E1 < 0, ]
  write.table(Acom, paste0(sm, ".100k.400k.Comb.CompartmentHomer.PC1.final.Acom.bedGraph"), sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(Bcom, paste0(sm, ".100k.400k.Comb.CompartmentHomer.PC1.final.Bcom.bedGraph"), sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}

############################################# Function analysis
# BiocManager::install("clusterProfiler",force = TRUE,version=4)
# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("pathview")
# devtools::install_github("datapplab/pathview")
# BiocManager::install("enrichplot")
# install.packages("/media/maoni/data/Rpackages/org.Mm.eg.db_3.18.0.tar.gz", repos = NULL, type = "source", dependencies = TRUE)
# install.packages("/media/maoni/data/Rpackages/org.Hs.eg.db_3.18.0.tar.gz", repos = NULL, type = "source", dependencies = TRUE)

library(clusterProfiler)
library(org.Mm.eg.db)
library(pathview)
library(enrichplot)
library(ggplot2)
library(ggpubr)
library(stringi)

source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_transcriptome/E13-23-RNA"
sample <- "E13-23-RNA"
setwd(paste0(work_path, "/6_seurat_bigSpotShape_new"))

DE_ABcom_allDEGs_final <- read.csv("Plot/DE_ABcom_allDEGs_final.csv", sep=",", header=TRUE)

DE_ABcom_allDEGs_final$class <- rep(NA, length(DE_ABcom_allDEGs_final$cell_type))
DE_ABcom_allDEGs_final$class[grep("Inhibitory neuron progenitors|Inhibitory interneurons|Excitatory/Inhibitory neurons|Cholinergic neurons", DE_ABcom_allDEGs_final$cell_type)] <- "CNS"
DE_ABcom_allDEGs_final$class[grep("Radial glia|Premature oligodendrocyte", DE_ABcom_allDEGs_final$cell_type)] <- "Radial_glia"
DE_ABcom_allDEGs_final$class[grep("Cardiac muscle lineages", DE_ABcom_allDEGs_final$cell_type)] <- "Heart"
DE_ABcom_allDEGs_final$class[grep("Hepatocytes & Primitive erythroid lineage", DE_ABcom_allDEGs_final$cell_type)] <- "Liver"

PCG_list <- list()
for (organ in c("CNS","Radial_glia","Heart", "Liver")){
  # organ <- "CNS"
  gene_list_PCG <- unique(DE_ABcom_allDEGs_final$gene[DE_ABcom_allDEGs_final$class==organ])
  PCG_list[[organ]] <- gene_list_PCG
  
  gene_list <- sort(gene_list_PCG, decreasing = TRUE)
  gene_id_list <- bitr(gene_list,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Mm.eg.db)
  gene_id_list <- dplyr::distinct(gene_id_list,SYMBOL,.keep_all=TRUE) # remove duplication
  gene_id_list_sort <- sort(gene_id_list$ENTREZID, decreasing = TRUE)
  
  go_enrich_BP <- enrichGO(gene = gene_id_list_sort,
                           OrgDb = org.Mm.eg.db,
                           keyType = 'ENTREZID',
                           readable = T,
                           ont = "BP",
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH")
  go_enrich_BP <- as.data.frame(go_enrich_BP)
  write.table(go_enrich_BP,paste0(organ, "_DE_ABcom_allDEG_enrichGO_GOBP.txt"),sep="\t",col.names=T,row.names=F,quote=F)
  
  go_enrich_KEGG <- enrichKEGG(gene = gene_id_list_sort,
                               organism = "hsa",
                               keyType = "kegg",
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH")
  go_enrich_KEGG <- as.data.frame(go_enrich_KEGG)
  write.table(go_enrich_KEGG,paste0(organ, "_DE_ABcom_allDEG_enrichKEGG_KEGG.txt"),sep="\t",col.names=T,row.names=F,quote=F)
}


### barplot for select GO terms
library(ggplot2)
library(forcats)

setwd(paste0(work_path, "/6_seurat_bigSpotShape_new"))

my_color <- c("#6766AD", brewer.pal(9, "Set1")[c(6,8,7)])
names(my_color) <- c("CNS","Radial_glia","Heart", "Liver")

GO_barplot_all <- c()
GO_barplot_list <- list()
for (organ in c("CNS","Radial_glia","Heart", "Liver")){
  # organ <- "Radial_glia"
  go_enrich_BP <- read.csv(paste0(organ, "_DE_ABcom_allDEG_enrichGO_GOBP.txt"), sep="\t", header=T)
  GO_barplot <- data.frame(go_enrich_BP[which(go_enrich_BP$highlight==1), ][, c("Description", "pvalue", "Count")], organ)[1:5,]
  GO_barplot$pvalue <- -log10(GO_barplot$pvalue)
  colnames(GO_barplot) <- c("GO_term", "pvalue", "Count", "organ")
  GO_barplot$GO_term <- factor(GO_barplot$GO_term, levels = rev(GO_barplot$GO_term))
  GO_barplot_all <- rbind(GO_barplot_all, GO_barplot)
}

GO_barplot_all$organ <- factor(GO_barplot_all$organ, levels = c("CNS","Radial_glia","Heart", "Liver"))
p <- ggplot(GO_barplot_all, aes(x=GO_term,y=pvalue, fill=organ)) + geom_bar(stat="identity", width=0.8) + coord_flip() +  xlab("GOBP") + ylab("-log10(pvalue)") +
        facet_wrap(~organ, scales = "free_y") +
        theme_bw() + 
        theme(legend.position="none") + 
        theme(panel.grid = element_blank()) + 
        scale_fill_manual(values=my_color)

pdf("/media/maoni/data/CZP/Figures/Spatial_cluster/Fig2_example_ABcom/DE_ABcom_allDEG_enrichGO_GOBP_barplot_all.pdf", width = 9, height = 4)
print(p)
dev.off()

########################################## revision (Rows and columns can be visualized on the UMAPs)
# cluster highlight
source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_hic/E1305-hic-24"
sample <- "E1305-hic-24"
higashi_out <- "8_higashi_out_higashi_500kb"
work_path_seurat <- "10_seurat_bigSpotShape_higashi_500kb"
embed_type <- "scAB_scale"
setwd(paste0(work_path,  "/", work_path_seurat, "/", embed_type))

load(paste0("Out/", sample, "_Seurat_object_anno.RData"))
Seurat_object <- SCTransform(Seurat_object, assay = "Spatial", verbose = T, variable.features.n = 2000)
Seurat_object <- RunPCA(Seurat_object, assay = "SCT", verbose = FALSE)
ElbowPlot(Seurat_object, ndims = 50)
DimPlot(Seurat_object, reduction = "umap", group.by = "seurat_clusters", label = FALSE)

Seurat_object$row <- Seurat_object@images$E1305.hic.24@coordinates$row
Seurat_object$col <- Seurat_object@images$E1305.hic.24@coordinates$col

for(i in 1:50){
  # i <- 1
  Seurat_object$row_number <- NA
  Seurat_object$row_number[which(Seurat_object$row==i)] <- "highlight"
  Seurat_object$row_number[which(Seurat_object$row!=i)] <- "other"
  Seurat_object[[paste0("row_", i)]] <- Seurat_object$row_number
}

for(i in 1:50){
  # i <- 1
  Seurat_object$col_number <- NA
  Seurat_object$col_number[which(Seurat_object$col==i)] <- "highlight"
  Seurat_object$col_number[which(Seurat_object$col!=i)] <- "other"
  Seurat_object[[paste0("col_", i)]] <- Seurat_object$col_number
}

### translate to adata (python)
dim(Seurat_object@meta.data)

library(reticulate)
use_python("/home/maoni/miniconda3/envs/r4.3/bin/python")
seurat2scanpy(Seurat_object, ann.X = "Spatial-data", ann.raw.X = NULL, h5ad_path = paste0("Out/", sample, "_Seurat_object_col_row_umap.h5ad"))

system(paste0("/home/maoni/miniconda3/envs/squidpy/bin/python /media/maoni/data/CZP/Figures/revision/UMAP_col_row/ARI_col.py ", sample))
system(paste0("/home/maoni/miniconda3/envs/squidpy/bin/python /media/maoni/data/CZP/Figures/revision/UMAP_col_row/ARI_row.py ", sample))
system(paste0("/home/maoni/miniconda3/envs/squidpy/bin/python /media/maoni/data/CZP/Figures/revision/UMAP_col_row/ARI_col_row.py ", sample))


DimPlot2_list <- list()

for(i in 1:50){
  
  #1) Dimplot
  p1 <- DimPlot(Seurat_object, reduction = "umap", group.by = paste0("row_", i), label = FALSE, order = TRUE, cells.highlight = which(Seurat_object$row==i), cols.highlight = "red", sizes.highlight = 1) + NoLegend() + ggtitle(paste0("row_", i)) 
  
  #2) get range
  getRange=function(x){
    #min(x) + 0.25 * diff(x)
    min(x) + 0.25 * (max(x)-min(x))
  }
  
  #3) set range
  DimPlot2_list[[i]] <- p1 + 
    scale_x_continuous(breaks = getRange(p1$data[,1]), guide = guide_axis(cap = 'upper')) +
    #scale_y_continuous(breaks = quantile(p1$data[,2], prob = 0.20), guide = guide_axis(cap = 'upper')) +
    scale_y_continuous(breaks = getRange(p1$data[,2]), guide = guide_axis(cap = 'upper')) +
    theme(aspect.ratio = 1,
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(arrow = arrow(type = "closed", length = unit(0.2, "cm"))),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(hjust = 0.025),
    )
}
pdf(paste0("/media/maoni/data/CZP/Figures/revision/UMAP_col_row/", sample, "_row_umap.pdf"), height = 9*3, width = 6*3)
ggarrange(plotlist = DimPlot2_list, ncol = 6, nrow = 9)
dev.off()


DimPlot2_list <- list()

for(i in 1:50){
  
  #1) Dimplot
  p1 <- DimPlot(Seurat_object, reduction = "umap", group.by = paste0("col_", i), label = FALSE, order = TRUE, cells.highlight = which(Seurat_object$col==i), cols.highlight = "red", sizes.highlight = 1) + NoLegend() + ggtitle(paste0("col_", i)) 
  
  #2) get range
  getRange=function(x){
    #min(x) + 0.25 * diff(x)
    min(x) + 0.25 * (max(x)-min(x))
  }
  
  #3) set range
  DimPlot2_list[[i]] <- p1 + 
    scale_x_continuous(breaks = getRange(p1$data[,1]), guide = guide_axis(cap = 'upper')) +
    #scale_y_continuous(breaks = quantile(p1$data[,2], prob = 0.20), guide = guide_axis(cap = 'upper')) +
    scale_y_continuous(breaks = getRange(p1$data[,2]), guide = guide_axis(cap = 'upper')) +
    theme(aspect.ratio = 1,
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(arrow = arrow(type = "closed", length = unit(0.2, "cm"))),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(hjust = 0.025),
    )
}
pdf(paste0("/media/maoni/data/CZP/Figures/revision/UMAP_col_row/", sample, "_col_umap.pdf"), height = 9*3, width = 6*3)
ggarrange(plotlist = DimPlot2_list, ncol = 6, nrow = 9)
dev.off()

pdf(paste0("/media/maoni/data/CZP/Figures/revision/UMAP_col_row/", sample, "_UMAP_col_row_ARI_plot.pdf"), width = 10, height = 5)
par(mfrow=c(1,2))
ARI_col <- read.csv(paste0("/media/maoni/data/CZP/Figures/revision/UMAP_col_row/", sample, "_ARI_col_vs_seurat_clusters.csv"), header = TRUE)
plot(ARI_col$ARI, type = "p", pch=16, col = "#E69F00", ylim = c(-0.1, 0.1), xlab = "col_number", ylab = "ARI_value", main = "ARI_col")
ARI_row <- read.csv(paste0("/media/maoni/data/CZP/Figures/revision/UMAP_col_row/", sample, "_ARI_row_vs_seurat_clusters.csv"), header = TRUE)
plot(ARI_row$ARI, type = "p", pch=16, col = "#56B4E9", ylim = c(-0.1, 0.1), xlab = "row_number", ylab = "ARI_value", main = "ARI_row")
dev.off()

