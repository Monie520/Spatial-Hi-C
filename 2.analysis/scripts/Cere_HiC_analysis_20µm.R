set.seed(123)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
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
# 0 file preparing
# source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_hic/Cere-hic-16"
sample <- "Cere-hic-16"
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

file.copy(from=paste0(work_path, "/", higashi_out, "/", sample, "_scAB_bin.csv"), 
          to=paste0(work_path,  "/", work_path_seurat, "/", embed_type, "/Matrix/"), overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)

#################################################################### E13-13 ##########################################################
# 1 Load data from st pipe
setwd(paste0(work_path, "/", work_path_seurat, "/", embed_type))
stat_cropped <- read.table("Matrix/stat_cropped.csv", sep =",", header = TRUE, dec =".", stringsAsFactors = F)
iA_iB <- paste(97-stat_cropped$iA, stat_cropped$iB, sep="x")
dim(stat_cropped)

# setwd("/media/maoni/data/CZP/spatial_hic/Cere-hic-13/8_higashi_out_higashi_500kb")
# aa <- read.table("Cere-hic-13_scAB_zscore.csv", sep =",", header = F, dec =".", stringsAsFactors = F)
# count_matrix <- aa

count_matrix <- read.table(paste0("Matrix/", sample, "_scAB_scale.csv"), sep =",", header = F, dec =".", stringsAsFactors = F)
# count_matrix <- matrix(runif(length(iA_iB)*4813), nrow = length(iA_iB), ncol = 4813)
rownames(count_matrix) <- iA_iB
count_matrix <- as.matrix(count_matrix)

bin_500kb <- read.table(paste0("Matrix/", sample, "_scAB_bin.csv"), sep =",", header = F, dec =".", stringsAsFactors = F)

All_ABcom <- bin_500kb$V1
All_ABcom_chr <- unlist(strsplit(All_ABcom, ":"))[seq(from = 1, to = length(unlist(strsplit(All_ABcom, ":"))), by = 2)]
All_ABcom_pos <- unlist(strsplit(All_ABcom, ":"))[seq(from = 2, to = length(unlist(strsplit(All_ABcom, ":"))), by = 2)]
All_ABcom_start <- unlist(strsplit(All_ABcom_pos, "-"))[seq(from = 2, to = length(unlist(strsplit(All_ABcom_pos, "-"))), by = 2)]
All_ABcom_end <- unlist(strsplit(All_ABcom_pos, "-"))[seq(from = 1, to = length(unlist(strsplit(All_ABcom_pos, "-"))), by = 2)]
All_ABcom_region <- data.frame(All_ABcom_chr, All_ABcom_start, All_ABcom_end, All_ABcom)
All_ABcom_region$All_ABcom <- paste(paste(All_ABcom_region$All_ABcom_chr,All_ABcom_region$All_ABcom_start, sep=":"), All_ABcom_region$All_ABcom_end, sep="-")
bin_500kb$V1 <- All_ABcom_region$All_ABcom

#bin_500kb <- read.table(paste0(work_path, "/Res_500kb_col_names.txt"), sep ="\t", header = FALSE, stringsAsFactors = F)
colnames(count_matrix) <- bin_500kb$V1
dim(count_matrix) # 9216 4812

emb.m <- as.matrix(read.table(paste0("Matrix/", sample, "_", embed_type, ".csv"), sep=",", header = F))
rownames(emb.m) <- rownames(count_matrix)
dim(emb.m)  # 1999  256

assay = "Spatial"
slice = sample
Seurat_object = CreateSeuratObject(counts = t(count_matrix), project = sample, assay = assay) # min.cells min.features
image.dir = "/media/maoni/data/CZP/spatial_hic/seurat/core/96*96_blank_img" # "./Img"
image.nam = paste0(sample,"_fix.png") # "grey_pixel_1080p.png"
coord.nam = "combine_barcode.round2round1_index1_index2.Seurat_2.txt"
image <- readPNG(source = file.path(image.dir, image.nam))[,,1:3]
scale.factors <- c("tissue_hires_scalef"=1, "fiducial_diameter_fullres"=1, "tissue_lowres_scalef"=1)
tissue.positions <- read.table(file = file.path(image.dir,coord.nam), col.names = c("barcodes", "tissue", "row", "col", "imagerow", "imagecol"), header = FALSE, as.is = TRUE, row.names = 1)
spot.radius <- 0.015 # estiamte:(0.13)*50/410/2
image <- new(Class = "VisiumV1", image = image, scale.factors = scalefactors(spot = scale.factors[1], fiducial = scale.factors[2], hires = scale.factors[1], lowres = scale.factors[3]), coordinates = tissue.positions, spot.radius = spot.radius)
image <- image[Cells(Seurat_object)]
DefaultAssay(object = image) <- assay
Seurat_object[[slice]] <- image

# 3 Data quality
pdf("Plot/0_nCounts.Seurat.pdf", width = 14, height = 7)
plot1 <- VlnPlot(Seurat_object, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(Seurat_object, features = "nCount_Spatial", max.cutoff = 2000, pt.size.factor = 0.8, stroke = NA) + theme(legend.position = "right")
plot2$layers[[1]]$aes_params=c(plot2$layers[[1]]$aes_params, shape=22)
plot_grid(plot1, plot2)
dev.off()

pdf("Plot/0_nCounts.Seurat_filter.pdf", width = 14, height = 7)
plot1 <- SpatialFeaturePlot(Seurat_object, features = "nCount_Spatial",crop = F,max.cutoff = 2000, pt.size.factor = 0.8, stroke = NA) + theme(legend.position = "right")
plot1$layers[[1]]$aes_params=c(plot1$layers[[1]]$aes_params, shape=22)
plot2 <- SpatialFeaturePlot(Seurat_object, features = "nFeature_Spatial",crop = F,max.cutoff = 1500, pt.size.factor = 0.8, stroke = NA) + theme(legend.position = "right")
plot2$layers[[1]]$aes_params=c(plot2$layers[[1]]$aes_params, shape=22)
plot_grid(plot1, plot2)
dev.off()

# 4 Run PCA
Seurat_object <- SCTransform(Seurat_object, assay = "Spatial", verbose = T, variable.features.n = 2000)
Seurat_object <- RunPCA(Seurat_object, assay = "SCT", npcs = 50, verbose = FALSE)
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
  spatial_p[[paste0("res",res)]] <- SpatialDimPlot(Seurat_object, label = TRUE, label.size = 3, pt.size.factor = 0.8, stroke = NA) + ggtitle(paste0("res",res)) # the spatial plot
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
for(pcNO in seq(10, 50, by=5)){
  Seurat_object <- FindNeighbors(Seurat_object, reduction = "pca", dims = 1:pcNO)
  Seurat_object <- FindClusters(Seurat_object, verbose = FALSE, resolution = res)
  Seurat_object <- RunUMAP(Seurat_object, reduction = "pca", dims = 1:pcNO)
  dim_p[[paste0("PC",pcNO)]] <- DimPlot(Seurat_object, reduction = "umap", label = TRUE) + ggtitle(paste0("PC",pcNO)) # the UMAP plot
  spatial_p[[paste0("PC",pcNO)]] <- SpatialDimPlot(Seurat_object, label = TRUE, label.size = 3, pt.size.factor = 0.8, stroke = NA) + ggtitle(paste0("PC",pcNO)) # the spatial plot
  spatial_p[[paste0("PC",pcNO)]]$layers[[1]]$aes_params=c(spatial_p[[paste0("PC",pcNO)]]$layers[[1]]$aes_params, shape=22)
}

pdf(file = paste("Plot/2_UMAP_clusters_res",res,".pdf",sep =""), width = 28, height = ceiling(length(dim_p)/4)*7)
ggarrange(plotlist = dim_p, ncol = 4, nrow = ceiling(length(dim_p)/4))
dev.off()

pdf(file = paste("Plot/2_Spatial_clusters_res",res,".pdf",sep =""), width = 28, height = ceiling(length(spatial_p)/4)*7)
ggarrange(plotlist = spatial_p, ncol = 4, nrow = ceiling(length(spatial_p)/4))
dev.off()


# Select PC and res
# for(res in c(1.2,1.5,2,1)){
  PC=36
  res=1
  Seurat_object <- FindNeighbors(Seurat_object, reduction = "pca", dims = 1:PC)
  Seurat_object <- FindClusters(Seurat_object, verbose = FALSE, resolution = res)
  Seurat_object <- RunUMAP(Seurat_object, reduction = "pca", dims = 1:PC)
  umap_Spa <- list()
  umap_Spa[["p1"]] <- DimPlot(Seurat_object, reduction = "umap", label = TRUE)
  umap_Spa[["p2"]] <- SpatialDimPlot(Seurat_object, label = TRUE, label.size = 3, pt.size.factor = 0.8, stroke = NA)
  umap_Spa[["p2"]]$layers[[1]]$aes_params=c(umap_Spa[["p2"]]$layers[[1]]$aes_params, shape=22)
  umap_Spa[["p3"]] <- DimPlot(Seurat_object, reduction = "umap", label = TRUE)
  umap_Spa[["p4"]] <- SpatialDimPlot(Seurat_object, label = TRUE, label.size = 3, pt.size.factor = 0.8, stroke = NA)
  umap_Spa[["p4"]]$layers[[1]]$aes_params=c(umap_Spa[["p4"]]$layers[[1]]$aes_params, shape=22)
  
  pdf(file = paste0("Plot/2_UMAP_clusters_PC",PC,"_res",res,".pdf"), width = 28, height = ceiling(length(umap_Spa)/4)*7)
  p <- ggarrange(plotlist = umap_Spa, ncol = 4, nrow = ceiling(length(umap_Spa)/4))
  print(p)
  dev.off()
  
  spatial_p_cluster <- list()
  for(i in sort(unique(Idents(Seurat_object)))){
    spatial_p_cluster[[paste0("idents_",i)]] <- SpatialDimPlot(Seurat_object, cells.highlight = CellsByIdentities(object = Seurat_object, idents = i), cols.highlight = c("#DE2D26", "grey90"), facet.highlight = TRUE, ncol = 1, pt.size.factor = 0.8, stroke = NA, label.size = 12) + ggtitle(paste0("idents",i)) # the spatial plot
    spatial_p_cluster[[paste0("idents_",i)]]$layers[[1]]$aes_params=c(spatial_p_cluster[[paste0("idents_",i)]]$layers[[1]]$aes_params, shape=22)
  }
  
  pdf(file = paste0("Plot/2_Spatial_clusters_PC",PC,"_res",res,".pdf"), width = 28, height = ceiling(length(spatial_p_cluster)/4)*7)
  p <- ggarrange(plotlist = spatial_p_cluster, ncol = 4, nrow = ceiling(length(spatial_p_cluster)/4))
  print(p)
  dev.off()
  
  pdf_combine(input = c(paste0("Plot/2_UMAP_clusters_PC",PC,"_res",res,".pdf"), paste0("Plot/2_Spatial_clusters_PC",PC,"_res",res,".pdf")), output = paste0("Plot/2_UMAP_Spatial_clusters_PC",PC,"_res",res,".pdf"))
# }


# 6.DE
de_markers_13 = FindAllMarkers(Seurat_object, test.use = "wilcox", only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
de_markers_13 %>% group_by(cluster) %>% slice_max(n=5, order_by = avg_log2FC) -> top5

DefaultAssay(Seurat_object) <- "Spatial"

pdf("Plot/3_findmarker_2D_impose.data_pc36_res1.pdf",height = length(table(top5$cluster))*7, width = 5*7)
plot <- SpatialFeaturePlot(object = Seurat_object, features = top5$gene, slot = "counts", alpha = c(0.1, 1), ncol = 5, pt.size.factor = 0.8)
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
DoHeatmap(Seurat_object, features = de_markers_13$gene, slot = "data", group.by = "seurat_clusters", size = 10) + scale_fill_gradientn(colors = c("white", "grey", "firebrick3")) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank())
dev.off()


write.csv(de_markers_13,row.names = F,file = "Out/1_de_ABcom.csv")
write.csv(Seurat_object@meta.data,row.names = T,file = "Out/1_cell_metadata.csv")
save(de_markers_13, top5, Seurat_object, file=paste0("Out/", sample, "_Seurat_object.RData"))


# 7. annotation
setwd(paste0(work_path,  "/", work_path_seurat, "/", embed_type))
PC=36
res=1
sample="Cere-hic-13"
load(paste0("Out/", sample, "_Seurat_object.RData"))

anno_marker <- rep(NA,nrow(de_markers_13))
anno_marker[de_markers_13$cluster %in% c(0)] <- "ML_1"
anno_marker[de_markers_13$cluster %in% c(1)] <- "ML_2"
anno_marker[de_markers_13$cluster %in% c(5)] <- "ML_3"
anno_marker[de_markers_13$cluster %in% c(9)] <- "ML_4"
anno_marker[de_markers_13$cluster %in% c(2,6)] <- "GL_1"
anno_marker[de_markers_13$cluster %in% c(11)] <- "GL_2"
anno_marker[de_markers_13$cluster %in% c(10)] <- "GL_3"
anno_marker[de_markers_13$cluster %in% c(4)] <- "GL_4"
anno_marker[de_markers_13$cluster %in% c(8)] <- "GL_5"
anno_marker[de_markers_13$cluster %in% c(7)] <- "PL_1"
anno_marker[de_markers_13$cluster %in% c(3)] <- "PL_2"
de_markers_13$cell_type <- factor(anno_marker,levels = c("GL_1", "GL_2", "GL_3", "GL_4", "GL_5", "PL_1", "PL_2", "ML_1", "ML_2", "ML_3", "ML_4"))

anno <- rep(NA,ncol(Seurat_object))
anno[Seurat_object$seurat_clusters %in% c(0)] <- "ML_1"
anno[Seurat_object$seurat_clusters %in% c(1)] <- "ML_2"
anno[Seurat_object$seurat_clusters %in% c(5)] <- "ML_3"
anno[Seurat_object$seurat_clusters %in% c(9)] <- "ML_4"
anno[Seurat_object$seurat_clusters %in% c(2,6)] <- "GL_1"
anno[Seurat_object$seurat_clusters %in% c(11)] <- "GL_2"
anno[Seurat_object$seurat_clusters %in% c(10)] <- "GL_3"
anno[Seurat_object$seurat_clusters %in% c(4)] <- "GL_4"
anno[Seurat_object$seurat_clusters %in% c(8)] <- "GL_5"
anno[Seurat_object$seurat_clusters %in% c(7)] <- "PL_1"
anno[Seurat_object$seurat_clusters %in% c(3)] <- "PL_2"
Seurat_object$cell_type <- factor(anno,levels = c("GL_1", "GL_2", "GL_3", "GL_4", "GL_5", "PL_1", "PL_2", "ML_1", "ML_2", "ML_3", "ML_4"))


source("/media/maoni/data/Zhanglin_Junjie/ColorPalettes.R") # 11 cell types
my_color_palette <- c("#E4A5F6", "#CE6CF5", "#BC18B5", "#C61385", "#791367", "#F3153C", "#FDBB15", "#7395BF", "#1380DD","#BE8476", "#584B5F")
scales::show_col(my_color_palette)
names(my_color_palette) <- c("GL_1", "GL_2", "GL_3", "GL_4", "GL_5", "PL_1", "PL_2", "ML_1", "ML_2", "ML_3", "ML_4")

umap_Spa <- list()
umap_Spa[["p1"]] <- DimPlot(Seurat_object, reduction = "umap", group.by = "cell_type", label = FALSE, pt.size = 0.5, cols = my_color_palette) + theme(aspect.ratio = 1)
umap_Spa[["p2"]] <- SpatialDimPlot(Seurat_object, group.by = "cell_type", label = FALSE, label.size = 3, pt.size.factor = 0.8, stroke = NA, cols = my_color_palette) + theme(aspect.ratio = 1)
umap_Spa[["p2"]]$layers[[1]]$aes_params=c(umap_Spa[["p2"]]$layers[[1]]$aes_params, shape=22)
umap_Spa[["p3"]] <- DimPlot(Seurat_object, reduction = "umap", group.by = "cell_type", label = FALSE, pt.size = 0.5, cols = my_color_palette)+ theme(aspect.ratio = 1)
umap_Spa[["p4"]] <- SpatialDimPlot(Seurat_object, group.by = "cell_type", label = FALSE, label.size = 3, pt.size.factor = 0.8, stroke = NA, cols = my_color_palette)+ theme(aspect.ratio = 1)
umap_Spa[["p4"]]$layers[[1]]$aes_params=c(umap_Spa[["p4"]]$layers[[1]]$aes_params, shape=22)

pdf(file = paste0("Plot/2_UMAP_clusters_PC",PC,"_res",res,"_final.pdf"), width = 28, height = ceiling(length(umap_Spa)/4)*7)
p <- ggarrange(plotlist = umap_Spa, ncol = 4, nrow = ceiling(length(umap_Spa)/4))
print(p)
dev.off()

spatial_p_cluster <- list()
for(i in names(table(Seurat_object$cell_type))){
  Seurat_object_subset <- Seurat_object[, Seurat_object$cell_type==i]
  spatial_p_cluster[[i]] <- SpatialDimPlot(Seurat_object_subset, group.by = "cell_type", crop = FALSE, label = FALSE, label.size = 3, pt.size.factor = 0.8, stroke = NA,  label.box = FALSE, cols = my_color_palette) + theme(aspect.ratio = 1) + ggtitle(i) + theme(legend.position="none") # the spatial plot
  spatial_p_cluster[[i]]$layers[[1]]$aes_params=c(spatial_p_cluster[[i]]$layers[[1]]$aes_params, shape=22)
}

pdf(file = paste0("Plot/2_Spatial_clusters_PC",PC,"_res",res,"_final.pdf"), width = 28, height = ceiling(length(spatial_p_cluster)/4)*7)
p <- ggarrange(plotlist = spatial_p_cluster, ncol = 4, nrow = ceiling(length(spatial_p_cluster)/4))
print(p)
dev.off()

pdf_combine(input = c(paste0("Plot/2_UMAP_clusters_PC",PC,"_res",res,"_final.pdf"), paste0("Plot/2_Spatial_clusters_PC",PC,"_res",res,"_final.pdf")), output = paste0("Plot/2_UMAP_Spatial_clusters_PC",PC,"_res",res,"_final.pdf"))

top10_anno <- de_markers_13 %>% group_by(cell_type) %>% top_n(n = 10, wt = avg_log2FC)  %>% arrange(cell_type, desc(avg_log2FC))
de_markers_13 <- de_markers_13 %>% group_by(cell_type) %>% arrange(cell_type, desc(avg_log2FC))

pdf(file = paste("Plot/3_top5_ABcom_heat_anno.pdf"), width = 20, height = 15)
DoHeatmap(Seurat_object, features = top10_anno$gene, slot = "data", group.by = "cell_type", size = 10) + scale_fill_gradientn(colors = c("white", "grey", "firebrick3")) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank())
dev.off()

pdf(file = paste("Plot/3_allDE_ABcom_heat_anno.pdf"), width = 20, height = 15)
DoHeatmap(Seurat_object, features = de_markers_13$gene, slot = "data", group.by = "cell_type", size = 10) + scale_fill_gradientn(colors = c("white", "grey", "firebrick3")) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank())
dev.off()

write.csv(de_markers_13,row.names = F,file = "Out/1_de_ABcom.csv")
write.csv(Seurat_object@meta.data,row.names = T,file = "Out/1_cell_metadata.csv")
save(de_markers_13, top10_anno, Seurat_object, file=paste0("Out/", sample, "_Seurat_object_anno.RData"))


# 8. hightlight
source("/media/maoni/data/R_functions_ST_E13.R")
Seurat_object_magic <- magic.Seurat(Seurat_object)
DefaultAssay(Seurat_object_magic) <- "MAGIC_Spatial"

png("Plot/3_findmarker_2D_impose.data_pc36_res1_magic.png",height = length(table(top10_anno$cell_type))*480, width = 5*480)
plot <- SpatialFeaturePlot(object = Seurat_object_magic, features = unique(top10_anno$gene), alpha = c(1, 1), ncol = 5, pt.size.factor = 0.8) # + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")))
plot$layers[[1]]$aes_params=c(plot$layers[[1]]$aes_params, shape=22)
plot
dev.off()


cell_type <- names(table(top10_anno$cell_type))
for(j in 1:length(cell_type)){
  spatial_p_feature <- list()
  for(i in 1:length(top10_anno$gene[top10_anno$cell_type==cell_type[j]])){
    DEAB_spe_select <- top10_anno$gene[top10_anno$cell_type==cell_type[j]][i]
    spatial_p_feature[[DEAB_spe_select]] <- SpatialFeaturePlot(object = Seurat_object_magic, features = DEAB_spe_select, alpha = c(1, 1), ncol = 1, pt.size.factor = 0.8, stroke = NA) + scale_fill_gradientn(values=c(0,0.3,0.5,0.6,0.7,0.8,1),colors = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), na.value = "white") + theme(aspect.ratio = 1)
    spatial_p_feature[[DEAB_spe_select]]$layers[[1]]$aes_params=c(spatial_p_feature[[DEAB_spe_select]]$layers[[1]]$aes_params, shape=22)
  }
  png(paste0("Plot/3_findmarker_2D_impose.data_pc36_res1_magic_", cell_type[j], ".png"),height = 2*480, width = 5*480)
  p <- ggarrange(plotlist = spatial_p_feature, ncol = 5, nrow = 2)
  print(p)
  dev.off()
}


de_spe_select <- top10_anno$gene[c(1,11,23,35,41,52,61,77,81,7,13,24,36,42,59,62,79,82)]
spatial_p_feature <- list()
for(i in 1:length(de_spe_select)){
  spatial_p_feature[[paste0("feature",i)]] <- SpatialFeaturePlot(object = Seurat_object_magic, features = de_spe_select[i], alpha = c(1, 1), ncol = 1, pt.size.factor = 0.8, stroke = NA) + scale_fill_gradientn(values=c(0,0.3,0.5,0.6,0.7,0.8,1),colors = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), na.value = "white") # + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) #  + theme(legend.position="right")
  spatial_p_feature[[paste0("feature",i)]]$layers[[1]]$aes_params=c(spatial_p_feature[[paste0("feature",i)]]$layers[[1]]$aes_params, shape=22)
}
pdf(file = paste("Plot/4_findmarker_2D_impose.data_pc36_res1_de_spe_select.pdf"), width = 9*5, height = ceiling(length(spatial_p_feature)/9)*7)
ggarrange(plotlist = spatial_p_feature, ncol = 9, nrow = ceiling(length(spatial_p_feature)/9))
dev.off()

save(Seurat_object_magic, file=paste0("Out/", sample, "_Seurat_object_magic.RData"))
load(paste0("Out/", sample, "_Seurat_object_magic.RData"))

########################################################################## Downstream analysis #################################################################

source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_hic/Cere-hic-13"
sample <- "Cere-hic-13"
higashi_out <- "8_higashi_out_higashi_500kb"
work_path_seurat <- "10_seurat_bigSpotShape_higashi_500kb"
embed_type <- "scAB_scale"

setwd(paste0(work_path,  "/", work_path_seurat, "/", embed_type))

load(paste0("Out/", sample, "_Seurat_object_anno.RData"))

###########################  prepare cropped_stat for each cluster
# hicup: prepare stat_cropped_barcode files for clusters
setwd(work_path)
stat_cropped <- read.table("6_crop/stat_cropped.csv", sep =",", header = TRUE, dec =".", stringsAsFactors = F)
iA_iB <- paste(97-stat_cropped$iA, stat_cropped$iB, sep="x")


for(i in names(table(Seurat_object$cell_type))){
  stat_cropped_barcode <- stat_cropped[match(colnames(Seurat_object)[Seurat_object$cell_type==i], iA_iB), ]
  stat_cropped_barcode <- data.frame(stat_cropped_barcode$pixel, paste0("Cereall.R1.", i, ".hicup"))
  write.table(stat_cropped_barcode, paste0("6_crop/stat_cropped_barcode_", i), sep="\t", row.names=F, col.names=F, quote=F)
}

# hicexplorer: prepare stat_cropped_barcode files for clusters
setwd(work_path)
stat_cropped <- read.table("6_crop/stat_cropped.csv", sep =",", header = TRUE, dec =".", stringsAsFactors = F)
iA_iB <- paste(97-stat_cropped$iA, stat_cropped$iB, sep="x")

for(i in names(table(Seurat_object$cell_type))){
  stat_cropped_barcode <- stat_cropped[match(colnames(Seurat_object)[Seurat_object$cell_type==i], iA_iB), ]
  
  stat_cropped_barcode_R1 <- data.frame(stat_cropped_barcode$pixel, paste0(sample, "_R1_", i))
  write.table(stat_cropped_barcode_R1, paste0("6_crop/stat_cropped_barcode_R1_", i), sep="\t", row.names=F, col.names=F, quote=F)
  
  stat_cropped_barcode_R2 <- data.frame(stat_cropped_barcode$pixel, paste0(sample, "_R2_", i))
  write.table(stat_cropped_barcode_R2, paste0("6_crop/stat_cropped_barcode_R2_", i), sep="\t", row.names=F, col.names=F, quote=F)
}



########################################################################## DE_ABcom heatmap and function analysis ##########################################################
########################### DE_ABcom heatmap 
# mean of each cell type
dim(Seurat_object@assays$Spatial@counts)
dim(de_markers_13)
count_matrix <- Seurat_object@assays$Spatial@counts


count_matrix_de_markers <- count_matrix[match(de_markers_13$gene, rownames(count_matrix)),]
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


count_matrix_de_markers_norm_mean$cell_type <- c("GL", "GL", "GL", "GL", "GL", "PL", "PL", "ML_1", "ML_2", "ML_3", "ML_4")
count_matrix_de_markers_norm_mean_mean <- aggregate(. ~ cell_type, data = count_matrix_de_markers_norm_mean, FUN = mean)
rownames(count_matrix_de_markers_norm_mean_mean) <- count_matrix_de_markers_norm_mean_mean$cell_type
count_matrix_de_markers_norm_mean_mean <- count_matrix_de_markers_norm_mean_mean[,-1]
count_matrix_de_markers_norm_mean_mean_trans <- t(count_matrix_de_markers_norm_mean_mean)

temp <- as.data.frame(count_matrix_de_markers_norm_mean_mean_trans)
temp <- temp[, c("GL", "PL", "ML_1", "ML_2", "ML_3", "ML_4")]

library(pheatmap)

annotation_col = data.frame(
  cluster = factor(c("GL_1", "GL_2", "GL_3", "GL_4", "GL_5", "PL_1", "PL_2", "ML_1", "ML_2", "ML_3", "ML_4"))
)
rownames(annotation_col) = c("GL_1", "GL_2", "GL_3", "GL_4", "GL_5", "PL_1", "PL_2", "ML_1", "ML_2", "ML_3", "ML_4")

ann_colors = list(
  cluster = c(GL_1 = "#E4A5F6", GL_2 = "#CE6CF5", GL_3 = "#BC18B5", GL_4 = "#C61385", GL_5 = "#791367", PL_1 = "#F3153C",
              PL_2 = "#FDBB15", ML_1 = "#7395BF", ML_2 = "#1380DD", ML_3 = "#BE8476", ML_4 = "#584B5F")
)

################################ Acom_spe
GL_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==1), which(apply(temp, 1, max) > 0.5))]
PL_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==2), which(apply(temp, 1, max) > 0.5))]
ML_1_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==3), which(apply(temp, 1, max) > 0.5))]
ML_2_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==4), which(apply(temp, 1, max) > 0.5))]
ML_3_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==5), which(apply(temp, 1, max) > 0.75))]
ML_4_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==6), which(apply(temp, 1, max) > 0.75))]

temp_Acom_spe <- c(GL_spe, PL_spe, ML_1_spe, ML_2_spe, ML_3_spe, ML_4_spe)
temp_Acom_spe_matrix <- count_matrix_de_markers_norm_mean_trans[match(temp_Acom_spe, rownames(count_matrix_de_markers_norm_mean_trans)), ]

temp_Acom_spe_2 <- c(GL_spe, PL_spe, ML_1_spe, ML_2_spe)
temp_Acom_spe_matrix_2 <- count_matrix_de_markers_norm_mean_trans[match(temp_Acom_spe_2, rownames(count_matrix_de_markers_norm_mean_trans)), ][, 1:9]

length(GL_spe) + length(PL_spe) + length(ML_1_spe) + length(ML_2_spe) + length(ML_3_spe) + length(ML_4_spe) # 213
length(GL_spe) # 68
length(PL_spe) # 33
length(ML_1_spe) # 13
length(ML_2_spe) # 1
length(ML_3_spe) # 49
length(ML_4_spe) # 49

pdf(file = paste("Plot/3_allDE_Acom_heat_anno_average.pdf"), width = 4, height = 4)
# pheatmap(temp_Acom_spe_matrix, border = F, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, color = colorRampPalette(c("#C44C9C", "#B96397", "#AA7591", "#968686", "#86917D", "#739D6F", "#5CA65F", "#1CB24B"))(100))
pheatmap(temp_Acom_spe_matrix, border = F, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, color = colorRampPalette(c("#36359D", "#0F74DC", "#05A5C6", "#4EBD90", "#A9BE69", "#F3B947", "#F4E61A"))(100))
dev.off()

pdf(file = paste("Plot/3_allDE_Acom_heat_anno_average_2.pdf"), width = 4, height = 4)
# pheatmap(temp_Acom_spe_matrix, border = F, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, color = colorRampPalette(c("#C44C9C", "#B96397", "#AA7591", "#968686", "#86917D", "#739D6F", "#5CA65F", "#1CB24B"))(100))
pheatmap(temp_Acom_spe_matrix_2, border = F, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, color = colorRampPalette(c("#36359D", "#0F74DC", "#05A5C6", "#4EBD90", "#A9BE69", "#F3B947", "#F4E61A"))(100))
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
                                features = list(GL_spe=GL_spe, PL_spe=PL_spe, ML_spe=c(ML_1_spe, ML_2_spe)),
                                name = c("GL", "PL", "ML"),
                                assay = "Spatial")

spatial_p_feature <- list()
for(feature in c("GL1", "PL2", "ML3")){
  spatial_p_feature[[paste0("umap_",feature)]] <- FeaturePlot(Seurat_object, features = feature, label = FALSE, repel = TRUE) + scale_color_gradientn(colors = c("#36359D", "#0F74DC", "#05A5C6", "#4EBD90", "#A9BE69", "#F3B947", "#F4E61A")) + theme(aspect.ratio = 1)  # c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407")
  spatial_p_feature[[paste0("spatial_",feature)]] <- SpatialFeaturePlot(object = Seurat_object, features = feature, alpha = c(1, 1), ncol = 1, pt.size.factor = 0.8, stroke = NA) + scale_fill_gradientn(colors = c("#36359D", "#0F74DC", "#05A5C6", "#4EBD90", "#A9BE69", "#F3B947", "#F4E61A"), na.value = "white") + theme(aspect.ratio = 1) #  + theme(legend.position="right")
  spatial_p_feature[[paste0("spatial_",feature)]]$layers[[1]]$aes_params=c(spatial_p_feature[[paste0("spatial_",feature)]]$layers[[1]]$aes_params, shape=22)
} 
pdf(file = paste("Plot/3_allDE_Acom_FeaturePlot_anno_average_2.pdf"), width = 4*7, height = 3*7)
p <- ggarrange(plotlist = spatial_p_feature, ncol = 4, nrow = 3)
print(p)
dev.off()


########################### Bcom_spe
GL_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.min)==1), which(apply(temp, 1, min) < -0.35))]
PL_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.min)==2), which(apply(temp, 1, min) < -0.35))]
ML_1_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.min)==3), which(apply(temp, 1, min) < -0.35))]
ML_2_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.min)==4), which(apply(temp, 1, min) < -0.35))]
ML_3_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.min)==5), which(apply(temp, 1, min) < -0.75))]
ML_4_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.min)==6), which(apply(temp, 1, min) < -0.75))]

temp_Acom_spe <- c(GL_spe, PL_spe, ML_1_spe, ML_2_spe, ML_3_spe, ML_4_spe)
temp_Acom_spe_matrix <- count_matrix_de_markers_norm_mean_trans[match(temp_Acom_spe, rownames(count_matrix_de_markers_norm_mean_trans)), ]

length(GL_spe) + length(PL_spe) + length(ML_1_spe) + length(ML_2_spe) + length(ML_3_spe) + length(ML_4_spe) # 184
length(GL_spe) # 69
length(PL_spe) # 15
length(ML_1_spe) # 16
length(ML_2_spe) # 3
length(ML_3_spe) # 23
length(ML_4_spe) # 58

pdf(file = paste("Plot/3_allDE_Bcom_heat_anno_average.pdf"), width = 4, height = 4)
#pheatmap(temp_Acom_spe_matrix, border = F, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, color = colorRampPalette(c("#C44C9C", "#B96397", "#AA7591", "#968686", "#86917D", "#739D6F", "#5CA65F", "#1CB24B"))(100))
pheatmap(temp_Acom_spe_matrix, border = F, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, color = colorRampPalette(c("#36359D", "#0F74DC", "#05A5C6", "#4EBD90", "#A9BE69", "#F3B947", "#F4E61A"))(100))
dev.off()

### AddModuleScore (Bcom)
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
                       features = list(GL_spe=GL_spe, PL_spe=PL_spe, ML_1_spe=ML_1_spe, ML_2_spe=ML_2_spe, ML_3_spe=ML_3_spe, ML_4_spe=ML_4_spe),
                       name = c("GL", "PL", "ML_1", "ML_2", "ML_3", "ML_4"),
                       assay = "Spatial")

require(stats)
x <- matrix(1:10, ncol = 2)
(centered.x <- scale(x, scale = FALSE))
cov(centered.scaled.x <- scale(x)) # all 1

spatial_p_feature <- list()
for(feature in c("GL1", "PL2", "ML_13", "ML_24", "ML_35", "ML_46")){
  spatial_p_feature[[paste0("umap_",feature)]] <- FeaturePlot(Seurat_object, features = feature, label = FALSE, repel = TRUE) + scale_color_gradientn(colors = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407")) + theme(aspect.ratio = 1)
  spatial_p_feature[[paste0("spatial_",feature)]] <- SpatialFeaturePlot(object = Seurat_object, features = feature, alpha = c(1, 1), ncol = 1, pt.size.factor = 0.8, stroke = NA) + scale_fill_gradientn(colors = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), na.value = "white") + theme(aspect.ratio = 1) #  + theme(legend.position="right")
  spatial_p_feature[[paste0("spatial_",feature)]]$layers[[1]]$aes_params=c(spatial_p_feature[[paste0("spatial_",feature)]]$layers[[1]]$aes_params, shape=22)
} 
pdf(file = paste("Plot/3_allDE_Bcom_FeaturePlot_anno_average.pdf"), width = 4*7, height = 3*7)
p <- ggarrange(plotlist = spatial_p_feature, ncol = 4, nrow = 3)
print(p)
dev.off()



###################################################  Identify DEGs located in significant ABcom regions ################################################### 
###########################  load Hi-C data
source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_hic/Cere-hic-13"
sample <- "Cere-hic-13"
higashi_out <- "8_higashi_out_higashi_500kb"
work_path_seurat <- "10_seurat_bigSpotShape_higashi_500kb"
embed_type <- "scAB_scale"

setwd(paste0(work_path,  "/", work_path_seurat, "/", embed_type))
load(paste0("Out/", sample, "_Seurat_object_anno.RData"))

# top10_hic <- de_markers_13 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)  %>% arrange(cluster, desc(avg_log2FC))
DE_ABcom <- unique(de_markers_13$gene)

DE_ABcom_chr <- unlist(strsplit(DE_ABcom, ":"))[seq(from = 1, to = length(unlist(strsplit(DE_ABcom, ":"))), by = 2)]
DE_ABcom_pos <- unlist(strsplit(DE_ABcom, ":"))[seq(from = 2, to = length(unlist(strsplit(DE_ABcom, ":"))), by = 2)]
DE_ABcom_start <- unlist(strsplit(DE_ABcom_pos, "-"))[seq(from = 1, to = length(unlist(strsplit(DE_ABcom_pos, "-"))), by = 2)]
DE_ABcom_end <- unlist(strsplit(DE_ABcom_pos, "-"))[seq(from = 2, to = length(unlist(strsplit(DE_ABcom_pos, "-"))), by = 2)]
DE_ABcom_region <- data.frame(DE_ABcom_chr, DE_ABcom_start, DE_ABcom_end, DE_ABcom)
write.table(DE_ABcom_region, "Plot/DE_ABcom_region.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

### command run
Reference=/media/maoni/data/Reference/mouse/GRCm38_mm10
AB_path=/media/maoni/data/CZP/spatial_hic/Cere-hic-13/10_seurat_bigSpotShape_higashi_500kb/scAB_scale/Plot
intersectBed -a $AB_path/DE_ABcom_region.bed -b $Reference/ensembl/All_gene_Anno_PCG.bed -wa -wb > $AB_path/DE_ABcom_region_PCGs.bed

DE_ABcom_region_PCGs <- read.table("Plot/DE_ABcom_region_PCGs.bed", sep = "\t", header = FALSE)
DE_ABcom_region_PCGs_unique <- unique(DE_ABcom_region_PCGs$V8)


########################### load RNA data
source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_transcriptome/Cere-RNA-6"
sample <- "Cere-RNA-6"
setwd(paste0(work_path, "/6_seurat_bigSpotShape_new"))
PC=20
res=2.5
load(paste0("Out/", sample, "_cropped_PC",PC,"_res",res,"_anno.RData"))
Seurat_object_cropped_magic <- readRDS(paste0("Out/", sample, "_cropped_magic.rds"))

top50_anno <- de_markers %>% group_by(cell_type) %>% top_n(n = 50, wt = avg_log2FC) %>% arrange(cell_type, desc(avg_log2FC))
all_anno <- de_markers %>% group_by(cell_type) %>% arrange(cell_type, desc(avg_log2FC))

### DEG DE_ABcom overlap
DE_ABcom_allDEGs <- data.frame(all_anno[which(all_anno$gene %in% DE_ABcom_region_PCGs_unique),])
DE_ABcom_allDEGs_final <- c()
for(i in 1:length(DE_ABcom_allDEGs$gene)){
  ABcom <- paste(DE_ABcom_region_PCGs$V4[which(DE_ABcom_region_PCGs$V8 %in% DE_ABcom_allDEGs$gene[i])], collapse=";")
  tmp <- cbind(DE_ABcom_allDEGs[i,], ABcom)
  DE_ABcom_allDEGs_final <- rbind(DE_ABcom_allDEGs_final, tmp)
}
write.csv(DE_ABcom_allDEGs_final, "Plot/DE_ABcom_allDEGs_final.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

pdf(paste0("Plot/3_findmarker_dotplot_data_PC",PC,"_res",res,"_DE_ABcom_allDEGs_final.pdf"), height = 30, width = 16)
p <- DotPlot(Seurat_object_cropped, features = unique(DE_ABcom_allDEGs_final$gene), cols = c("grey", "red"), group.by = 'cell_type') + 
  coord_flip() + 
  theme(axis.text.x = element_text(hjust = 1, angle = 45), axis.title = element_blank()) +
  ggtitle(label = 'DE_ABcom_allDEGs_final')
print(p)
dev.off()


################################# DEG DE_ABcom highlight
# hic magic
source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_hic/Cere-hic-13"
sample <- "Cere-hic-13"
PC <- 36
res <- 1
higashi_out <- "8_higashi_out_higashi_500kb"
work_path_seurat <- "10_seurat_bigSpotShape_higashi_500kb"
embed_type <- "scAB_scale"
setwd(paste0(work_path,  "/", work_path_seurat, "/", embed_type))
load(paste0("Out/", sample, "_Seurat_object_anno.RData"))
load(paste0("Out/", sample, "_Seurat_object_magic.RData"))

# RNA magic
source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_transcriptome/Cere-RNA-6"
sample <- "Cere-RNA-6"
PC <- 20
res <- 2.5
setwd(paste0(work_path, "/6_seurat_bigSpotShape_new"))
load(paste0("Out/", sample, "_cropped_PC",PC,"_res",res,"_anno.RData"))
Seurat_object_cropped_magic <- readRDS(paste0("Out/", sample, "_cropped_magic.rds"))

# highlight
cell_type <- unique(DE_ABcom_allDEGs_final$cell_type)

for(n in 1:length(cell_type)){
  
  # i <- match("Slc17a6", DE_ABcom_allDEGs_final$gene)
  
  DE_ABcom_allDEGs_final_subset <- DE_ABcom_allDEGs_final[DE_ABcom_allDEGs_final$cell_type==cell_type[n],]
  
  
  spatial_p_feature <- list()
  for(i in 1:length(DE_ABcom_allDEGs_final_subset$gene)){
    DEG_spe_select <- DE_ABcom_allDEGs_final_subset$gene[i]
    DEAB_spe_select <- unlist(strsplit(DE_ABcom_allDEGs_final_subset$ABcom[i], ";"))
    
    spatial_p_feature[[DEG_spe_select]] <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = DEG_spe_select, alpha = c(1, 1), ncol = 1, pt.size.factor = 0.8, stroke = NA) + theme(aspect.ratio = 1) + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) + theme(aspect.ratio = 1) #  + theme(legend.position="right")
    spatial_p_feature[[DEG_spe_select]]$layers[[1]]$aes_params=c(spatial_p_feature[[DEG_spe_select]]$layers[[1]]$aes_params, shape=22)
    
    for(j in 1:length(DEAB_spe_select)){
      spatial_p_feature[[DEAB_spe_select[j]]] <- SpatialFeaturePlot(object = Seurat_object_magic, features = DEAB_spe_select[j], alpha = c(1, 1), ncol = 1, pt.size.factor = 0.8, stroke = NA) + scale_fill_gradientn(values=c(0,0.3,0.5,0.6,0.7,0.8,1),colors = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), na.value = "white") + theme(aspect.ratio = 1) # + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) #  + theme(legend.position="right")
      spatial_p_feature[[DEAB_spe_select[j]]]$layers[[1]]$aes_params=c(spatial_p_feature[[DEAB_spe_select[j]]]$layers[[1]]$aes_params, shape=22)
    }
  }
  
  png(file = paste0("Plot/4_findmarker_2D_impose.data_pc36_res1_de_spe_select_RNA_", cell_type[n], ".png"), width = 4*480, height = ceiling(length(spatial_p_feature)/4)*480)
  p <- ggarrange(plotlist = spatial_p_feature, ncol = 4, nrow = ceiling(length(spatial_p_feature)/4))
  print(p)
  dev.off()
}


###################################### ABcom analysis example (Figure 3)

# cluster highlight
source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_hic/Cere-hic-13"
sample <- "Cere-hic-13"
higashi_out <- "8_higashi_out_higashi_500kb"
work_path_seurat <- "10_seurat_bigSpotShape_higashi_500kb"
embed_type <- "scAB_scale"

setwd(paste0(work_path,  "/", work_path_seurat, "/", embed_type))
load(paste0("Out/", sample, "_Seurat_object_anno.RData"))

source("/media/maoni/data/Zhanglin_Junjie/ColorPalettes.R") # 11 cell types
my_color_palette <- c("#E4A5F6", "#CE6CF5", "#BC18B5", "#C61385", "#791367", "#F3153C", "#FDBB15", "#7395BF", "#1380DD","#BE8476", "#584B5F")
scales::show_col(my_color_palette)
names(my_color_palette) <- c("GL_1", "GL_2", "GL_3", "GL_4", "GL_5", "PL_1", "PL_2", "ML_1", "ML_2", "ML_3", "ML_4")

cluster <- c("GL_1|GL_2|GL_3|GL_4|GL_5", "PL_1|PL_2", "ML_1|ML_2")
celltype <- c("GL", "PL", "ML")
spatial_p_cluster <- list()
for(i in 1:length(cluster)){
  Seurat_object_subset <- Seurat_object[ ,grep(cluster[i], Seurat_object$cell_type)]
  spatial_p_cluster[[celltype[i]]] <- SpatialDimPlot(Seurat_object_subset, group.by = "cell_type", crop = FALSE, label = FALSE, label.size = 3, pt.size.factor = 0.8, stroke = NA,  label.box = FALSE, cols = my_color_palette) + theme(aspect.ratio = 1) + ggtitle(celltype[i]) + theme(legend.position="none") + theme(plot.title = element_text(size = 10)) # the spatial plot
  spatial_p_cluster[[celltype[i]]]$layers[[1]]$aes_params=c(spatial_p_cluster[[celltype[i]]]$layers[[1]]$aes_params, shape=22)
}


pdf(file = paste0("Plot/5_Spatial_clusters_PC",PC,"_res",res,"_final_GLPLML.pdf"), width = 3, height = 7)
p_spatial <- spatial_p_cluster[[celltype[1]]]/spatial_p_cluster[[celltype[2]]]/spatial_p_cluster[[celltype[3]]]/spatial_p_cluster[[celltype[4]]]
print(p_spatial)
dev.off()

## Dotplot example
source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_transcriptome/Cere-RNA-6"
sample <- "Cere-RNA-6"
setwd(paste0(work_path, "/6_seurat_bigSpotShape_new"))
PC=20
res=2.5
load(paste0("Out/", sample, "_cropped_PC",PC,"_res",res,"_anno.RData"))
DE_ABcom_allDEGs_final <- read.csv("Plot/DE_ABcom_allDEGs_final.csv", sep=",", header=TRUE)

gene_list_example <- list("GL"=c("Kcnip4", "Kcnd2", "Grik2", "Cntn6", "Car10", "Nav2"),
                          "PL"=c("Grid2", "Grip1", "Ptprm", "Sox5", "Mybpc1"),
                          "ML"=c("Slc24a2", "Asic2"))

Seurat_object_cropped_subset <- Seurat_object_cropped[, grep(pattern = "Granule_1|Granule_2|Granule_3|Granule_4|Purkinje|Bergmann|MLI_1" ,Seurat_object_cropped$cell_type)]
Seurat_object_cropped_subset$cell_type <- factor(Seurat_object_cropped_subset$cell_type, levels = c("Granule_1", "Granule_2", "Granule_3", "Granule_4", "Purkinje", "Bergmann", "MLI_1"))

p_list <- list()
for(i in 1:length(gene_list_example)){
  p_list[[i]] <- DotPlot(Seurat_object_cropped_subset, features = rev(gene_list_example[[i]]), cols = c("grey95", "#B30000"), group.by = 'cell_type') + 
    coord_flip() + 
    theme(axis.text.x = element_blank(), axis.title.y = element_text(size = 10), axis.title = element_blank()) +
    # theme(axis.text.x = element_text(hjust = 1, angle = 45), axis.title = element_blank()) +
    theme(legend.position = "right", legend.key.size = unit(0.3, "cm")) + 
    theme(panel.border = element_rect(color = "black", size = 0.5), legend.text = element_text(size = 10), legend.title = element_text(size = 10))
}

pdf(file = paste0("Plot/5_Spatial_DotPlot_PC",PC,"_res",res,"_final_GLPLML.pdf"), width = 5, height = 7)
p_dotplot <- p_list[[1]]/p_list[[2]]/p_list[[3]]
print(p_dotplot)
dev.off()

source("/media/maoni/data/Zhanglin_Junjie/ColorPalettes.R")
my_color_palette <- c("#E4A5F6", "#CE6CF5", "#BC18B5", "#1312AC", "#F3153C", "#FDBB15", "#7395BF")
scales::show_col(my_color_palette)
names(my_color_palette) <- c("Granule_1", "Granule_2", "Granule_3", "Granule_4", "Purkinje", "Bergmann", "MLI_1")

pdf(file = paste0("Plot/5_Spatial_DotPlot_PC",PC,"_res",res,"_final_GLPLML_color.pdf"), width = 10, height = 1)
scales::show_col(my_color_palette, labels = FALSE, borders = NA, cex_label = 1, ncol = length(my_color_palette))
dev.off()


# DEG DE_ABcom hightlight example 
source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_transcriptome/Cere-RNA-6"
sample <- "Cere-RNA-6"
setwd(paste0(work_path, "/6_seurat_bigSpotShape_new"))

DE_ABcom_allDEGs_final <- read.csv("Plot/DE_ABcom_allDEGs_final.csv", sep=",", header=TRUE)

gene_example_final <- c("Grik2", "Ptprm", "Slc24a2")

DE_ABcom_allDEGs_final <- DE_ABcom_allDEGs_final[match(gene_example_final, DE_ABcom_allDEGs_final$gene), ]

spatial_p_feature <- list()
for(i in 1:length(DE_ABcom_allDEGs_final$gene)){
  DEG_spe_select <- DE_ABcom_allDEGs_final$gene[i]
  DEAB_spe_select <- unlist(strsplit(DE_ABcom_allDEGs_final$ABcom[i], ";"))[[1]]
  
  spatial_p_feature[[DEG_spe_select]] <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = DEG_spe_select, alpha = c(1, 1), ncol = 1, pt.size.factor = 0.8, stroke = NA) + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) + theme(aspect.ratio = 1)
  spatial_p_feature[[DEG_spe_select]]$layers[[1]]$aes_params=c(spatial_p_feature[[DEG_spe_select]]$layers[[1]]$aes_params, shape=22)
  
  spatial_p_feature[[DEAB_spe_select]] <- SpatialFeaturePlot(object = Seurat_object_magic, features = DEAB_spe_select, alpha = c(1, 1), ncol = 1, pt.size.factor = 0.8, stroke = NA) + scale_fill_gradientn(values=c(0,0.3,0.5,0.6,0.7,0.8,1),colors = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), na.value = "white") + theme(aspect.ratio = 1)
  spatial_p_feature[[DEAB_spe_select]]$layers[[1]]$aes_params=c(spatial_p_feature[[DEAB_spe_select]]$layers[[1]]$aes_params, shape=22)
}

pdf(file = paste0("Plot/5_Spatial_PC",PC,"_res",res,"_final_GLPLML_DEG_DEAB_example_final.pdf"), width = 2*7, height = 4*7)
p <- ggarrange(plotlist = spatial_p_feature, ncol = 2, nrow = 4)
print(p)
dev.off()


### gene_example_track_region
setwd("/media/maoni/data/CZP/Figures/ABcom")
gene_list_example <- list("GL"=c("Kcnip4", "Kcnd2", "Grik2", "Cntn6", "Car10", "Nav2"),
                          "PL"=c("Grid2", "Grip1", "Ptprm", "Sox5", "Mybpc1"),
                          "ML"=c("Slc24a2", "Asic2"))
PCG_anno <- read.table("/media/maoni/data/Reference/mouse/GRCm38_mm10/ensembl/All_gene_Anno_PCG.bed", sep="\t", header=FALSE)
gene_list_example_anno <- PCG_anno[match(unlist(gene_list_example), PCG_anno$V4), ]
gene_list_example_anno$middle <- (gene_list_example_anno$V3 + gene_list_example_anno$V2)/2
gene_list_example_anno$middle_floor <- floor(gene_list_example_anno$middle / 100000) * 100000
gene_list_example_anno$region_start <- gene_list_example_anno$middle_floor-3500000
gene_list_example_anno$region_end <- gene_list_example_anno$middle_floor+3500000
gene_list_example_anno$region <- apply(gene_list_example_anno, 1, function(x) paste(x[1], paste(x[11], x[12], sep = "-"), sep = ":"))
gene_list_example_anno$region <- gsub(" ", "", gene_list_example_anno$region)
write.table(gene_list_example_anno, "Cere_gene_list_example_final.txt", sep = "\t", quote = FALSE,row.names = FALSE, col.names = FALSE)


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
work_path <- "/media/maoni/data/CZP/spatial_transcriptome/Cere-RNA-6"
sample <- "Cere-RNA-6"
setwd(paste0(work_path, "/6_seurat_bigSpotShape_new"))

DE_ABcom_allDEGs_final <- read.csv("Plot/DE_ABcom_allDEGs_final.csv", sep=",", header=TRUE)

DE_ABcom_allDEGs_final$class <- rep(NA, length(DE_ABcom_allDEGs_final$cell_type))
DE_ABcom_allDEGs_final$class[grep("Granule_1|Granule_2|Granule_3|Granule_4", DE_ABcom_allDEGs_final$cell_type)] <- "GL"
DE_ABcom_allDEGs_final$class[grep("Purkinje|Bergmann", DE_ABcom_allDEGs_final$cell_type)] <- "PL"
DE_ABcom_allDEGs_final$class[grep("MLI_1|MLI_2", DE_ABcom_allDEGs_final$cell_type)] <- "ML"

PCG_list <- list()
for (organ in c("GL","PL","ML")){
  # organ <- "GL"
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

my_color <- c(brewer.pal(9, "Set1")[c(6,8,7)])
names(my_color) <- c("GL","PL","ML")

GO_barplot_all <- c()
GO_barplot_list <- list()
for (organ in c("GL","PL","ML")){
  # organ <- "Radial_glia"
  go_enrich_BP <- read.csv(paste0(organ, "_DE_ABcom_allDEG_enrichGO_GOBP.txt"), sep="\t", header=T)
  GO_barplot <- data.frame(go_enrich_BP[1:10, ][, c("Description", "pvalue", "Count")], organ)
  GO_barplot$pvalue <- -log10(GO_barplot$pvalue)
  colnames(GO_barplot) <- c("GO_term", "pvalue", "Count", "organ")
  GO_barplot$GO_term <- factor(GO_barplot$GO_term, levels = rev(GO_barplot$GO_term))
  GO_barplot_all <- rbind(GO_barplot_all, GO_barplot)
}

GO_barplot_all$organ <- factor(GO_barplot_all$organ, levels = c("GL","PL","ML"))
p <- ggplot(GO_barplot_all, aes(x=GO_term,y=pvalue, fill=organ)) + geom_bar(stat="identity", width=0.8) + coord_flip() +  xlab("GOBP") + ylab("-log10(pvalue)") +
  facet_wrap(~organ, scales = "free_y") +
  theme_bw() + 
  theme(legend.position="none") + 
  theme(panel.grid = element_blank()) + 
  scale_fill_manual(values=my_color)

pdf("Plot/DE_ABcom_allDEG_enrichGO_GOBP_barplot_GLPLML.pdf", width = 25, height = 3)
print(p)
dev.off()

###########################################################################################################################################
############################################################## GL subset analysis #########################################################
###########################################################################################################################################
# RNA
source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_transcriptome/Cere-RNA-6"
sample <- "Cere-RNA-6"
PC <- 20
res <- 2.5
setwd(paste0(work_path, "/6_seurat_bigSpotShape_new"))
load(paste0("Out/", sample, "_cropped_PC",PC,"_res",res,"_anno.RData"))
Seurat_object_cropped_magic <- readRDS(paste0("Out/", sample, "_cropped_magic.rds"))

# hic
source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_hic/Cere-hic-13"
sample <- "Cere-hic-13"
PC <- 36
res <- 1
higashi_out <- "8_higashi_out_higashi_500kb"
work_path_seurat <- "10_seurat_bigSpotShape_higashi_500kb"
embed_type <- "scAB_scale"
setwd(paste0(work_path,  "/", work_path_seurat, "/", embed_type))
load(paste0("Out/", sample, "_Seurat_object_anno.RData"))
load(paste0("Out/", sample, "_Seurat_object_magic.RData"))

# cluster highlight
cluster <- c("GL_1|GL_2", "GL_3|GL_4", "GL_5")
celltype <- c("GL_12", "GL_34", "GL_5")
spatial_p_cluster <- list()
for(i in 1:length(cluster)){
  Seurat_object_subset <- Seurat_object[ ,grep(cluster[i], Seurat_object$cell_type)]
  spatial_p_cluster[[celltype[i]]] <- SpatialDimPlot(Seurat_object_subset, group.by = "cell_type", crop = FALSE, label = FALSE, 
                                                     label.size = 3, pt.size.factor = 0.8, stroke = NA,  label.box = FALSE, cols = my_color_palette) + 
    theme(aspect.ratio = 1) + ggtitle(celltype[i]) + theme(legend.position="none") + theme(plot.title = element_text(size = 10)) # the spatial plot
  spatial_p_cluster[[celltype[i]]]$layers[[1]]$aes_params=c(spatial_p_cluster[[celltype[i]]]$layers[[1]]$aes_params, shape=22)
}

pdf(file = paste0("Plot/5_Spatial_clusters_PC",PC,"_res",res,"_final_GL12_GL34_GL5.pdf"), width = 3, height = 8)
p_spatial <- spatial_p_cluster[[celltype[1]]]/spatial_p_cluster[[celltype[2]]]/spatial_p_cluster[[celltype[3]]]
print(p_spatial)
dev.off()

Seurat_object_subset <- Seurat_object[, grep("GL_1|GL_2|GL_3|GL_4|GL_5", Seurat_object$cell_type)]
Seurat_object_cropped_subset <- Seurat_object_cropped[, grep("Granule_1|Granule_2|Granule_3", Seurat_object_cropped$cell_type)]

# DE_ABcom_allDEGs_final <- read.csv("Plot/DE_ABcom_allDEGs_final.csv", sep=",", header=TRUE)
# gene_example_final <- c("Kcnip4", "Gprin3", "Galntl6", "Kcnd2")
# DE_ABcom_allDEGs_final_example <- DE_ABcom_allDEGs_final[match(gene_example_final, DE_ABcom_allDEGs_final$gene), ]

DEG_spe_select_all <- c("Kcnip4", "Kcnd2", "Dpp6", "Gprin3", "Galntl6")
DEAB_spe_select_all <- c("chr5:48500000-49000000", "chr6:21000000-21500000", "chr5:26500000-27000000","chr6:59000000-59500000", "chr8:58500000-59000000")
DEG_spe_select_box_RNA_all <- c("Kcnip4-RNA", "Kcnd2-RNA", "Dpp6-RNA", "Gprin3-RNA", "Galntl6-RNA")
DEAB_spe_select_box_hic_all <- c("chr5:48500000-49000000-hic", "chr6:21000000-21500000-hic", "chr5:26500000-27000000-hic","chr6:59000000-59500000-hic", "chr8:58500000-59000000-hic")

# ggboxplot (pre data for hic)
cluster = c("GL_12", "GL_34", "GL_5")
Cluster_spe_region_Normalized_ABscore <- c()
for(i in 1:length(DEAB_spe_select_all)){
  geneExpression <- data.frame(DEAB_spe_select_all[i], Seurat_object_subset@assays$Spatial@data[which(rownames(Seurat_object_subset)==DEAB_spe_select_all[i]),], Seurat_object_subset@meta.data$cell_type)
  colnames(geneExpression) <- c("group", "Normalized_ABscore", "cell_type")
  Cluster_spe_region_Normalized_ABscore <- rbind(Cluster_spe_region_Normalized_ABscore, geneExpression)
}

Cluster_spe_region_Normalized_ABscore$cell_type <- as.character(Cluster_spe_region_Normalized_ABscore$cell_type)

Cluster_spe_region_Normalized_ABscore[grep("GL_1|GL_2", Cluster_spe_region_Normalized_ABscore$cell_type),3] <- "GL_12"
Cluster_spe_region_Normalized_ABscore[grep("GL_3|GL_4", Cluster_spe_region_Normalized_ABscore$cell_type),3] <- "GL_34"
Cluster_spe_region_Normalized_ABscore$cell_type <- factor(Cluster_spe_region_Normalized_ABscore$cell_type, levels = cluster)
Cluster_spe_region_Normalized_ABscore$group <- factor(Cluster_spe_region_Normalized_ABscore$group, levels = DEAB_spe_select_all)

# ggboxplot (pre data for RNA)
Cluster_spe_region_Normalized_Expression <- c()
for(i in 1:length(DEG_spe_select_all)){
  geneExpression <- data.frame(DEG_spe_select_all[i], Seurat_object_cropped_subset@assays$SCT@data[which(rownames(Seurat_object_cropped_subset)==DEG_spe_select_all[i]),], Seurat_object_cropped_subset@meta.data$cell_type)
  colnames(geneExpression) <- c("group", "Normalized_Expression", "cell_type")
  Cluster_spe_region_Normalized_Expression <- rbind(Cluster_spe_region_Normalized_Expression, geneExpression)
}
Cluster_spe_region_Normalized_Expression$group <- factor(Cluster_spe_region_Normalized_Expression$group, levels = DEG_spe_select_all)

# plot
color_exp <- list(Kcnip4 = c(1:6, 10),
                  Kcnd2 = c(1:6, 9,10),
                  Dpp6 = c(1:6, 9,10),
                  Gprin3 = c(1:6, 10:11),
                  Galntl6 = c(1:6, 10:11))

color_value_AB <- list(Kcnip4 = c(0,0.3,0.5,0.6,0.7,0.8,1),
                       Kcnd2 = c(0,0.3,0.5,0.6,0.7,0.8,1),
                       Dpp6 = c(0,0.3,0.5,0.6,0.7,0.8,1),
                       Gprin3 = c(0,0.3,0.5,0.6,0.7,0.8,1),
                       Galntl6 = c(0,0.3,0.5,0.6,0.7,0.8,1))

spatial_p_feature <- list()

for(i in 1:length(DEG_spe_select_all)){
  
  # i <- 2
  
  DEG_spe_select <- DEG_spe_select_all[i]
  DEAB_spe_select <- DEAB_spe_select_all[i]
  DEAB_spe_select_box_hic <- DEAB_spe_select_box_hic_all[i]
  DEG_spe_select_box_RNA <- DEG_spe_select_box_RNA_all[i]
  
  ## AB region
  spatial_p_feature[[DEAB_spe_select]] <- SpatialFeaturePlot(object = Seurat_object_magic, features = DEAB_spe_select, alpha = c(1, 1), ncol = 1, pt.size.factor = 0.8, stroke = NA) + 
    scale_fill_gradientn(values=color_value_AB[[i]],colors = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), na.value = "white") + theme(aspect.ratio = 1) + 
    theme(legend.position = "top")
  spatial_p_feature[[DEAB_spe_select]]$layers[[1]]$aes_params=c(spatial_p_feature[[DEAB_spe_select]]$layers[[1]]$aes_params, shape=22)
  
  ## gene expression
  spatial_p_feature[[DEG_spe_select]] <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = DEG_spe_select, alpha = c(1, 1), ncol = 1, pt.size.factor = 0.8, stroke = NA) + 
    scale_fill_gradientn(colours = olcol <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))[color_exp[[i]]]) + theme(aspect.ratio = 1) + 
    theme(legend.position = "top")
  spatial_p_feature[[DEG_spe_select]]$layers[[1]]$aes_params=c(spatial_p_feature[[DEG_spe_select]]$layers[[1]]$aes_params, shape=22)
  
  
  ## compare AB region
  Cluster_spe_region_Normalized_ABscore_spe <- Cluster_spe_region_Normalized_ABscore[Cluster_spe_region_Normalized_ABscore$group==DEAB_spe_select_all[i],]
  
  stat.test <- Cluster_spe_region_Normalized_ABscore_spe %>%
    group_by(group) %>%
    wilcox_test(Normalized_ABscore ~ cell_type)
  stat.test
  
  stat.test <- stat.test %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  stat.test <- stat.test %>% add_y_position()
  
  my_color_hic <- ann_colors$cluster[c(1,3,5)]
  names(my_color_hic) <- c("GL_12", "GL_34", "GL_5")
  
  spatial_p_feature[[DEAB_spe_select_box_hic]] <- ggboxplot(Cluster_spe_region_Normalized_ABscore_spe, x = "cell_type", y = "Normalized_ABscore", fill = "cell_type", facet.by = "group", 
                  ncol=1, font.label = list(size = 3, color = "black"), outlier.shape = NA) +
    theme(legend.position = "top") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + 
    stat_summary(fun.y = median, geom = "line", aes(group = 1), color = "black") +
    stat_summary(fun.y = median, geom = "point", color = "black", size = 2) + 
    scale_fill_manual(values=my_color_hic) + 
    theme(aspect.ratio = 1)
  
  ## compare RNA
  Cluster_spe_region_Normalized_Expression_spe <- Cluster_spe_region_Normalized_Expression[Cluster_spe_region_Normalized_Expression$group==DEG_spe_select_all[i],]
  
  stat.test <- Cluster_spe_region_Normalized_Expression_spe %>%
    group_by(group) %>%
    wilcox_test(Normalized_Expression ~ cell_type)
  stat.test
  
  stat.test <- stat.test %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  stat.test <- stat.test %>% add_y_position()
  
  my_color_RNA <- c("#E4A5F6", "#CE6CF5", "#BC18B5", "#1312AC", "#F3153C", "#FDBB15", "#7395BF")
  names(my_color_RNA) <- c("Granule_1", "Granule_2", "Granule_3", "Granule_4", "Purkinje", "Bergmann", "MLI_1")
  
  spatial_p_feature[[DEG_spe_select_box_RNA]] <- ggboxplot(Cluster_spe_region_Normalized_Expression_spe, x = "cell_type", y = "Normalized_Expression", fill = "cell_type", facet.by = "group", 
                                                           ncol=1, font.label = list(size = 3, color = "black"), outlier.shape = NA) +
    theme(legend.position = "top") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + 
    stat_summary(fun.y = median, geom = "line", aes(group = 1), color = "black") +
    stat_summary(fun.y = median, geom = "point", color = "black", size = 2) + 
    scale_fill_manual(values=my_color_RNA) + 
    theme(aspect.ratio = 1)
}

pdf(file = paste0("Plot/5_Spatial_PC",PC,"_res",res,"_final_GL_DEG_DEAB_example_final_5_7.pdf"), width = 4*5, height = 5*5)
p1 <- ggarrange(plotlist = spatial_p_feature, ncol = 4, nrow = 5)
print(p1)
dev.off()



################################################### DE_ABcom scale
# mean of each cell type
dim(Seurat_object@assays$Spatial@counts)
dim(de_markers_13)
count_matrix <- Seurat_object@assays$Spatial@counts

de_markers_13_GL <- de_markers_13[grep("GL_1|GL_2|GL_3|GL_4|GL_5", de_markers_13$cell_type), ]

length(unique(de_markers_13_GL$gene))  # 186
length(unique(de_markers_13_GL$gene[grep("GL_1", de_markers_13_GL$cell_type)]))  # 88
length(unique(de_markers_13_GL$gene[grep("GL_2", de_markers_13_GL$cell_type)])) # 105
length(unique(de_markers_13_GL$gene[grep("GL_3", de_markers_13_GL$cell_type)])) # 89
length(unique(de_markers_13_GL$gene[grep("GL_4", de_markers_13_GL$cell_type)])) # 80
length(unique(de_markers_13_GL$gene[grep("GL_5", de_markers_13_GL$cell_type)])) # 16

count_matrix_de_markers <- count_matrix[match(de_markers_13_GL$gene, rownames(count_matrix)),]
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

temp <- as.data.frame(count_matrix_de_markers_norm_mean_trans[, 1:5])
temp <- temp[, c("GL_1", "GL_2", "GL_3", "GL_4", "GL_5")]

library(pheatmap)

annotation_col = data.frame(
  cluster = factor(c("GL_1", "GL_2", "GL_3", "GL_4", "GL_5"))
)
rownames(annotation_col) = c("GL_1", "GL_2", "GL_3", "GL_4", "GL_5")

ann_colors = list(
  cluster = c(GL_1 = "#E4A5F6", GL_2 = "#CE6CF5", GL_3 = "#BC18B5", GL_4 = "#C61385", GL_5 = "#791367")
)


########################################################## Venn plot    
# 读取数据文件
library(venn)
venn_list <- list(de_markers_13_GL$gene[de_markers_13_GL$cell_type=="GL_1"], de_markers_13_GL$gene[de_markers_13_GL$cell_type=="GL_2"], de_markers_13_GL$gene[de_markers_13_GL$cell_type=="GL_3"], de_markers_13_GL$gene[de_markers_13_GL$cell_type=="GL_4"], de_markers_13_GL$gene[de_markers_13_GL$cell_type=="GL_5"])   # 制作韦恩图搜所需要的列表文件
names(venn_list) <- c("GL_1", "GL_2", "GL_3", "GL_4", "GL_5")   # 把列名赋值给列表的key值
venn_list = purrr::map(venn_list,na.omit)      # 删除列表中每个向量中的NA

mycolor <- ann_colors$cluster

pdf("Plot/4_de_markers_13_GL_vennplot.pdf",width=7,heigh=7)
venn(venn_list,
     zcolor = mycolor, # 调整颜色，style是默认颜色的，bw是无颜色，当然也可以自定义颜色
     opacity = 0.5,  # 调整颜色透明度
     box = F,        # 是否添加边框
     ilcs = 2,     # 数字大小
     sncs = 2        # 组名字体大小
)
dev.off()


# 查看交集详情,并导出结果
inter <- get.venn.partitions(venn_list)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = '|')
inter <- subset(inter, select = -..values.. )
inter <- subset(inter, select = -..set.. )
write.table(inter, "Plot/4_de_markers_13_GL_vennplot_intersect.csv", row.names = FALSE, sep = ',', quote = FALSE)

##################################################### Acom_spe
GL_1_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==1), which(apply(temp, 1, max) > 0.5))]
GL_2_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==2), which(apply(temp, 1, max) > 0.5))]
GL_3_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==3), which(apply(temp, 1, max) > 0.5))]
GL_4_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==4), which(apply(temp, 1, max) > 0.5))]
GL_5_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==5), which(apply(temp, 1, max) > 0.5))]

# GL_1_spe <- unlist(strsplit(inter$values[31], "\\|"))
# GL_2_spe <- unlist(strsplit(inter$values[30], "\\|"))
# GL_3_spe <- unlist(strsplit(inter$values[28], "\\|"))
# GL_4_spe <- unlist(strsplit(inter$values[24], "\\|"))
# GL_5_spe <- unlist(strsplit(inter$values[16], "\\|"))

GL_1_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==1), which(apply(temp, 1, max) > 0))]
GL_2_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==2), which(apply(temp, 1, max) > 0))]
GL_3_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==3), which(apply(temp, 1, max) > 0))]
GL_4_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==4), which(apply(temp, 1, max) > 0))]
GL_5_spe <- rownames(temp)[intersect(which(apply(temp, 1, which.max)==5), which(apply(temp, 1, max) > 0))]

length(GL_1_spe) # 7
length(GL_2_spe) # 48
length(GL_3_spe) # 23
length(GL_4_spe) # 29
length(GL_5_spe) # 9

temp_Acom_spe <- c(GL_1_spe, GL_2_spe, GL_3_spe, GL_4_spe, GL_5_spe) # 163

temp_Acom_spe_matrix12 <- count_matrix_de_markers_norm_mean_trans[match(c(GL_1_spe, GL_2_spe), rownames(count_matrix_de_markers_norm_mean_trans)), ][,1:5]
temp_Acom_spe_matrix12 <- data.frame(rowMeans(temp_Acom_spe_matrix12[, c("GL_1", "GL_2")]), rowMeans(temp_Acom_spe_matrix12[, c("GL_3", "GL_4")]), temp_Acom_spe_matrix12[, c("GL_5")])
temp_Acom_spe_matrix12 <- temp_Acom_spe_matrix12[order(temp_Acom_spe_matrix12[,1], decreasing = TRUE), ]
colnames(temp_Acom_spe_matrix12) <- c("GL_12", "GL_34", "GL_5")

temp_Acom_spe_matrix34 <- count_matrix_de_markers_norm_mean_trans[match(c(GL_3_spe, GL_4_spe), rownames(count_matrix_de_markers_norm_mean_trans)), ][,1:5]
temp_Acom_spe_matrix34 <- data.frame(rowMeans(temp_Acom_spe_matrix34[, c("GL_1", "GL_2")]), rowMeans(temp_Acom_spe_matrix34[, c("GL_3", "GL_4")]), temp_Acom_spe_matrix34[, c("GL_5")])
temp_Acom_spe_matrix34 <- temp_Acom_spe_matrix34[order(temp_Acom_spe_matrix34[,2], decreasing = TRUE), ]
colnames(temp_Acom_spe_matrix34) <- c("GL_12", "GL_34", "GL_5")

temp_Acom_spe_matrix5 <- count_matrix_de_markers_norm_mean_trans[match(c(GL_5_spe), rownames(count_matrix_de_markers_norm_mean_trans)), ][,1:5]
temp_Acom_spe_matrix5 <- data.frame(rowMeans(temp_Acom_spe_matrix5[, c("GL_1", "GL_2")]), rowMeans(temp_Acom_spe_matrix5[, c("GL_3", "GL_4")]), temp_Acom_spe_matrix5[, c("GL_5")])
temp_Acom_spe_matrix5 <- temp_Acom_spe_matrix5[order(temp_Acom_spe_matrix5[,3], decreasing = TRUE), ]
colnames(temp_Acom_spe_matrix5) <- c("GL_12", "GL_34", "GL_5")

dim(temp_Acom_spe_matrix12) # 87
dim(temp_Acom_spe_matrix34) # 67
dim(temp_Acom_spe_matrix5) # 9

temp_Acom_spe_matrix <- rbind(temp_Acom_spe_matrix12[1:20,], temp_Acom_spe_matrix34[1:30,], temp_Acom_spe_matrix5)

annotation_col = data.frame(
  cluster = factor(c("GL_12", "GL_34", "GL_5"))
)
rownames(annotation_col) = c("GL_12", "GL_34", "GL_5")

ann_colors = list(
  cluster = c(GL_12 = "#E4A5F6", GL_34 = "#BC18B5", GL_5 = "#791367")
)

pdf(file = paste("Plot/4_allDE_Acom_heat_anno_average_GL_spe.pdf"), width = 4, height = 4)
# pheatmap(temp_Acom_spe_matrix, border = F, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, color = colorRampPalette(c("#C44C9C", "#B96397", "#AA7591", "#968686", "#86917D", "#739D6F", "#5CA65F", "#1CB24B"))(100))
pheatmap(temp_Acom_spe_matrix, border = F, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, color = colorRampPalette(c("#36359D", "#0F74DC", "#05A5C6", "#4EBD90", "#A9BE69", "#F3B947", "#F4E61A")[c(1,2,3,4,6,7)])(100))
dev.off()

DE_ABcom_GL <- data.frame(temp_Acom_spe, c(rep("GL_1", length(GL_1_spe)), rep("GL_2", length(GL_2_spe)), rep("GL_3", length(GL_3_spe)), rep("GL_4", length(GL_4_spe)), rep("GL_5", length(GL_5_spe))))
colnames(DE_ABcom_GL) <- c("gene", "cell_type")

cell_type <- names(table(DE_ABcom_GL$cell_type))
for(i in 1:length(table(DE_ABcom_GL$cell_type))){
  png(paste0("Plot/4_findmarker_2D_impose.data_pc36_res1_magic_", cell_type[i], ".png"),height = table(DE_ABcom_GL$cell_type)[i]/5*480, width = 5*480)
  plot <- SpatialFeaturePlot(object = Seurat_object_magic, features = DE_ABcom_GL$gene[DE_ABcom_GL$cell_type==cell_type[i]], alpha = c(1, 1), ncol = 5, pt.size.factor = 0.8) # + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")))
  plot$layers[[1]]$aes_params=c(plot$layers[[1]]$aes_params, shape=22)
  print(plot)
  dev.off()
}


## top10_anno_vlnplot
Seurat_object_subset <- Seurat_object[, grep("GL_1|GL_2|GL_3|GL_4|GL_5", Seurat_object$cell_type)]

# ggviolin(five cluster)
Seurat_object_subset <- AddModuleScore(Seurat_object_subset,
                                features = list(GL_1=top10_anno$gene[top10_anno$cell_type=="GL_1"], GL_2=top10_anno$gene[top10_anno$cell_type=="GL_2"], GL_3=top10_anno$gene[top10_anno$cell_type=="GL_3"], GL_4=top10_anno$gene[top10_anno$cell_type=="GL_4"], GL_5=top10_anno$gene[top10_anno$cell_type=="GL_5"]),
                                name = c("GL_1", "GL_2", "GL_3", "GL_4", "GL_5"),
                                assay = "Spatial")

modules <- c("GL_11", "GL_22", "GL_33", "GL_44", "GL_55")
cluster = c("GL_1", "GL_2", "GL_3", "GL_4", "GL_5")
Cluster_spe_region_module_score <- c()
for(i in 1:5){
  module_score <- data.frame(cluster[i], Seurat_object_subset@meta.data[,c(modules[i], "cell_type")])
  colnames(module_score) <- c("group", "AverageExpression", "cell_type")
  Cluster_spe_region_module_score <- rbind(Cluster_spe_region_module_score, module_score)
}

stat.test <- Cluster_spe_region_module_score %>%
  group_by(group) %>%
  t_test(AverageExpression ~ cell_type)
stat.test

stat.test <- stat.test %>% add_y_position()

pdf(paste0("Plot/4_top_anno_GL_VlnPlot.pdf"), height = 4, width = 15)
ggviolin(Cluster_spe_region_module_score, x = "cell_type", y = "AverageExpression", fill = "cell_type", facet.by = "group", 
         ncol = 5, font.label = list(size = 3, color = "black"),
         add = "boxplot",
         add.params = list(fill = "white", width = 0.05,linetype = 1)) +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + 
  scale_fill_manual(values=ann_colors$cluster)
dev.off()

# ggviolin (three cluster)
Seurat_object_subset <- AddModuleScore(Seurat_object_subset,
                                       features = list(GL_12=top10_anno$gene[grep("GL_1|GL_2", top10_anno$cell_type)], GL_34=top10_anno$gene[grep("GL_3|GL_4", top10_anno$cell_type)], GL_5=top10_anno$gene[top10_anno$cell_type=="GL_5"]),
                                       name = c("GL_12", "GL_34", "GL_5"),
                                       assay = "Spatial")

modules <- c("GL_121", "GL_342", "GL_53")
cluster = c("GL_12", "GL_34", "GL_5")
Cluster_spe_region_module_score <- c()
for(i in 1:3){
  module_score <- data.frame(cluster[i], Seurat_object_subset@meta.data[,c(modules[i], "cell_type")])
  colnames(module_score) <- c("group", "AverageExpression", "cell_type")
  Cluster_spe_region_module_score <- rbind(Cluster_spe_region_module_score, module_score)
}

Cluster_spe_region_module_score$cell_type <- as.character(Cluster_spe_region_module_score$cell_type)

Cluster_spe_region_module_score[grep("GL_1|GL_2", Cluster_spe_region_module_score$cell_type),3] <- "GL_12"
Cluster_spe_region_module_score[grep("GL_3|GL_4", Cluster_spe_region_module_score$cell_type),3] <- "GL_34"
Cluster_spe_region_module_score$cell_type <- factor(Cluster_spe_region_module_score$cell_type, levels = cluster)

stat.test <- Cluster_spe_region_module_score %>%
  group_by(group) %>%
  t_test(AverageExpression ~ cell_type)
stat.test

stat.test <- stat.test %>% add_y_position()

my_color <- ann_colors$cluster[c(1,3,5)]
names(my_color) <- c("GL_12", "GL_34", "GL_5")

pdf(paste0("Plot/4_top_anno_GL_VlnPlot_2.pdf"), height = 8, width = 3)
ggviolin(Cluster_spe_region_module_score, x = "cell_type", y = "AverageExpression", fill = "cell_type", facet.by = "group", 
         ncol = 1, font.label = list(size = 3, color = "black"),
         add = "boxplot",
         add.params = list(fill = "white", width = 0.05,linetype = 1)) +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + 
  scale_fill_manual(values=my_color)
dev.off()



# DotPlot
pdf(paste0("Plot/4_top10_anno_GL_DotPlot.pdf"), height = 5, width = 15)
p <- DotPlot(Seurat_object_subset, features = unique(top10_anno$gene[grep("GL_1|GL_2|GL_3|GL_4|GL_5", top10_anno$cell_type)]), cols = c("grey", "red"), group.by = 'cell_type') + 
  # coord_flip() + 
  scale_color_viridis() +
  theme(axis.text.x = element_text(hjust = 1, angle = 45), axis.title = element_blank()) +
  ggtitle(label = 'top10_anno_GL_DotPlot')
print(p)
dev.off()

top10_anno_select <- top10_anno[c(1,7,11,13,16,23,24,35,36,38,41,42,46,49,50),]
write.table(top10_anno_select, "/media/maoni/data/CZP/Figures/ABcom/Cere_PC1_example/GL_spe/top10_anno_select.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# DE_ABcom_select
top10_anno_select_final <- top10_anno_select[c(1,8,13), ]
spatial_p_feature <- list()
for(i in 1:length(top10_anno_select_final$gene)){
  DEAB_spe_select <- top10_anno_select_final$gene[i]

  spatial_p_feature[[DEAB_spe_select]] <- SpatialFeaturePlot(object = Seurat_object_magic, features = DEAB_spe_select, alpha = c(1, 1), ncol = 1, pt.size.factor = 0.8, stroke = NA) + scale_fill_gradientn(values=c(0,0.3,0.5,0.6,0.7,0.8,1),colors = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), na.value = "white") + theme(aspect.ratio = 1)
  spatial_p_feature[[DEAB_spe_select]]$layers[[1]]$aes_params=c(spatial_p_feature[[DEAB_spe_select]]$layers[[1]]$aes_params, shape=22)
}
pdf(file = paste0("Plot/5_Spatial_PC",PC,"_res",res,"_final_GL_DEAB_example_final.pdf"), width = 3*7, height = 1*7)
p <- ggarrange(plotlist = spatial_p_feature, ncol = 3, nrow = 1)
print(p)
dev.off()


# DE_ABcom_DEGs_select_gene
DE_ABcom_allDEGs_final_GL_gene <- c("Kcnip4", "Gprin3", "Galntl6", "Kcnd2", "Rnf152", "Syt1", "Gpr158", "Cadm2", "Ralyl")
DE_ABcom_allDEGs_final_GL_gene <- c("Kcnip4", "Gprin3", "Galntl6")
spatial_p_feature <- list()
for(i in 1:length(DE_ABcom_allDEGs_final_GL_gene)){
  DEAB_spe_select <- DE_ABcom_allDEGs_final_GL_gene[i]
  
  spatial_p_feature[[DEAB_spe_select]] <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = DEAB_spe_select, alpha = c(1, 1), ncol = 1, pt.size.factor = 0.8, stroke = NA) + 
    # scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))[c(1:5, 9,10,11)]) + theme(aspect.ratio = 1)
    scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))[c(1:5, 11)]) + theme(aspect.ratio = 1)
  spatial_p_feature[[DEAB_spe_select]]$layers[[1]]$aes_params=c(spatial_p_feature[[DEAB_spe_select]]$layers[[1]]$aes_params, shape=22)
}
pdf(file = paste0("Plot/5_Spatial_PC",PC,"_res",res,"_final_GL_DEAB_DEGs_example_final_2.pdf"), width = 5*7, height = 2*7)
p <- ggarrange(plotlist = spatial_p_feature, ncol = 5, nrow = 2)
print(p)
dev.off()

# DE_ABcom_DEGs_select_region
DE_ABcom_allDEGs_final <- read.csv("Plot/DE_ABcom_allDEGs_final.csv", sep=",", header=TRUE)
DE_ABcom_allDEGs_final_GL_region <- DE_ABcom_allDEGs_final$ABcom[match(DE_ABcom_allDEGs_final_GL_gene, DE_ABcom_allDEGs_final$gene)]
DE_ABcom_allDEGs_final_GL_region <- c("chr5:48000000-48500000;chr5:48500000-49000000", DE_ABcom_allDEGs_final_GL_region)

DE_ABcom_allDEGs_final_GL_region_list_final <- data.frame(DE_ABcom_allDEGs_final_GL_gene, DE_ABcom_allDEGs_final_GL_region_list)
write.table(DE_ABcom_allDEGs_final_GL_region_list_final, "/media/maoni/data/CZP/Figures/ABcom/Cere_PC1_example/GL_spe/DE_ABcom_allDEGs_final_GL_region_list_final.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

######################################################## DE_Acom_GL fuctional analysis
library(clusterProfiler)
library(org.Mm.eg.db)
library(pathview)
library(enrichplot)
library(ggplot2)
library(ggpubr)
library(stringi)

# de_markers_13_GL <- de_markers_13[grep("GL_1|GL_2|GL_3|GL_4|GL_5", de_markers_13$cell_type), ]
# DE_ABcom_region_PCGs <- read.table("Plot/DE_ABcom_region_PCGs.bed", sep="\t", header=FALSE)
# Acom_GL <- intersect(de_markers_13_GL$gene, DE_ABcom_region_PCGs$V4)
# DE_Acom_GL_PCGs <- DE_ABcom_region_PCGs[which(DE_ABcom_region_PCGs$V4 %in% Acom_GL), ]
# DE_Acom_GL_PCGs_cell_type <- data.frame(DE_Acom_GL_PCGs$V4, de_markers_13_GL$cell_type[match(DE_Acom_GL_PCGs$V4, de_markers_13_GL$gene)], DE_Acom_GL_PCGs$V8)
# 
# colnames(DE_Acom_GL_PCGs_cell_type) <- c("DE_Acom", "cell_type", "gene")
# head(DE_Acom_GL_PCGs_cell_type)
# 
# write.table(DE_Acom_GL_PCGs_cell_type, "DE_Acom_GL_PCGs_cell_type.txt", sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
# 
# DE_Acom_GL_PCGs_cell_type <- read.table("DE_Acom_GL_PCGs_cell_type.txt", sep="\t", header=TRUE)

DE_Acom_GL <- read.table("Plot/DE_Acom_GL.txt", sep="\t", header = TRUE)
DE_ABcom_region_PCGs <- read.table("Plot/DE_ABcom_region_PCGs.bed", sep="\t", header=FALSE)
DE_Acom_GL_PCGs <- data.frame(DE_Acom_GL, DE_ABcom_region_PCGs$V8[match(DE_Acom_GL$gene, DE_ABcom_region_PCGs$V4)])
colnames(DE_Acom_GL_PCGs) <- c("DE_Acom", "cell_type", "gene")
head(DE_Acom_GL_PCGs)

PCG_list <- list()
for (organ in c("GL_1|GL_2","GL_3|GL_4","GL_5")){
  # organ <- "GL_5"
  gene_list_PCG <- unique(DE_Acom_GL_PCGs$gene[grep(organ, DE_Acom_GL_PCGs$cell_type)])
  PCG_list[[organ]] <- gene_list_PCG
  
  gene_list <- sort(gene_list_PCG, decreasing = TRUE)
  gene_id_list <- bitr(gene_list,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Mm.eg.db)
  gene_id_list <- dplyr::distinct(gene_id_list,SYMBOL,.keep_all=TRUE) # remove duplication
  gene_id_list_sort <- sort(gene_id_list$ENTREZID, decreasing = TRUE)
  
  go_enrich_BP <- enrichGO(gene = gene_id_list_sort,
                           OrgDb = org.Mm.eg.db,
                           keyType = 'ENTREZID',
                           readable = T,
                           ont = "ALL",
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

cell_type <- c("GL_1|GL_2","GL_3|GL_4","GL_5")
my_color <- ann_colors$cluster[c(1,4,5)]
names(my_color) <- cell_type

for(i in 1:3){
  i <- 3
  GO_barplot_all <- c()
  GO_barplot_list <- list()
  for (organ in cell_type[i]){
    # organ <- "Radial_glia"
    go_enrich_BP <- read.csv(paste0(organ, "_DE_ABcom_allDEG_enrichGO_GOBP.txt"), sep="\t", header=T)
    go_enrich_BP$Description <- paste(go_enrich_BP$ONTOLOGY,go_enrich_BP$Description, sep = ":")
    GO_barplot <- data.frame(go_enrich_BP[1:6, ][, c("Description", "pvalue", "Count")], organ)
    
    GO_barplot$pvalue <- -log10(GO_barplot$pvalue)
    colnames(GO_barplot) <- c("GO_term", "pvalue", "Count", "organ")
    GO_barplot$GO_term <- factor(GO_barplot$GO_term, levels = rev(GO_barplot$GO_term))
    GO_barplot_all <- rbind(GO_barplot_all, GO_barplot)
  }
  
  GO_barplot_all$organ <- factor(GO_barplot_all$organ, levels = cell_type[i])
  p <- ggplot(GO_barplot_all, aes(x=GO_term,y=pvalue, fill=organ)) + geom_bar(stat="identity", width=0.8) + coord_flip() +  xlab("GO terms") + ylab("-log10(pvalue)") +
    facet_wrap(~organ, scales = "free_y") +
    theme_bw() + 
    theme(legend.position="none") + 
    theme(panel.grid = element_blank()) + 
    scale_fill_manual(values=my_color)
  
  pdf(paste0("DE_Acom_allPCG_enrichGO_GOBP_barplot_", cell_type[i], ".pdf"), width = 10, height = 3)
  print(p)
  dev.off()
}

PCG_list <- list()
for (organ in c("GL_1|GL_2","GL_3|GL_4","GL_5")){
  
  organ <- "GL_1|GL_2"
  
  gene_list_PCG <- unique(DE_Acom_GL_PCGs$gene[grep(organ, DE_Acom_GL_PCGs$cell_type)])
  PCG_list[[organ]] <- gene_list_PCG
  
  gene_list <- sort(gene_list_PCG, decreasing = TRUE)
  gene_id_list <- bitr(gene_list,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Mm.eg.db)
  gene_id_list <- dplyr::distinct(gene_id_list,SYMBOL,.keep_all=TRUE) # remove duplication
  gene_id_list_sort <- sort(gene_id_list$ENTREZID, decreasing = TRUE)
  
  go_enrich_BP <- enrichGO(gene = gene_id_list_sort,
                           OrgDb = org.Mm.eg.db,
                           keyType = 'ENTREZID',
                           readable = T,
                           ont = "ALL",
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH")
  
  pdf(paste0("DE_Acom_allPCG_enrichGO_GOBP_cnetplot_CC_", organ, ".pdf"), width = 5, height = 5)
  cnetplot(go_enrich_BP, showCategory = go_enrich_BP@result$Description[go_enrich_BP@result$ONTOLOGY=="CC"], foldChange=gene_list)
  dev.off()
}


########################################## revision (Rows and columns can be visualized on the UMAPs)
# cluster highlight
source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_hic/Cere-hic-13"
sample <- "Cere-hic-13"
PC <- 36
res <- 1
higashi_out <- "8_higashi_out_higashi_500kb"
work_path_seurat <- "10_seurat_bigSpotShape_higashi_500kb"
embed_type <- "scAB_scale"
setwd(paste0(work_path,  "/", work_path_seurat, "/", embed_type))
load(paste0("Out/", sample, "_Seurat_object_anno.RData"))

library(future)
options(future.globals.maxSize = 2 * 1024^3)  # 设置为 2GB
Seurat_object <- SCTransform(Seurat_object, assay = "Spatial", verbose = T, variable.features.n = 2000)
Seurat_object <- RunPCA(Seurat_object, assay = "SCT", verbose = FALSE)
ElbowPlot(Seurat_object, ndims = 50)
DimPlot(Seurat_object, reduction = "umap", group.by = "seurat_clusters", label = FALSE)

Seurat_object$row <- Seurat_object@images$Cere.hic.13@coordinates$row
Seurat_object$col <- Seurat_object@images$Cere.hic.13@coordinates$col

for(i in 1:96){
  # i <- 1
  Seurat_object$row_number <- NA
  Seurat_object$row_number[which(Seurat_object$row==i)] <- "highlight"
  Seurat_object$row_number[which(Seurat_object$row!=i)] <- "other"
  Seurat_object[[paste0("row_", i)]] <- Seurat_object$row_number
}

for(i in 1:96){
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

for(i in 1:96){
  
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
pdf(paste0("/media/maoni/data/CZP/Figures/revision/UMAP_col_row/", sample, "_row_umap.pdf"), height = 12*3, width = 8*3)
ggarrange(plotlist = DimPlot2_list, ncol = 8, nrow = 12)
dev.off()


DimPlot2_list <- list()

for(i in 1:96){
  
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
pdf(paste0("/media/maoni/data/CZP/Figures/revision/UMAP_col_row/", sample, "_col_umap.pdf"), height = 12*3, width = 8*3)
ggarrange(plotlist = DimPlot2_list, ncol = 8, nrow = 12)
dev.off()

pdf(paste0("/media/maoni/data/CZP/Figures/revision/UMAP_col_row/", sample, "_UMAP_col_row_ARI_plot.pdf"), width = 10, height = 5)
par(mfrow=c(1,2))
ARI_col <- read.csv(paste0("/media/maoni/data/CZP/Figures/revision/UMAP_col_row/", sample, "_ARI_col_vs_seurat_clusters.csv"), header = TRUE)
plot(ARI_col$ARI, type = "p", pch=16, col = "#E69F00", ylim = c(-0.1, 0.1), xlab = "col_number", ylab = "ARI_value", main = "ARI_col")
ARI_row <- read.csv(paste0("/media/maoni/data/CZP/Figures/revision/UMAP_col_row/", sample, "_ARI_row_vs_seurat_clusters.csv"), header = TRUE)
plot(ARI_row$ARI, type = "p", pch=16, col = "#56B4E9", ylim = c(-0.1, 0.1), xlab = "row_number", ylab = "ARI_value", main = "ARI_row")
dev.off()

############################### revision homer pearson plot ##############################################3
library(pheatmap)
library(gridExtra)
# 自定义调色板
my_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# 读取并生成 pheatmap grobs

# 定义缩放函数到 [-1, 1]

scale_to_minus1_1 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) {
    return(matrix(0, nrow = nrow(x), ncol = ncol(x)))  # 所有值相等时，返回0矩阵
  } else {
    return(2 * (x - rng[1]) / (rng[2] - rng[1]) - 1)
  }
}


region <- "chr16:50000000-90000000"
clusters <- c("GL_1", "GL_2", "GL_3", "GL_4", "GL_5", "PL_1", "PL_2")
heatmaps <- list()
for (i in 1:length(clusters)) {
  i <- 1
  mat <- read.table(paste0("/media/maoni/data/CZP/spatial_hic/hicup/Cere/Cereall/R1.", clusters[i], "/homer/Pearson/cis.Brainall.chrall.R1.", clusters[i], ".homer.pearson.res100k.win400k.", region, ".txt"), sep = "\t", header = TRUE)
  mat <- mat[, -1]
  rownames(mat) <- mat[,1]
  mat <- mat[, -1]
  mat <- as.matrix(mat)
  
  mat_scaled <- scale_to_minus1_1(mat)
  write.table(mat_scaled, paste0("/media/maoni/data/CZP/spatial_hic/hicup/Cere/Cereall/R1.", clusters[i], "/homer/Pearson/cis.Brainall.chrall.R1.", clusters[i], ".homer.pearson.res100k.win400k.", region, "_scale_to_minus1_1.txt"), sep = "\t", col.names = T, row.names = T, quote = F)
  
  # 生成热图 grob（绘图对象），不立即绘制
  heatmaps[[i]] <- pheatmap(mat_scaled,
                            color = my_colors,
                            cluster_rows = FALSE,
                            cluster_cols = FALSE,
                            show_rownames = FALSE,
                            show_colnames = FALSE,
                            silent = TRUE)
}

# 输出到 PDF，宽度设置为足够容纳7个图
pdf("combined_heatmaps.pdf", width = 14, height = 4)

# 1行7列排布
grid.arrange(grobs = heatmaps, ncol = 7)

dev.off()


