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
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(png)
library(tidyr)
library(cowplot)
library(viridis)

source("/media/maoni/data/R_functions_ST_E13.R")
work_path <- "/media/maoni/data/CZP/spatial_transcriptome/E13-23-RNA"
sample <- "E13-23-RNA"

# 0 file preparing
set.seed(123)
setwd(work_path)
dir.create("6_seurat_bigSpotShape_new")
setwd(paste0(work_path, "/6_seurat_bigSpotShape_new"))
dir.create("Img")
dir.create("Matrix")
dir.create("Out")
dir.create("Plot")
file.copy(from=paste0(work_path, "/3_umi_tools/", sample, ".debarcoded_passed_reads_stat.csv"), 
          to=paste0(work_path, "/6_seurat_bigSpotShape_new/Matrix/"), overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
file.copy(from=paste0(work_path, "/4_zUMIs/", sample, ".debarcoded_passed_reads_stat_cropped.csv"), 
          to=paste0(work_path, "/6_seurat_bigSpotShape_new/Matrix/"), overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
file.copy(from=paste0(work_path, "/4_zUMIs/", sample, ".dgecounts.rds"), 
          to=paste0(work_path, "/6_seurat_bigSpotShape_new/Matrix/"), overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)

# 0 Load data from st pipe
setwd(paste0(work_path, "/6_seurat_bigSpotShape_new"))
# run this perl cmd
# perl -lane 'if(/gene_id\s.(.+?).\;.+?gene_name\s.(.+?).\;.+?gene_biotype\s.(.+?).\;/){print "$1\t$2\t$3"}' Mus_musculus.GRCm38.102.chr.gtf | sort -u > gGRCm38.102.gene.ids.names
gene.names <- read.table("/media/maoni/data/Reference/mouse/GRCm38_mm10/ensembl/gGRCm38.102.gene.ids.names", sep = "\t")

# remove unnecessary pseudogene, rRNA, snoRNA, tRNA, ribozyme
removed.genes1 <- gene.names$V2[grep("pseudogene|tRNA|rRNA|snoRNA|ribozyme|snRNA|scaRNA|misc_RNA|scRNA|sRNA", gene.names$V3)]
removed.genes2 <- gene.names$V2[grep("mt-|Rps|Rpl|Hba|Hbb", gene.names$V2)]  # 1829
removed.genes <- unique(c(removed.genes1, removed.genes2))  # 38300

allCounts <- readRDS(paste0("Matrix/", sample, ".dgecounts.rds"))
count_matrix <- allCounts$umicount$inex$all
count_matrix <- as.matrix(count_matrix)

raw_stat = read.table(paste0("Matrix/", sample, ".debarcoded_passed_reads_stat.csv"),sep =",", header = TRUE, dec =".", stringsAsFactors = F)
raw_stat_cell = paste(51-raw_stat$iA, raw_stat$iB, sep="x")
raw_stat_barcode = paste0(raw_stat$bc_B, raw_stat$bc_A)
colnames(count_matrix) <- raw_stat_cell[match(colnames(count_matrix),raw_stat_barcode)]
rownames(count_matrix) <- gene.names$V2[match(rownames(count_matrix), gene.names[,1])]  # 31867  2401

ID <- as.character(rownames(count_matrix))
d <- duplicated(ID)
ID <- factor(ID, levels = unique(ID))
count_matrix <- rowsum(as.matrix(count_matrix), ID, reorder = FALSE, na.rm = TRUE)
removed.gene <- intersect(rownames(count_matrix), removed.genes)
count_matrix <- count_matrix[-match(removed.gene, rownames(count_matrix)),]  #27336  2401
dim(count_matrix)

# 1 prepare cropped matrix
stat_cropped = read.table(paste0("Matrix/", sample, ".debarcoded_passed_reads_stat_cropped.csv"),sep =",", header = TRUE, dec =".", stringsAsFactors = F)
iA_iB <- paste(51-stat_cropped$iA, stat_cropped$iB, sep="x")
count_matrix_cropped <- count_matrix[,match(iA_iB,colnames(count_matrix))]
dim(count_matrix_cropped)  # 30027  2045

# 2 Seurat object
# See Load10X_Spatial; Read10X_Image
# See satijalab/seurat/issue/3539 4993
assay = "Spatial"
slice = sample
image.dir = "/media/maoni/data/CZP/spatial_transcriptome/seurat_yuhao/core/50*50_blank_img" # "./Img"
image.nam = paste0(sample,"_fix_1080.png") # "grey_pixel_1080p.png"
coord.nam = "combine_barcode.round2round1_index1_index2.Seurat.txt"
Seurat_object_cropped = CreateSeuratObject(counts = count_matrix_cropped, project = sample,assay = assay, min.cells = 3, min.features = 0) # min.cells min.features
image <- readPNG(source = file.path(image.dir, image.nam))[,,1:3]
scale.factors <- c("tissue_hires_scalef"=1, "fiducial_diameter_fullres"=1, "tissue_lowres_scalef"=1)
tissue.positions <- read.table(file = file.path(image.dir,coord.nam), col.names = c("barcodes", "tissue", "row", "col", "imagerow", "imagecol"), header = FALSE, as.is = TRUE, row.names = 1)
spot.radius <- 0.015 # estiamte:(0.13)*50/410/2
image <- new(Class = "VisiumV1", image = image, scale.factors = scalefactors(spot = scale.factors[1], fiducial = scale.factors[2], hires = scale.factors[1], lowres = scale.factors[3]), coordinates = tissue.positions, spot.radius = spot.radius)
image <- image[Cells(Seurat_object_cropped)]
DefaultAssay(object = image) <- assay
Seurat_object_cropped[[slice]] <- image

# 3 Data quality
Seurat_object_cropped[["percent.mMt"]] <- PercentageFeatureSet(Seurat_object_cropped, pattern = "^mt-")
Seurat_object_cropped[["percent.mRp"]] <- PercentageFeatureSet(Seurat_object_cropped, pattern = "^Rp[sl]")
Seurat_object_cropped[["percent.mHb"]] <- PercentageFeatureSet(Seurat_object_cropped, pattern = "^Hba|^Hbb")

plot0 <- VlnPlot(Seurat_object_cropped, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mMt", "percent.mRp", "percent.mHb"), ncol = 5, pt.size = 0)
plot0

plot1 <- FeatureScatter(Seurat_object_cropped, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plot1

plot3 <- SpatialFeaturePlot(Seurat_object_cropped, features = "nCount_Spatial", pt.size.factor = 1.6, stroke = NA) + theme(legend.position = "top")
plot3$layers[[1]]$aes_params=c(plot3$layers[[1]]$aes_params, shape=22)
plot4 <- SpatialFeaturePlot(Seurat_object_cropped, features = "nFeature_Spatial", pt.size.factor = 1.6, stroke = NA) + theme(legend.position = "top")
plot4$layers[[1]]$aes_params=c(plot4$layers[[1]]$aes_params, shape=22)
plot3 + plot4

pdf("Plot/1_QC_check_vlnplot.pdf", width = 10, height = 4)
print(plot0)
dev.off()
pdf("Plot/1_QC_check_spatial.pdf", width = 10, height = 4)
print(plot1 + plot3 + plot4)
dev.off()

# 4 Run PCA
Seurat_object_cropped <- SCTransform(Seurat_object_cropped, assay = "Spatial", verbose = FALSE)
Seurat_object_cropped <- RunPCA(Seurat_object_cropped, assay = "SCT", verbose = FALSE)
ElbowPlot(Seurat_object_cropped, ndims=50)

# 5 Cluster
dim_p <- list()
spatial_p <- list()
pcNO <- 30
for(res in c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5, 2)){
  Seurat_object_cropped <- FindNeighbors(Seurat_object_cropped, reduction = "pca", dims = 1:pcNO)
  Seurat_object_cropped <- FindClusters(Seurat_object_cropped, verbose = FALSE, resolution = res)
  Seurat_object_cropped <- RunUMAP(Seurat_object_cropped, reduction = "pca", dims = 1:pcNO)
  dim_p[[paste0("res",res)]] <- DimPlot(Seurat_object_cropped, reduction = "umap", label = TRUE) + ggtitle(paste0("res",res)) # the UMAP plot
  spatial_p[[paste0("res",res)]] <- SpatialDimPlot(Seurat_object_cropped, label = FALSE, label.size = 3, pt.size.factor = 1.6, stroke = NA) + ggtitle(paste0("res",res)) # the spatial plot
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
res <- 0.5
for(pcNO in seq(10,50,by=5)){
  Seurat_object_cropped <- FindNeighbors(Seurat_object_cropped, reduction = "pca", dims = 1:pcNO)
  Seurat_object_cropped <- FindClusters(Seurat_object_cropped, verbose = FALSE, resolution = res)
  Seurat_object_cropped <- RunUMAP(Seurat_object_cropped, reduction = "pca", dims = 1:pcNO)
  dim_p[[paste0("PC",pcNO)]] <- DimPlot(Seurat_object_cropped, reduction = "umap", label = TRUE) + ggtitle(paste0("PC",pcNO)) # the UMAP plot
  spatial_p[[paste0("PC",pcNO)]] <- SpatialDimPlot(Seurat_object_cropped, label = FALSE, label.size = 3, pt.size.factor = 1.6, stroke = NA) + ggtitle(paste0("PC",pcNO)) # the spatial plot
  spatial_p[[paste0("PC",pcNO)]]$layers[[1]]$aes_params=c(spatial_p[[paste0("PC",pcNO)]]$layers[[1]]$aes_params, shape=22)
}

pdf(file = paste("Plot/2_UMAP_clusters_res",res,".pdf",sep =""), width = 28, height = ceiling(length(dim_p)/4)*7)
ggarrange(plotlist = dim_p, ncol = 4, nrow = ceiling(length(dim_p)/4))
dev.off()

pdf(file = paste("Plot/2_Spatial_clusters_res",res,".pdf",sep =""), width = 28, height = ceiling(length(spatial_p)/4)*7)
ggarrange(plotlist = spatial_p, ncol = 4, nrow = ceiling(length(spatial_p)/4))
dev.off()

# Select PC and res
PC=30
res=0.8
Seurat_object_cropped <- FindNeighbors(Seurat_object_cropped, reduction = "pca", dims = 1:PC)
Seurat_object_cropped <- FindClusters(Seurat_object_cropped, verbose = FALSE, resolution = res)
Seurat_object_cropped <- RunUMAP(Seurat_object_cropped, reduction = "pca", dims = 1:PC)
umap_Spa <- list()
umap_Spa[["p1"]] <- DimPlot(Seurat_object_cropped, reduction = "umap", label = TRUE)
umap_Spa[["p2"]] <- SpatialDimPlot(Seurat_object_cropped, label = FALSE, label.size = 3, pt.size.factor = 1.6, stroke = NA)
umap_Spa[["p2"]]$layers[[1]]$aes_params=c(umap_Spa[["p2"]]$layers[[1]]$aes_params, shape=22)
umap_Spa[["p3"]] <- DimPlot(Seurat_object_cropped, reduction = "umap", label = TRUE)
umap_Spa[["p4"]] <- SpatialDimPlot(Seurat_object_cropped, label = FALSE, label.size = 3, pt.size.factor = 1.6, stroke = NA)
umap_Spa[["p4"]]$layers[[1]]$aes_params=c(umap_Spa[["p4"]]$layers[[1]]$aes_params, shape=22)

pdf(file = paste0("Plot/2_UMAP_clusters_PC",PC,"_res",res,".pdf"), width = 28, height = ceiling(length(umap_Spa)/4)*7)
p <- ggarrange(plotlist = umap_Spa, ncol = 4, nrow = ceiling(length(umap_Spa)/4))
print(p)
dev.off()

spatial_p_cluster <- list()
for(i in sort(unique(Idents(Seurat_object_cropped)))){
  spatial_p_cluster[[paste0("idents_",i)]] <- SpatialDimPlot(Seurat_object_cropped, cells.highlight = CellsByIdentities(object = Seurat_object_cropped, idents = i), cols.highlight = c("#DE2D26", "grey90"), facet.highlight = TRUE, ncol = 1, pt.size.factor = 1.6, stroke = NA, label.size = 12) + ggtitle(paste0("idents",i)) # the spatial plot
  spatial_p_cluster[[paste0("idents_",i)]]$layers[[1]]$aes_params=c(spatial_p_cluster[[paste0("idents_",i)]]$layers[[1]]$aes_params, shape=22)
}

pdf(file = paste0("Plot/2_Spatial_clusters_PC",PC,"_res",res,".pdf"), width = 28, height = ceiling(length(spatial_p_cluster)/4)*7)
p <- ggarrange(plotlist = spatial_p_cluster, ncol = 4, nrow = ceiling(length(spatial_p_cluster)/4))
print(p)
dev.off()

pdf_combine(input = c(paste0("Plot/2_UMAP_clusters_PC",PC,"_res",res,".pdf"), paste0("Plot/2_Spatial_clusters_PC",PC,"_res",res,".pdf")), output = paste0("Plot/2_UMAP_Spatial_clusters_PC",PC,"_res",res,".pdf"))


# 6 DE
# de_markers <- FindMarkers(Seurat_object_cropped, ident.1 = 5, ident.2 = 6)
de_markers <- FindAllMarkers(Seurat_object_cropped, only.pos = T, min.pct = 0.1, logfc.threshold = 0.25) %>% arrange(cluster, desc(avg_log2FC))  #8434
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)  %>% arrange(cluster, desc(avg_log2FC))

png(paste0("Plot/3_findmarker_UMAP_data_PC",PC,"_res",res,".png"), height = ceiling(length(top5$cluster)/5)*480, width = 5*480)
p <- FeaturePlot(Seurat_object_cropped, features = top5$gene, ncol = 5,pt.size = 0.1,order = T)
print(p)
dev.off()

Seurat_object_cropped <- ScaleData(Seurat_object_cropped, features = top5$gene)
pdf(paste0("Plot/3_findmarker_heatmap_PC",PC,"_res",res,".pdf"),height = 14, width = 28)
p <- DoHeatmap(Seurat_object_cropped, features = top5$gene,size = 5.5)+
        theme(plot.title = element_text(hjust = 0.5, size = 5, face = "bold"),
        axis.text=element_text(size=5,face="bold"),
        axis.title=element_text(size=5,face="bold"),
        legend.text=element_text(size=10),
        legend.title = element_blank())
print(p)
dev.off()

pdf(paste0("Plot/3_findmarker_vlnPlot_data_PC",PC,"_res",res,".pdf"), height = 48, width = 14)
p <- VlnPlot(Seurat_object_cropped, split.by = "seurat_clusters", features = top5$gene, flip = T, stack = T)
print(p)
dev.off()

write.csv(top5,row.names = F,file = paste0("Out/1_de_genes_top5_PC",PC,"_res",res,".csv"))
write.csv(de_markers,row.names = F,file = paste0("Out/1_de_genes_PC",PC,"_res",res,".csv"))
write.csv(Seurat_object_cropped@meta.data,row.names = T,file = paste0("Out/1_cell_metadata_PC",PC,"_res",res,".csv"))
save(top5, de_markers, Seurat_object_cropped, file = paste0("Out/", sample, "_cropped_PC",PC,"_res",res,".RData"))

# 7 magic
source("/media/maoni/data/R_functions_ST_E13.R")
Seurat_object_cropped_magic <- magic.Seurat(Seurat_object_cropped)
DefaultAssay(object = Seurat_object_cropped_magic) <- "MAGIC_SCT"

png(paste0("Plot/3_findmarker_Spatial_data_PC", PC, "_res", res,".png"), height = ceiling(length(top5$cluster)/5)*480, width = 5*480)
plot <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = top5$gene, alpha = c(0.1, 1), ncol = 5, pt.size.factor = 1.6) # + scale_fill_gradientn(colours = viridis::inferno(100))
plot$layers[[1]]$aes_params=c(plot$layers[[1]]$aes_params, shape=22)
plot
dev.off()

saveRDS(Seurat_object_cropped, file = paste0("Out/", sample, "_cropped", ".rds"))
saveRDS(Seurat_object_cropped_magic, file = paste0("Out/", sample, "_cropped_magic", ".rds"))
save(Seurat_object_cropped, Seurat_object_cropped_magic, file = paste0("Out/", sample, "_cropped_and_", sample, "_cropped_magic.RData"))

# 8 annotation
setwd("/media/maoni/data/CZP/spatial_transcriptome/E13-23-RNA/6_seurat_bigSpotShape_new")
PC=30
res=0.8
load(paste0("Out/", sample, "_cropped_PC",PC,"_res",res,".RData"))

anno_marker <- rep(NA,nrow(de_markers))
anno_marker[de_markers$cluster %in% c(0)] <- "Excitatory/Inhibitory neurons"
anno_marker[de_markers$cluster %in% c(1)] <- "Connective tissue progenitors"
anno_marker[de_markers$cluster %in% c(2)] <- "Inhibitory interneurons"
anno_marker[de_markers$cluster %in% c(3)] <- "Ependymal cell"
anno_marker[de_markers$cluster %in% c(4)] <- "Hepatocytes & Primitive erythroid lineage"
anno_marker[de_markers$cluster %in% c(5,15)] <- "Chondrocytes & osteoblasts"
anno_marker[de_markers$cluster %in% c(6,14)] <- "Epithelial cells"
anno_marker[de_markers$cluster %in% c(7,10)] <- "Myocytes"
anno_marker[de_markers$cluster %in% c(8)] <- "Cardiac muscle lineages"
anno_marker[de_markers$cluster %in% c(9)] <- "Radial glia"
anno_marker[de_markers$cluster %in% c(11)] <- "Osteoblasts"
anno_marker[de_markers$cluster %in% c(12)] <- "Other"
anno_marker[de_markers$cluster %in% c(13)] <- "Cholinergic neurons"
anno_marker[de_markers$cluster %in% c(16)] <- "Inhibitory neuron progenitors"
anno_marker[de_markers$cluster %in% c(17)] <- "Premature oligodendrocyte"
anno_marker[de_markers$cluster %in% c(18)] <- "Endothelial cells"

de_markers$cell_type <- factor(anno_marker,levels = c("Inhibitory neuron progenitors", "Inhibitory interneurons", "Excitatory/Inhibitory neurons", "Cholinergic neurons", "Radial glia", "Premature oligodendrocyte",
                                                      "Cardiac muscle lineages", "Hepatocytes & Primitive erythroid lineage", 
                                                      "Connective tissue progenitors", "Ependymal cell", "Chondrocytes & osteoblasts", "Epithelial cells", "Myocytes", "Osteoblasts", "Endothelial cells", "Other")) # mysemchyme

anno <- rep(NA,ncol(Seurat_object_cropped))
anno[Seurat_object_cropped$seurat_clusters %in% c(0)] <- "Excitatory/Inhibitory neurons"
anno[Seurat_object_cropped$seurat_clusters %in% c(1)] <- "Connective tissue progenitors"
anno[Seurat_object_cropped$seurat_clusters %in% c(2)] <- "Inhibitory interneurons"
anno[Seurat_object_cropped$seurat_clusters %in% c(3)] <- "Ependymal cell"
anno[Seurat_object_cropped$seurat_clusters %in% c(4)] <- "Hepatocytes & Primitive erythroid lineage"
anno[Seurat_object_cropped$seurat_clusters %in% c(5,15)] <- "Chondrocytes & osteoblasts"
anno[Seurat_object_cropped$seurat_clusters %in% c(6,14)] <- "Epithelial cells"
anno[Seurat_object_cropped$seurat_clusters %in% c(7,10)] <- "Myocytes"
anno[Seurat_object_cropped$seurat_clusters %in% c(8)] <- "Cardiac muscle lineages"
anno[Seurat_object_cropped$seurat_clusters %in% c(9)] <- "Radial glia"
anno[Seurat_object_cropped$seurat_clusters %in% c(11)] <- "Osteoblasts"
anno[Seurat_object_cropped$seurat_clusters %in% c(12)] <- "Other"
anno[Seurat_object_cropped$seurat_clusters %in% c(13)] <- "Cholinergic neurons"
anno[Seurat_object_cropped$seurat_clusters %in% c(16)] <- "Inhibitory neuron progenitors"
anno[Seurat_object_cropped$seurat_clusters %in% c(17)] <- "Premature oligodendrocyte"
anno[Seurat_object_cropped$seurat_clusters %in% c(18)] <- "Endothelial cells"

Seurat_object_cropped$cell_type <- factor(anno,levels = c("Inhibitory neuron progenitors", "Inhibitory interneurons", "Excitatory/Inhibitory neurons", "Cholinergic neurons", "Radial glia", "Premature oligodendrocyte",
                                                          "Cardiac muscle lineages", "Hepatocytes & Primitive erythroid lineage", 
                                                          "Connective tissue progenitors", "Ependymal cell", "Chondrocytes & osteoblasts", "Epithelial cells", "Myocytes", "Osteoblasts", "Endothelial cells", "Other")) # mysemchyme

aa <- c("Inhibitory neuron progenitors", "Inhibitory interneurons", "Excitatory/Inhibitory neurons", "Cholinergic neurons", "Radial glia", "Premature oligodendrocyte",
        "Cardiac muscle lineages", "Hepatocytes & Primitive erythroid lineage", 
        "Connective tissue progenitors", "Ependymal cell", "Chondrocytes & osteoblasts", "Epithelial cells", "Myocytes", "Osteoblasts", "Endothelial cells", "Other")

aa <- paste0("R", seq_along(aa), "_", aa)
print(aa)
bb <- data.frame(aa)
write.csv(bb,row.names = F, file = paste0("Out/Cell_type_anno_table.csv"))

# 16 cell types
source("/media/maoni/data/Zhanglin_Junjie/ColorPalettes.R")
my_color_palette <- ArchRPalettes$bear
names(my_color_palette) <- c("Inhibitory neuron progenitors", "Inhibitory interneurons", "Excitatory/Inhibitory neurons", "Cholinergic neurons", "Radial glia", "Premature oligodendrocyte",
                             "Cardiac muscle lineages", "Hepatocytes & Primitive erythroid lineage", 
                             "Connective tissue progenitors", "Ependymal cell", "Chondrocytes & osteoblasts", "Epithelial cells", "Myocytes", "Osteoblasts", "Endothelial cells", "Other")

umap_Spa <- list()
umap_Spa[["p1"]] <- DimPlot(Seurat_object_cropped, reduction = "umap", group.by = "cell_type", label = FALSE, pt.size = 0.5, cols = my_color_palette) + theme(aspect.ratio = 1)
umap_Spa[["p2"]] <- SpatialDimPlot(Seurat_object_cropped, group.by = "cell_type", label = FALSE, label.size = 3, pt.size.factor = 1.6, stroke = NA, cols = my_color_palette) + theme(aspect.ratio = 1)
umap_Spa[["p2"]]$layers[[1]]$aes_params=c(umap_Spa[["p2"]]$layers[[1]]$aes_params, shape=22)
umap_Spa[["p3"]] <- DimPlot(Seurat_object_cropped, reduction = "umap", group.by = "cell_type", label = FALSE, pt.size = 0.5, cols = my_color_palette)+ theme(aspect.ratio = 1)
umap_Spa[["p4"]] <- SpatialDimPlot(Seurat_object_cropped, group.by = "cell_type", label = FALSE, label.size = 3, pt.size.factor = 1.6, stroke = NA, cols = my_color_palette)+ theme(aspect.ratio = 1)
umap_Spa[["p4"]]$layers[[1]]$aes_params=c(umap_Spa[["p4"]]$layers[[1]]$aes_params, shape=22)

pdf(file = paste0("Plot/2_UMAP_clusters_PC",PC,"_res",res,"_anno.pdf"), width = 28, height = ceiling(length(umap_Spa)/4)*7)
p <- ggarrange(plotlist = umap_Spa, ncol = 4, nrow = ceiling(length(umap_Spa)/4))
print(p)
dev.off()

spatial_p_cluster <- list()
for(i in names(table(Seurat_object_cropped$cell_type))){
  Seurat_object_cropped_subset <- Seurat_object_cropped[, Seurat_object_cropped$cell_type==i]
  spatial_p_cluster[[i]] <- SpatialDimPlot(Seurat_object_cropped_subset, group.by = "cell_type", crop = FALSE, label = FALSE, label.size = 3, pt.size.factor = 1.6, stroke = NA,  label.box = FALSE, cols = my_color_palette) + theme(aspect.ratio = 1) + ggtitle(i) + theme(legend.position="none") # the spatial plot
  spatial_p_cluster[[i]]$layers[[1]]$aes_params=c(spatial_p_cluster[[i]]$layers[[1]]$aes_params, shape=22)
}

pdf(file = paste0("Plot/2_Spatial_clusters_PC",PC,"_res",res,"_anno.pdf"), width = 28, height = ceiling(length(spatial_p_cluster)/4)*7)
p <- ggarrange(plotlist = spatial_p_cluster, ncol = 4, nrow = ceiling(length(spatial_p_cluster)/4))
print(p)
dev.off()

pdf_combine(input = c(paste0("Plot/2_UMAP_clusters_PC",PC,"_res",res,"_anno.pdf"), paste0("Plot/2_Spatial_clusters_PC",PC,"_res",res,"_anno.pdf")), output = paste0("Plot/2_UMAP_Spatial_clusters_PC",PC,"_res",res,"_anno.pdf"))

top3_anno <- de_markers %>% group_by(cell_type) %>% top_n(n = 3, wt = avg_log2FC)  %>% arrange(cell_type, desc(avg_log2FC))
top5_anno <- de_markers %>% group_by(cell_type) %>% top_n(n = 5, wt = avg_log2FC)  %>% arrange(cell_type, desc(avg_log2FC))
top10_anno <- de_markers %>% group_by(cell_type) %>% top_n(n = 10, wt = avg_log2FC)  %>% arrange(cell_type, desc(avg_log2FC))

de_markers <- de_markers %>% group_by(cell_type) %>% arrange(cell_type, desc(avg_log2FC))
de_markers_HBFB <- de_markers[grep("Inhibitory neuron progenitors|Inhibitory interneurons|Excitatory/Inhibitory neurons|Cholinergic neurons|Radial glia|Premature oligodendrocyte", de_markers$cell_type), ]
de_markers_HBHT <- de_markers[grep("Inhibitory neuron progenitors|Excitatory/Inhibitory neurons|Cardiac muscle lineages", de_markers$cell_type), ]


Idents(Seurat_object_cropped) <- Seurat_object_cropped$cell_type
de_markers_NeuMes <- FindMarkers(Seurat_object_cropped, ident.1 = c("Inhibitory neuron progenitors", "Inhibitory interneurons", "Excitatory/Inhibitory neurons", "Cholinergic neurons", "Radial glia", "Premature oligodendrocyte"), ident.2 = c("Connective tissue progenitors", "Chondrocytes & osteoblasts", "Myocytes"))
de_markers_NeuMes_2 <- de_markers[sort(na.omit(match(rownames(de_markers_NeuMes), de_markers$gene))), ]
write.csv(de_markers_NeuMes_2,row.names = T,file = paste0("Out/1_de_genes_PC",PC,"_res",res,"_NeuMes.csv"))


All_gene_Anno_PCG <- read.table("/media/maoni/data/Reference/mouse/GRCm38_mm10/ensembl/All_gene_Anno_PCG.bed", sep = "\t", header = FALSE)

de_markers_NeuMes_bed <- All_gene_Anno_PCG[na.omit(match(unique(de_markers_NeuMes_2$gene), All_gene_Anno_PCG$V4)), ]
de_markers_NeuMes_bed_2 <- data.frame(de_markers_NeuMes_bed, de_markers_NeuMes_2$cell_type[na.omit(match(de_markers_NeuMes_bed$V4, de_markers_NeuMes_2$gene))])

de_markers_HBFB_bed <- All_gene_Anno_PCG[na.omit(match(unique(de_markers_HBFB$gene), All_gene_Anno_PCG$V4)), ]
de_markers_HBFB_bed_2 <- data.frame(de_markers_HBFB_bed, de_markers_HBFB$cell_type[na.omit(match(de_markers_HBFB_bed$V4, de_markers_HBFB$gene))])

de_markers_NeuMes_bed <- All_gene_Anno_PCG[na.omit(match(unique(de_markers_HBFB$gene), All_gene_Anno_PCG$V4)), ]
de_markers_NeuMes_bed_2 <- data.frame(de_markers_HBFB_bed, de_markers_HBFB$cell_type[na.omit(match(de_markers_HBFB_bed$V4, de_markers_HBFB$gene))])

de_markers_HBHT_bed <- All_gene_Anno_PCG[na.omit(match(unique(de_markers_HBHT$gene), All_gene_Anno_PCG$V4)), ]
de_markers_HBHT_bed_2 <- data.frame(de_markers_HBHT_bed, de_markers_HBHT$cell_type[na.omit(match(de_markers_HBHT_bed$V4, de_markers_HBHT$gene))])

write.table(de_markers_HBFB_bed_2, "Out/de_markers_HBFB.bed",sep="\t",quote = F,row.names = F,col.names = F)
write.table(de_markers_HBHT_bed_2, "Out/de_markers_HBHT.bed",sep="\t",quote = F,row.names = F,col.names = F)
write.table(de_markers_HBFB_bed_2, "Out/de_markers_FBHB.bed",sep="\t",quote = F,row.names = F,col.names = F)
write.table(de_markers_HBHT_bed_2, "Out/de_markers_HTHB.bed",sep="\t",quote = F,row.names = F,col.names = F)

write.table(de_markers_NeuMes_bed_2, "Out/de_markers_NeuMes.bed",sep="\t",quote = F,row.names = F,col.names = F)
write.table(de_markers_NeuMes_bed_2, "Out/de_markers_MesNeu.bed",sep="\t",quote = F,row.names = F,col.names = F)

de_markers_NeuMes_2 <- read.csv(paste0("Out/1_de_genes_PC",PC,"_res",res,"_NeuMes.csv"), header = T, row.names = 1)
de_markers_Neu_2 <- de_markers_NeuMes_2[grep("Inhibitory neuron progenitors|Inhibitory interneurons|Excitatory/Inhibitory neurons|Cholinergic neurons|Radial glia|Premature oligodendrocyte", de_markers_NeuMes_2$cell_type), ]
de_markers_Neu_3 <- de_markers_Neu_2[de_markers_Neu_2$avg_log2FC > 1, ]
de_markers_Neu_bed <- All_gene_Anno_PCG[na.omit(match(unique(de_markers_Neu_3$gene), All_gene_Anno_PCG$V4)), ]
de_markers_Neu_bed_2 <- data.frame(de_markers_Neu_bed, de_markers_Neu_3$cell_type[na.omit(match(de_markers_Neu_bed$V4, de_markers_Neu_3$gene))])
write.table(de_markers_Neu_bed_2, "Out/de_markers_Neu_topFC1.bed",sep="\t",quote = F,row.names = F,col.names = F)

write.csv(top10_anno,row.names = F,file = paste0("Out/1_de_genes_top10_PC",PC,"_res",res,"_anno.csv"))
write.csv(de_markers,row.names = F,file = paste0("Out/1_de_genes_PC",PC,"_res",res,"_anno.csv"))
write.csv(de_markers_HBFB,row.names = F,file = paste0("Out/1_de_genes_PC",PC,"_res",res,"_anno_HBFB.csv"))
write.csv(de_markers_HBHT,row.names = F,file = paste0("Out/1_de_genes_PC",PC,"_res",res,"_anno_HBHT.csv"))

write.csv(de_markers,row.names = F,file = paste0("Out/1_de_genes_PC",PC,"_res",res,"_anno.csv"))
qsave(Seurat_object_cropped, paste0("Out/", sample, "_cropped_PC",PC,"_res",res,"_anno.qs"))
save(top10_anno, de_markers, Seurat_object_cropped, file = paste0("Out/", sample, "_cropped_PC",PC,"_res",res,"_anno.RData"))


setwd("/media/maoni/data/CZP/spatial_transcriptome/E13-23-RNA/6_seurat_bigSpotShape_new")
sample="E13-23-RNA"
PC=30
res=0.8
load(paste0("Out/", sample, "_cropped_PC",PC,"_res",res,"_anno.RData"))

top3_anno <- de_markers %>% group_by(cell_type) %>% top_n(n = 3, wt = avg_log2FC)  %>% arrange(cell_type, desc(avg_log2FC))
top5_anno <- de_markers %>% group_by(cell_type) %>% top_n(n = 5, wt = avg_log2FC)  %>% arrange(cell_type, desc(avg_log2FC))
top10_anno <- de_markers %>% group_by(cell_type) %>% top_n(n = 10, wt = avg_log2FC)  %>% arrange(cell_type, desc(avg_log2FC))


# scRNAtoolVis (averageHeatmap)
# devtools::install_github("sajuukLyu/ggunchull", type = "source")
# devtools::install_github('junjunlab/scRNAtoolVis')

library(ggunchull)
library(scRNAtoolVis)
Idents(Seurat_object_cropped) <- Seurat_object_cropped$cell_type
pdf(file = paste0("Plot/3_findmarker_heatmap_PC",PC,"_res",res,"_anno_top3.pdf"), width = 6, height = 10)
averageHeatmap(object = Seurat_object_cropped, markerGene = top3_anno$gene, assays = "SCT", annoCol = TRUE, myanCol = my_color_palette, clusterAnnoName = F)
dev.off()

top5_anno_select <- c("Ebf1", "Robo3", "Islr2", "Gad2", "Adarb2", "Mgat4c", "Nefm",
                      "Lsamp", "Adgrv1", "Sox1ot", "Pax6", "Cdon", "Myh7", "Myh6", "Afp", "Alb", 
                      "Ptn", "Igf1", "Trpm3", "Ranbp3l", "Cxcl14", "Pax3", "Epcam", "Eya1", "Myh3", "Myh8", 
                      "Col2a1", "Col11a2", "Eln", "Fbln5", "Synpo2", "Dock9")


pdf(file = paste0("Plot/3_findmarker_heatmap_PC",PC,"_res",res,"_anno_top5_select.pdf"), width = 6, height = 8)
averageHeatmap(object = Seurat_object_cropped, markerGene = top5_anno_select, assays = "SCT", annoCol = TRUE, myanCol = my_color_palette, clusterAnnoName = F)
dev.off()

annoGene <- c("Sox2ot","Sox1ot")
averageHeatmap(object = Seurat_object_cropped, markerGene = top5_anno$gene, assays = "SCT", annoCol = TRUE, myanCol = my_color_palette, clusterAnnoName = F, showRowNames = F, markGenes = annoGene)

htCol = c("#0099CC", "white", "#CC0033")

pdf(file = paste0("Plot/3_findmarker_heatmap_PC",PC,"_res",res,"_anno_top5_select_Dotplot.pdf"), width = 20, height = 15)
DotPlot(Seurat_object_cropped, features = rev(top5_anno_select), cols = c("grey", "red"), group.by = 'cell_type') + 
  coord_flip() + 
  theme(axis.text.x = element_text(hjust = 1, angle = 45), axis.title = element_blank()) +
  ggtitle(label = 'top5_anno_select')
dev.off()


### cell_cycle
cell_cycle_genes <- read.csv("/media/maoni/data/Reference/cell_cycle/tinyatlas-master/cell_cycle/Mus_musculus.csv")
all_gene_anno <- read.table("/media/maoni/data/Reference/mouse/GRCm38_mm10/ensembl/All_gene_Anno.bed", sep="\t", header = FALSE)
cell_cycle_genes_s <- all_gene_anno$V4[match(cell_cycle_genes$geneID[cell_cycle_genes$phase=="S"], all_gene_anno$V7)]
cell_cycle_genes_g2m <- all_gene_anno$V4[match(cell_cycle_genes$geneID[cell_cycle_genes$phase=="G2/M"], all_gene_anno$V7)]
s.genes <- cell_cycle_genes_s
g2m.genes <- cell_cycle_genes_g2m

score_cc <- function(Seurat_object_cropped) {
  Seurat_object_cropped <- CellCycleScoring(Seurat_object_cropped, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  Seurat_object_cropped@meta.data$CC.Diff <- Seurat_object_cropped@meta.data$S.Score - Seurat_object_cropped@meta.data$G2M.Score
  return(Seurat_object_cropped)
}
Seurat_object_cropped_cc <- score_cc(Seurat_object_cropped)
FeatureScatter(Seurat_object_cropped_cc, "G2M.Score", "S.Score", group.by = "Phase", pt.size = .1) + coord_fixed(ratio = 1)
Seurat_object_cropped_cc <- ScaleData(Seurat_object_cropped_cc, features = c(s.genes, g2m.genes))
Seurat_object_cropped_cc <- RunPCA(Seurat_object_cropped_cc, features = c(s.genes, g2m.genes))

pdf(paste0("Plot/E13_cell_cycle_and_DotPlot.pdf"), height = 7, width = 20)
p1 <- DimPlot(Seurat_object_cropped_cc) + theme(aspect.ratio = 1)
p2 <- SpatialDimPlot(Seurat_object_cropped_cc, label = TRUE, label.size = 3, pt.size.factor = 1.6, stroke = NA) + theme(aspect.ratio = 1) # the spatial plot
p3 <- DotPlot(Seurat_object_cropped_cc, features = c("S.Score", "G2M.Score"), cols = c("grey", "red"), group.by = 'cell_type') + 
  coord_flip() + 
  theme(axis.text.x = element_text(hjust = 1, angle = 45), axis.title = element_blank()) +
  ggtitle(label = 'G2M.S.Score_DotPlot')
print(p1 + p2 + p3)
dev.off()


#### highlight markers in Fig1
gene_example_fig1 <- c("Sorcs3", "Zfpm2")
spatial_p_feature <- list()
for(i in 1:length(gene_example_fig1)){
  DEG_spe_select <- gene_example_fig1[i]
  spatial_p_feature[[DEG_spe_select]] <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = DEG_spe_select, alpha = c(1, 1), ncol = 1, pt.size.factor = 1.6, stroke = NA) + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))[c(1:6,9:11)]) + theme(aspect.ratio = 1)
  spatial_p_feature[[DEG_spe_select]]$layers[[1]]$aes_params=c(spatial_p_feature[[DEG_spe_select]]$layers[[1]]$aes_params, shape=22)
}
pdf(file = paste0("Plot/5_Spatial_PC",PC,"_res",res,"_final_gene_example_fig1.pdf"), width = 2*7, height = 4*7)
p <- ggarrange(plotlist = spatial_p_feature, ncol = 2, nrow = 4)
print(p)
dev.off()






# highlight marker genes 
source("/media/maoni/data/R_functions_ST_E13.R")
setwd("/media/maoni/data/CZP/spatial_transcriptome/E13-23-RNA/6_seurat_bigSpotShape_new")
PC=30
res=0.8
sample="E13-23-RNA"
Seurat_object_cropped_magic <- readRDS(paste0("Out/", sample, "_cropped_magic.rds"))

# overlap with markers in paper
gene_list_paper <- c("Il1rapl2", "Meox2", "Tgfb2", "Adamts9", "Postn", "Ror1", "Runx2", "Twist2", "Prrx1", "Wt1", "Mylk", "Ednra", "Sox9", 
                     "Foxp2", "Col2a1", "Col9a1", "Col11a1", "Pax9", "Ntng1", "Car10", "Epcam", "Trp63", "Grhl2", "Pth2r", "Fabp7", "Pax3", 
                     "Fzd10", "Hes5", "Gpc5", "Smoc1", "Prmt8", "Gadd45g", "Cdkn1c", "Btg2", "Nkx6-3", "Nrn1", "Slc17a6", "Grem2", "Slc17a6", 
                     "Sox1", "Olig2", "Nkx2-1", "Fgf8", "En2", "Fgf15", "Fgf17", "Pax5", "Neb", "Myh3", "Tpm2", "Acta2", "Foxb1", "Scube2", "Prtg", 
                     "Pax2", "Slc6a5", "Il23a", "Bmpr1a", "Prtg", "Col1a1", "Camk1d", "Rbm8a", "Pax2", "Slc6a5", "Fut9", "Id4", "Pcdh19", "Cdon", 
                     "Emx1", "Ptprb", "Pecam1", "Vwf", "Klhl4", "Hbegf", "Egfl7", "Itga11", "Atp1a2", "Lamc3", "Epha7", "Snca", 
                     "Slc4a1", "Kel", "Plp1", "Cdh19", "Syt13", "Shox2", "Ptprr", "Pcbp3", "Msx1", "Fgf10", "Wnt5a", "Lmx1b", 
                     "Abcb4", "Dlx1", "Dlx2", "Neurod2", "Tiam2", "Afp", "Alb", "Apoa2", "Pik3c2g", "Hoga1", 
                     "Tfeb", "Shh", "Slit1", "Ntn1", "Apoe", "Lyz2", "Selenop", "Ptprc", "Ly86", "Ctss", "Sostdc1", "Htr2c", "Kcne2", "Ttr", "Slit2", 
                     "Slit3", "Chat", "Myl2", "Myocd", "Hcn4", "Ctnna3", "Ryr2", "Tbx20", "Pf4", "Itgb3", "Itga2b", "Ppbp", "Cd226", "Tyr", "Trpm1", 
                     "Pmel", "Cryba1", "Ngp", "S100a8")
# Fan Rong
gene_list_paper <- c("Hand2", "Ascl1", "Kcnq3", "Foxa2", "Hoxa9", "Cacna2d1", "Nr2e1", "Six1", "Dlx5", "Otx2", "Six2", "Foxg1", "Ackr3", "Ina", "Pou3f2", 
                     "Otx2", "Ptprn2", "En2", "Pou3f3", "Foxg1", "Neurog1", "Neurod6", "Pou2f2", "Snhg11", "Kndc1", "Mef2c", "Thra", "Fam155a", "Six6", 
                     "Six2", "Foxl1", "Adgra2", "Pax6", "Myt1l", "Nrxn2", "Rbfox3", "Myh6", "Myh7", "Gata4", "Hemgn", "Fgf9", "Fgfr1", "Fgfr2", 
                     "Fgfr3", "Ror2", "Apoa2", "Sptb", "Gata2", "Nfe2", "Pax6", "Sox6", "Olig2", "Ascl1", "Ascl1", "Neurog2", "Sox2", "Sox1", 
                     "Sox10", "Sox2", "Pax6", "Ntng1", "Car10", "Car3", "Syt8", "Foxc2", "Hoxd11", "Nkx6-1", "Hoxa7", "Hoxa9", "Hoxa10", "Hoxa3", "Hoxa4", 
                     "Hoxa5", "Tbx3", "Tbx5")
# Peng Guangdun
gene_list_paper <- c("Mef2c", "Thra", "Nfix", "Foxp1", "Pou3f1", "Sox11", "Pou3f3", "Satb2", "Neurod2", "Neurod6", "Pou2f2", "Bcl11b", "Nfib", "Pou2f1", 
                     "Tbr1", "Meis2", "Nfia", "Notch2", "Eomes", "Notch3", "Meis1", "Sox6", "Pax6", "Notch1", "Neurog1", "Sox9", "Smad3", "Sox2", "Dlx2")



# add marker genes
gene_list_paper <- c("Nrxn1", "Syt1", "Lingo2", "Col11a1", "Foxj1", "Gstm1", "Ebf2", "Nrk")

group_index <- ceiling(seq_along(gene_list_paper) / 10)
gene_list <- split(gene_list_paper, group_index)

for(i in 1:length(gene_list)){
  
  cluster_genes <- gene_list[[i]]
  spatial_p_feature <- list()
  
  for(j in 1:length(cluster_genes)){
    
    spatial_p_feature[[j]] <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = cluster_genes[j], alpha = c(1, 1), ncol = 1, pt.size.factor = 1.6, stroke = NA) + theme(aspect.ratio = 1) + 
      scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")))
    spatial_p_feature[[j]]$layers[[1]]$aes_params=c(spatial_p_feature[[j]]$layers[[1]]$aes_params, shape=22)
    spatial_p_feature[[j]]
  }
  pdf(file = paste0("Plot/markers_highlight/marker_genes_in_Addhighlightgenes_cluster_", i, ".pdf"), width = 5*5, height = 2*5)
  p <- ggarrange(plotlist = spatial_p_feature, ncol = 5, nrow = 2)
  print(p)
  dev.off()
}



# adjust color
gene_list_paper <- c("Afp", "Alb", "Apoa2", "Snhg11")

group_index <- ceiling(seq_along(gene_list_paper) / 10)
gene_list <- split(gene_list_paper, group_index)

for(i in 1:length(gene_list)){
  # i <- 1
  
  cluster_genes <- gene_list[[i]]
  spatial_p_feature <- list()
  for(j in 1:length(cluster_genes)){
    j <- 4
    pos <- which(rownames(Seurat_object_cropped_magic)==cluster_genes[j])
    exp <- Seurat_object_cropped_magic@assays$MAGIC_SCT$data[pos,]
    Seurat_object_cropped_magic@assays$MAGIC_SCT$data[pos,which(exp > 2.5)] <- 2.5
    
    spatial_p_feature[[j]] <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = cluster_genes[j], alpha = c(1, 1), ncol = 1, pt.size.factor = 1.6, stroke = NA) + theme(aspect.ratio = 1) + 
      scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))[c(1:3, 6:11)])
      #scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")))
    spatial_p_feature[[j]]$layers[[1]]$aes_params=c(spatial_p_feature[[j]]$layers[[1]]$aes_params, shape=22)
    spatial_p_feature[[j]]
  }
  
  pdf(file = paste0("Plot/markers_highlight/marker_genes_in_Adjustcolor_cluster_", i, ".pdf"), width = 5*5, height = 2*5)
  p <- ggarrange(plotlist = spatial_p_feature, ncol = 5, nrow = 2)
  print(p)
  dev.off()
}








