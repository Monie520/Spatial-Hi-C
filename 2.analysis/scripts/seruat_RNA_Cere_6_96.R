source("/media/maoni/data/R_functions_ST_lung.R")
work_path <- "/media/maoni/data/CZP/spatial_transcriptome/Cere-RNA-6"
sample <- "Cere-RNA-6"

# 0 file preparing
set.seed(123)
setwd(work_path)
dir.create("6_seurat_bigSpotShape_new")
setwd(paste0(work_path, "/6_seurat_bigSpotShape_new"))
dir.create("Img")
dir.create("Matrix")
dir.create("Out")
dir.create("Plot")

# 0 Load data from st pipe
# run this perl cmd
# perl -lane 'if(/gene_id\s.(.+?).\;.+?gene_name\s.(.+?).\;.+?gene_biotype\s.(.+?).\;/){print "$1\t$2\t$3"}' Mus_musculus.GRCm38.102.chr.gtf | sort -u > gGRCm38.102.gene.ids.names
gene.names <- read.table("/media/maoni/data/Reference/mouse/GRCm38_mm10/ensembl/gGRCm38.102.gene.ids.names", sep = "\t")
# remove unnecessary pseudogene, rRNA, snoRNA, tRNA, ribozyme
removed.genes1 <- gene.names$V2[grep("pseudogene|tRNA|rRNA|snoRNA|ribozyme|snRNA|scaRNA|misc_RNA|scRNA|sRNA", gene.names$V3)]
removed.genes2 <- gene.names$V2[grep("mt-|Rps|Rpl|Hba|Hbb", gene.names$V2)]  # 1829
removed.genes <- unique(c(removed.genes1, removed.genes2))  # 38300

setwd(work_path)
allCounts <- readRDS(paste0("4_zUMIs/", sample, ".dgecounts.rds"))
count_matrix <- allCounts$umicount$inex$all
count_matrix <- as.matrix(count_matrix)
dim(count_matrix) # 35193  9216

raw_stat = read.table(paste0("3_umi_tools/", sample, ".debarcoded_passed_reads_stat.csv"),sep =",", header = TRUE, dec =".", stringsAsFactors = F)
raw_stat_cell = paste(raw_stat$iA, raw_stat$iB, sep="x")
raw_stat_barcode = paste0(raw_stat$bc_B, raw_stat$bc_A)
colnames(count_matrix) <- raw_stat_cell[match(colnames(count_matrix),raw_stat_barcode)]
rownames(count_matrix) <- gene.names$V2[match(rownames(count_matrix), gene.names[,1])]
dim(count_matrix)

ID <- as.character(rownames(count_matrix))
d <- duplicated(ID)
ID <- factor(ID, levels = unique(ID))
count_matrix <- rowsum(as.matrix(count_matrix), ID, reorder = FALSE, na.rm = TRUE)
removed.gene <- intersect(rownames(count_matrix), removed.genes)
count_matrix <- count_matrix[-match(removed.gene, rownames(count_matrix)),]  #25356  9216

# 1 prepare cropped matrix
stat_cropped = read.table(paste0("4_zUMIs/", sample,".debarcoded_passed_reads_stat_cropped.csv"),sep =",", header = TRUE, dec =".", stringsAsFactors = F)
iA_iB <- paste(stat_cropped$iA, stat_cropped$iB, sep="x")
count_matrix_cropped <- count_matrix[,na.omit(match(iA_iB,colnames(count_matrix)))]
dim(count_matrix_cropped)

# 2 Seurat object
# See Load10X_Spatial; Read10X_Image
# See satijalab/seurat/issue/3539 4993的
assay = "Spatial"
slice = sample
Seurat_object_cropped = CreateSeuratObject(counts = count_matrix_cropped, project = sample, assay = assay, min.cells = 3, min.features = 0) # min.cells min.features
image.dir = "/media/maoni/data/CZP/spatial_transcriptome/seurat_yuhao/core/96*96_blank_img" # "./Img"
image.nam = paste0(sample,"_fix.png") # "grey_pixel_1080p.png"
coord.nam = "combine_barcode.round2round1_index1_index2.Seurat_2.txt"
image <- readPNG(source = file.path(image.dir, image.nam))[,,1:3]
scale.factors <- c("tissue_hires_scalef"=1, "fiducial_diameter_fullres"=1, "tissue_lowres_scalef"=1)
tissue.positions <- read.table(file = file.path(image.dir,coord.nam), col.names = c("barcodes", "tissue", "row", "col", "imagerow", "imagecol"), header = FALSE, as.is = TRUE, row.names = 1)
spot.radius <- 0.015 # estiamte:(0.13)*50/410/2
image <- new(Class = "VisiumV1", image = image, scale.factors = scalefactors(spot = scale.factors[1], fiducial = scale.factors[2], hires = scale.factors[1], lowres = scale.factors[3]), coordinates = tissue.positions, spot.radius = spot.radius)
image <- image[Cells(Seurat_object_cropped)]
DefaultAssay(object = image) <- assay
Seurat_object_cropped[[slice]] <- image


# 3 Data quality
setwd(paste0(work_path, "/6_seurat_bigSpotShape_new"))

Seurat_object_cropped[["percent.MT"]] <- PercentageFeatureSet(Seurat_object_cropped, pattern = "^mt-")
Seurat_object_cropped[["percent.RP"]] <- PercentageFeatureSet(Seurat_object_cropped, pattern = "^Rps|^Rpl")
Seurat_object_cropped[["percent.HB"]] <- PercentageFeatureSet(Seurat_object_cropped, pattern = "^Hba|^Hbb")

plot0 <- VlnPlot(Seurat_object_cropped, features = c("nCount_Spatial", "nFeature_Spatial", "percent.MT", "percent.RP", "percent.HB"), ncol = 5, pt.size = 0)
plot0

plot1 <- FeatureScatter(Seurat_object_cropped, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plot2 <- FeatureScatter(Seurat_object_cropped, feature1 = "percent.MT", feature2 = "nFeature_Spatial")
plot1 + plot2

plot3 <- SpatialFeaturePlot(Seurat_object_cropped, features = "nCount_Spatial", pt.size.factor = 0.9, stroke = NA) + theme(legend.position = "top")
plot3$layers[[1]]$aes_params=c(plot3$layers[[1]]$aes_params, shape=22)
plot4 <- SpatialFeaturePlot(Seurat_object_cropped, features = "nFeature_Spatial", pt.size.factor = 0.9, stroke = NA) + theme(legend.position = "top")
plot4$layers[[1]]$aes_params=c(plot4$layers[[1]]$aes_params, shape=22)
plot3 + plot4

pdf("Plot/1_QC_check_vlnplot.pdf", width = 10, height = 4)
print(plot0)
dev.off()
pdf("Plot/1_QC_check_spatial.pdf", width = 10, height = 8)
print(plot1 + plot2 + plot3 + plot4)
dev.off()

save(Seurat_object_cropped, file = paste0("Matrix/", sample, ".cropped.RData"))
load(paste0("Matrix/", sample, ".cropped.RData"))


# 4 Run PCA
# min_nFeature_Spatial <-quantile(Seurat_object_cropped$nFeature_Spatial, probs = 0.01, na.rm = TRUE)
# max_nFeature_Spatial <-quantile(Seurat_object_cropped$nFeature_Spatial, probs = 0.99, na.rm = TRUE)
# Seurat_object_cropped <- subset(Seurat_object_cropped, subset = nFeature_Spatial > 250)
Seurat_object_cropped <- SCTransform(Seurat_object_cropped, assay = "Spatial", verbose = FALSE)

Seurat_object_cropped <- RunPCA(Seurat_object_cropped, assay = "SCT", verbose = FALSE)
ElbowPlot(Seurat_object_cropped, ndims=50)

#### PCA拐点定量识别
pct <- Seurat_object_cropped[["pca"]]@stdev / sum(Seurat_object_cropped[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
ElbowPlot(Seurat_object_cropped)$data %>% ggplot() +
  geom_point(aes(x = dims,y = stdev)) +
  geom_vline(xintercept = pc.use, color = "darkred") +
  theme_bw() + labs(title = "Elbow plot: quantitative approach")

# 5 Cluster
dim_p <- list()
spatial_p <- list()
pcNO <- 20
for(res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1, 1.2, 1.5, 2)){
  Seurat_object_cropped <- FindNeighbors(Seurat_object_cropped, reduction = "pca", dims = 1:pcNO)
  Seurat_object_cropped <- FindClusters(Seurat_object_cropped, verbose = FALSE, resolution = res)
  Seurat_object_cropped <- RunUMAP(Seurat_object_cropped, reduction = "pca", dims = 1:pcNO)
  dim_p[[paste0("res",res)]] <- DimPlot(Seurat_object_cropped, reduction = "umap", label = TRUE) + ggtitle(paste0("res",res)) # the UMAP plot
  spatial_p[[paste0("res",res)]] <- SpatialDimPlot(Seurat_object_cropped, label = TRUE, label.size = 3, pt.size.factor = 0.9, stroke = NA) + ggtitle(paste0("res",res)) # the spatial plot
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
resNO <- 0.5
for(pcNO in seq(10,50,by=5)){
  Seurat_object_cropped <- FindNeighbors(Seurat_object_cropped, reduction = "pca", dims = 1:pcNO)
  Seurat_object_cropped <- FindClusters(Seurat_object_cropped, verbose = FALSE, resolution = resNO)
  Seurat_object_cropped <- RunUMAP(Seurat_object_cropped, reduction = "pca", dims = 1:pcNO)
  dim_p[[paste0("PC",pcNO)]] <- DimPlot(Seurat_object_cropped, reduction = "umap", label = TRUE) + ggtitle(paste0("PC",pcNO)) # the UMAP plot
  spatial_p[[paste0("PC",pcNO)]] <- SpatialDimPlot(Seurat_object_cropped, label = TRUE, label.size = 3, pt.size.factor = 0.9, stroke = NA) + ggtitle(paste0("PC",pcNO)) # the spatial plot
  spatial_p[[paste0("PC",pcNO)]]$layers[[1]]$aes_params=c(spatial_p[[paste0("PC",pcNO)]]$layers[[1]]$aes_params, shape=22)
}

pdf(file = paste("Plot/2_UMAP_clusters_res",resNO,".pdf",sep =""), width = 28, height = ceiling(length(dim_p)/4)*7)
ggarrange(plotlist = dim_p, ncol = 4, nrow = ceiling(length(dim_p)/4))
dev.off()

pdf(file = paste("Plot/2_Spatial_clusters_res",resNO,".pdf",sep =""), width = 28, height = ceiling(length(spatial_p)/4)*7)
ggarrange(plotlist = spatial_p, ncol = 4, nrow = ceiling(length(spatial_p)/4))
dev.off()

# Select PC and res
PC=20
res=2.5
Seurat_object_cropped <- FindNeighbors(Seurat_object_cropped, reduction = "pca", dims = 1:PC)
Seurat_object_cropped <- FindClusters(Seurat_object_cropped, verbose = FALSE, resolution = res)
Seurat_object_cropped <- RunUMAP(Seurat_object_cropped, reduction = "pca", dims = 1:PC)
umap_Spa <- list()
umap_Spa[["p1"]] <- DimPlot(Seurat_object_cropped, reduction = "umap", label = TRUE)
umap_Spa[["p2"]] <- SpatialDimPlot(Seurat_object_cropped, label = TRUE, label.size = 3, pt.size.factor = 0.9, stroke = NA)
umap_Spa[["p2"]]$layers[[1]]$aes_params=c(umap_Spa[["p2"]]$layers[[1]]$aes_params, shape=22)
umap_Spa[["p3"]] <- DimPlot(Seurat_object_cropped, reduction = "umap", label = TRUE)
umap_Spa[["p4"]] <- SpatialDimPlot(Seurat_object_cropped, label = TRUE, label.size = 3, pt.size.factor = 0.9, stroke = NA)
umap_Spa[["p4"]]$layers[[1]]$aes_params=c(umap_Spa[["p4"]]$layers[[1]]$aes_params, shape=22)

pdf(file = paste0("Plot/2_UMAP_clusters_PC",PC,"_res",res,".pdf"), width = 28, height = ceiling(length(umap_Spa)/4)*7)
p <- ggarrange(plotlist = umap_Spa, ncol = 4, nrow = ceiling(length(umap_Spa)/4))
print(p)
dev.off()

spatial_p_cluster <- list()
for(i in sort(unique(Idents(Seurat_object_cropped)))){
  spatial_p_cluster[[paste0("idents_",i)]] <- SpatialDimPlot(Seurat_object_cropped, cells.highlight = CellsByIdentities(object = Seurat_object_cropped, idents = i), cols.highlight = c("#DE2D26", "grey90"), facet.highlight = TRUE, ncol = 1, pt.size.factor = 0.9, stroke = NA, label.size = 12) + ggtitle(paste0("idents",i)) # the spatial plot
  spatial_p_cluster[[paste0("idents_",i)]]$layers[[1]]$aes_params=c(spatial_p_cluster[[paste0("idents_",i)]]$layers[[1]]$aes_params, shape=22)
}

pdf(file = paste0("Plot/2_Spatial_clusters_PC",PC,"_res",res,".pdf"), width = 28, height = ceiling(length(spatial_p_cluster)/4)*7)
p <- ggarrange(plotlist = spatial_p_cluster, ncol = 4, nrow = ceiling(length(spatial_p_cluster)/4))
print(p)
dev.off()

pdf_combine(input = c(paste0("Plot/2_UMAP_clusters_PC",PC,"_res",res,".pdf"), paste0("Plot/2_Spatial_clusters_PC",PC,"_res",res,".pdf")), output = paste0("Plot/2_UMAP_Spatial_clusters_PC",PC,"_res",res,".pdf"))


# 6 DE
de_markers <- FindAllMarkers(Seurat_object_cropped, only.pos = T, min.pct = 0.1, logfc.threshold = 0.25) %>% arrange(cluster, desc(avg_log2FC))  #8434
top10 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)  %>% arrange(cluster, desc(avg_log2FC))

Seurat_object_cropped <- ScaleData(Seurat_object_cropped, features = top10$gene)
pdf(paste0("Plot/3_findmarker_heatmap_PC",PC,"_res",res,".pdf"),height = 7,width = 14)
p <- DoHeatmap(Seurat_object_cropped, features = top10$gene,size = 5.5)+
  theme(plot.title = element_text(hjust = 0.5, size = 5, face = "bold"),
        axis.text=element_text(size=5,face="bold"),
        axis.title=element_text(size=5,face="bold"),
        legend.text=element_text(size=10),
        legend.title = element_blank())
print(p)
dev.off()

# # VlnPlot
# pdf(paste0("Plot/3_findmarker_vlnPlot_data_PC",PC,"_res",res,".pdf"), height = 15, width = 10)
# p <- VlnPlot(Seurat_object_cropped, split.by = "seurat_clusters", features = top10$gene, flip = T, stack = T)
# print(p)
# dev.off()

# DotPlot
pdf(paste0("Plot/3_findmarker_dotplot_data_PC",PC,"_res",res,".pdf"), height = 15, width = 15)
p <- DotPlot(Seurat_object_cropped, features = unique(top10$gene), group.by = 'seurat_clusters') + 
  coord_flip() + 
  scale_color_viridis() +
  ggtitle(label = 'Top10_cellmarkers')
print(p)
dev.off()

# FeaturePlot
png(paste0("Plot/3_findmarker_UMAP_data_PC",PC,"_res",res,".png"), height = ceiling(length(table(top10$cluster)))*2*480, width = 5*480)
p <- FeaturePlot(Seurat_object_cropped, features = unique(top10$gene),ncol = 5,pt.size = 0.1,order = T)
print(p)
dev.off()

write.csv(top10,row.names = F,file = paste0("Out/1_de_genes_top10_PC",PC,"_res",res,".csv"))
write.csv(de_markers,row.names = F,file = paste0("Out/1_de_genes_PC",PC,"_res",res,".csv"))
write.csv(Seurat_object_cropped@meta.data,row.names = T,file = paste0("Out/1_cell_metadata_PC",PC,"_res",res,".csv"))
save(top10, de_markers, Seurat_object_cropped, file = paste0("Out/", sample, "_cropped_PC",PC,"_res",res,".RData"))
load(paste0("Out/", sample, "_cropped_PC",PC,"_res",res,".RData"))

# 7 magic
source("/media/maoni/data/R_functions_ST_lung.R")
Seurat_object_cropped_magic <- magic.Seurat(Seurat_object_cropped)
DefaultAssay(object = Seurat_object_cropped_magic) <- "MAGIC_SCT"

png(paste0("Plot/3_findmarker_Spatial_data_PC",PC,"_res",res,".png"), height = ceiling(length(top10$cluster)/5)*480, width = 5*480)
plot <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = unique(top10$gene), alpha = c(0.1, 1), ncol = 5, pt.size.factor = 0.9) # + scale_fill_gradientn(colours = viridis::inferno(100))
plot$layers[[1]]$aes_params=c(plot$layers[[1]]$aes_params, shape=22)
plot
dev.off()

saveRDS(Seurat_object_cropped, file = paste0("Out/", sample, "_cropped_PC", PC, "_res", res, ".RDS"))
saveRDS(Seurat_object_cropped_magic, file = paste0("Out/", sample, "_cropped_magic_PC", PC, "_res", res, ".RDS"))
Seurat_object_cropped_magic <- readRDS(paste0("Out/", sample, "_cropped_magic_PC", PC, "_res", res, ".RDS"))

SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = c("Galntl6", "Kcnd2", "Gpm6b", "Nrg1", "Ttr", "Enpp2", "Lgi2", "Gfap"), alpha = c(0.1, 1), ncol = 4, pt.size.factor = 0.9)
SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = c("Prkch", "Bc1", "Itga9", "Tlk1", "C130073E24Rik", "Enpp2", "Cspp1", "Plp1", "Mbp", "Trim30c"), alpha = c(0.1, 1), ncol = 4, pt.size.factor = 0.9)


# 8 highlight genes in papers
features_highlight <- list("Purkinje" = c("Ppp1r17", "Pcp4"),
                           "Granular" = c("Gabra6", "Slc17a7", "Rbfox3"),
                           "Golgi" = c("Slc6a5", "Grm2", "Sst", "Lgi2"),
                           "MLI1" = c("Prkcd", "Sorcs3", "Ptprk", "Lypd6"),
                           "MLI2" = c("Nxph1", "Cdh22"),
                           "PLI_3 (putative Lugaro)" = c("Htr2a", "Edil3"),
                           "PLI_1_PLI_2" = c("Aldh1a3", "Slc6a5", "Klhl1"),
                           "UBC" = c("Eomes"),
                           "Bergmann" = c("Gdf10", "Gabra2"),
                           "Astrocyte" = c("Aqp4", "Gfap"),
                           "ODC" = c("Mobp", "Ppfibp1"),
                           "other" = c("Dcn", "Kcnj8", "Ttr", "Mrc1", "C1qa", "Flt1", "Foxj1"))

png(paste0("Plot/3_findmarker_Spatial_data_PC",PC,"_res",res,"_selectgene_inpapers.png"), height = ceiling(length(unique(unlist(features_highlight)))/5)*480, width = 5*480)
plot <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = unique(unlist(features_highlight)), alpha = c(1, 1), ncol = 5, pt.size.factor = 0.9, stroke = NA) + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) + theme(aspect.ratio = 1) # + scale_fill_gradientn(colours = viridis::inferno(100))
plot$layers[[1]]$aes_params=c(plot$layers[[1]]$aes_params, shape=22)
plot
dev.off()

pdf(paste0("Plot/3_findmarker_Spatial_data_PC",PC,"_res",res,"_selectgene_inpapers.pdf"), height = ceiling(length(unique(unlist(features_highlight)))/5)*7, width = 5*7)
plot <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = unique(unlist(features_highlight)), alpha = c(1, 1), ncol = 5, pt.size.factor = 0.9, stroke = NA) + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) + theme(aspect.ratio = 1) # + scale_fill_gradientn(colours = viridis::inferno(100))
plot$layers[[1]]$aes_params=c(plot$layers[[1]]$aes_params, shape=22)
plot
dev.off()



# 9 annotation
setwd("/media/maoni/data/CZP/spatial_transcriptome/Cere-RNA-6/6_seurat_bigSpotShape_new")
PC=20
res=2.5
sample="Cere-RNA-6"
load(paste0("Out/", sample, "_cropped_PC",PC,"_res",res,".RData"))
Seurat_object_cropped_magic <- readRDS(paste0("Out/", sample, "_cropped_magic.rds"))


marker_gene_list <- read.table("/media/maoni/data/CZP/spatial_transcriptome/Cere_marker_list.txt", sep = "\t", header = TRUE)
marker_gene_list <- marker_gene_list[marker_gene_list$logFC > 1, ]

top10_marker_anno <- c()
for(i in 1:length(table(Seurat_object_cropped$seurat_clusters))){
  marker_gene_cluster <- marker_gene_list[which(marker_gene_list$gene %in% top10$gene[top10$cluster==i-1]), ]
  pos <- match(marker_gene_cluster$gene, top10$gene[top10$cluster==i-1])
  marker_gene_cluster$avg_log2FC <- top10$avg_log2FC[top10$cluster==i-1][pos]
  marker_gene_cluster$cluster <- i-1
  top10_marker_anno <- rbind(top10_marker_anno, marker_gene_cluster)
}

top10_marker_anno <- top10_marker_anno %>% arrange(cluster, desc(avg_log2FC))
write.csv(top10_marker_anno,row.names = F,file = paste0("Out/1_de_genes_top10_PC",PC,"_res",res,"_anno.csv"))


anno_marker <- rep(NA,nrow(de_markers))
anno_marker[de_markers$cluster %in% c(0)] <- "undefined"
anno_marker[de_markers$cluster %in% c(1)] <- "MLI_2"
anno_marker[de_markers$cluster %in% c(2)] <- "MLI_1"
anno_marker[de_markers$cluster %in% c(3,15)] <- "Granule_1"
anno_marker[de_markers$cluster %in% c(4)] <- "Granule_2"
anno_marker[de_markers$cluster %in% c(5)] <- "Granule_4"
anno_marker[de_markers$cluster %in% c(6)] <- "ODC_1"
anno_marker[de_markers$cluster %in% c(7,9)] <- "Bergmann"
anno_marker[de_markers$cluster %in% c(8,14)] <- "Purkinje"
anno_marker[de_markers$cluster %in% c(10)] <- "ODC_2"
anno_marker[de_markers$cluster %in% c(11)] <- "Fibroblast"
anno_marker[de_markers$cluster %in% c(12)] <- "Granule_3"
anno_marker[de_markers$cluster %in% c(13)] <- "Endothelial_stalk"

de_markers$cell_type <- factor(anno_marker,levels = c("Purkinje", "Bergmann", "Granule_1", "Granule_2", "Granule_3", "Granule_4", "ODC_1", "ODC_2", "MLI_1", "MLI_2", "Fibroblast", "Endothelial_stalk", "undefined"))

anno <- rep(NA,ncol(Seurat_object_cropped))
anno[Seurat_object_cropped$seurat_clusters %in% c(0)] <- "undefined"
anno[Seurat_object_cropped$seurat_clusters %in% c(1)] <- "MLI_2"
anno[Seurat_object_cropped$seurat_clusters %in% c(2)] <- "MLI_1"
anno[Seurat_object_cropped$seurat_clusters %in% c(3,15)] <- "Granule_1"
anno[Seurat_object_cropped$seurat_clusters %in% c(4)] <- "Granule_2"
anno[Seurat_object_cropped$seurat_clusters %in% c(5)] <- "Granule_4"
anno[Seurat_object_cropped$seurat_clusters %in% c(6)] <- "ODC_1"
anno[Seurat_object_cropped$seurat_clusters %in% c(7,9)] <- "Bergmann"
anno[Seurat_object_cropped$seurat_clusters %in% c(8,14)] <- "Purkinje"
anno[Seurat_object_cropped$seurat_clusters %in% c(10)] <- "ODC_2"
anno[Seurat_object_cropped$seurat_clusters %in% c(11)] <- "Fibroblast"
anno[Seurat_object_cropped$seurat_clusters %in% c(12)] <- "Granule_3"
anno[Seurat_object_cropped$seurat_clusters %in% c(13)] <- "Endothelial_stalk"

Seurat_object_cropped$cell_type <- factor(anno, levels = c("Purkinje", "Bergmann", "Granule_1", "Granule_2", "Granule_3", "Granule_4", "ODC_1", "ODC_2", "MLI_1", "MLI_2", "Fibroblast", "Endothelial_stalk", "undefined"))

# 13 cell types
source("/media/maoni/data/Zhanglin_Junjie/ColorPalettes.R")
my_color_palette <- c("#F3153C", "#FDBB15", "#E4A5F6", "#CE6CF5", "#BC18B5", "#1312AC", "#4F9D16", "#20B971", "#7395BF", "#1380DD", "#584B5F", "#BE8476", "#CDC5C1")
scales::show_col(my_color_palette)
names(my_color_palette) <- c("Purkinje", "Bergmann", "Granule_1", "Granule_2", "Granule_3", "Granule_4", "ODC_1", "ODC_2", "MLI_1", "MLI_2", "Fibroblast", "Endothelial_stalk", "undefined")

umap_Spa <- list()
umap_Spa[["p1"]] <- DimPlot(Seurat_object_cropped, reduction = "umap", group.by = "cell_type", label = FALSE, pt.size = 0.5, cols = my_color_palette) + theme(aspect.ratio = 1)
umap_Spa[["p2"]] <- SpatialDimPlot(Seurat_object_cropped, group.by = "cell_type", label = FALSE, label.size = 3, pt.size.factor = 0.9, stroke = NA, cols = my_color_palette) + theme(aspect.ratio = 1)
umap_Spa[["p2"]]$layers[[1]]$aes_params=c(umap_Spa[["p2"]]$layers[[1]]$aes_params, shape=22)
umap_Spa[["p3"]] <- DimPlot(Seurat_object_cropped, reduction = "umap", group.by = "cell_type", label = FALSE, pt.size = 0.5, cols = my_color_palette)+ theme(aspect.ratio = 1)
umap_Spa[["p4"]] <- SpatialDimPlot(Seurat_object_cropped, group.by = "cell_type", label = FALSE, label.size = 3, pt.size.factor = 0.9, stroke = NA, cols = my_color_palette)+ theme(aspect.ratio = 1)
umap_Spa[["p4"]]$layers[[1]]$aes_params=c(umap_Spa[["p4"]]$layers[[1]]$aes_params, shape=22)

pdf(file = paste0("Plot/2_UMAP_clusters_PC",PC,"_res",res,"_anno.pdf"), width = 28, height = ceiling(length(umap_Spa)/4)*7)
p <- ggarrange(plotlist = umap_Spa, ncol = 4, nrow = ceiling(length(umap_Spa)/4))
print(p)
dev.off()

spatial_p_cluster <- list()
for(i in names(table(Seurat_object_cropped$cell_type))){
  Seurat_object_cropped_subset <- Seurat_object_cropped[, Seurat_object_cropped$cell_type==i]
  spatial_p_cluster[[i]] <- SpatialDimPlot(Seurat_object_cropped_subset, group.by = "cell_type", crop = FALSE, label = FALSE, label.size = 3, pt.size.factor = 0.9, stroke = NA,  label.box = FALSE, cols = my_color_palette) + theme(aspect.ratio = 1) + ggtitle(i) + theme(legend.position="none") # the spatial plot
  spatial_p_cluster[[i]]$layers[[1]]$aes_params=c(spatial_p_cluster[[i]]$layers[[1]]$aes_params, shape=22)
}

pdf(file = paste0("Plot/2_Spatial_clusters_PC",PC,"_res",res,"_anno.pdf"), width = 28, height = ceiling(length(spatial_p_cluster)/4)*7)
p <- ggarrange(plotlist = spatial_p_cluster, ncol = 4, nrow = ceiling(length(spatial_p_cluster)/4))
print(p)
dev.off()

pdf_combine(input = c(paste0("Plot/2_UMAP_clusters_PC",PC,"_res",res,"_anno.pdf"), paste0("Plot/2_Spatial_clusters_PC",PC,"_res",res,"_anno.pdf")), output = paste0("Plot/2_UMAP_Spatial_clusters_PC",PC,"_res",res,"_anno.pdf"))


de_markers <- de_markers %>% group_by(cell_type) %>% arrange(cell_type, desc(avg_log2FC))
de_markers_GLPL <- de_markers[grep("Granule_1|Granule_2|Granule_3|Granule_4|Purkinje|Bergmann", de_markers$cell_type), ]

write.csv(de_markers,row.names = F,file = paste0("Out/1_de_genes_PC",PC,"_res",res,"_anno.csv"))
write.csv(de_markers_GLPL,row.names = F,file = paste0("Out/1_de_genes_PC",PC,"_res",res,"_anno_GLPL.csv"))

All_gene_Anno_PCG <- read.table("/media/maoni/data/Reference/mouse/GRCm38_mm10/ensembl/All_gene_Anno_PCG.bed", sep = "\t", header = FALSE)
de_markers_GLPL_bed <- All_gene_Anno_PCG[na.omit(match(unique(de_markers_GLPL$gene), All_gene_Anno_PCG$V4)), ]
de_markers_GLPL_bed_2 <- data.frame(de_markers_GLPL_bed, de_markers_GLPL$cell_type[na.omit(match(de_markers_GLPL_bed$V4, de_markers_GLPL$gene))])
de_markers_GLPL_bed_2_chr16 <- de_markers_GLPL_bed_2[de_markers_GLPL_bed_2$V1=="chr16", ]

write.table(de_markers_GLPL_bed_2, "Out/de_markers_GLPL.bed",sep="\t",quote = F,row.names = F,col.names = F)
write.table(de_markers_GLPL_bed_2, "Out/de_markers_PLGL.bed",sep="\t",quote = F,row.names = F,col.names = F)
write.table(de_markers_GLPL_bed_2_chr16, "Out/de_markers_GLPL_chr16.bed",sep="\t",quote = F,row.names = F,col.names = F)

top_anno <- de_markers %>% group_by(cell_type) %>% top_n(n = 20, wt = avg_log2FC)  %>% arrange(cell_type, desc(avg_log2FC))
save(top_anno, de_markers, Seurat_object_cropped, file = paste0("Out/", sample, "_cropped_PC",PC,"_res",res,"_anno.RData"))
load(paste0("Out/", sample, "_cropped_PC",PC,"_res",res,"_anno.RData"))

# lncRNA(de_markers_GLPL_bed_2_chr16)
png("Plot/Cere_gene_highlight_de_markers_GLPL_bed_2_chr16.png", width = 5*480, height = 7*480)
p1 <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = de_markers_GLPL_bed_2_chr16$V4, alpha = c(1, 1), ncol = 5, pt.size.factor = 0.9, stroke = NA) + theme(aspect.ratio = 1) + theme(aspect.ratio = 1) # + scale_fill_gradientn(colours = viridis::inferno(100))
p1
dev.off()

# lncRNA(Gm2694)
pdf("/media/maoni/data/Gm2694/Cere_gene_highlight.pdf", width = 4, height = 12)
p1 <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = c("Gm2694", "Cbln1", "Calb1"), alpha = c(1, 1), ncol = 1, pt.size.factor = 0.9, stroke = NA) + theme(aspect.ratio = 1) + theme(aspect.ratio = 1) # + scale_fill_gradientn(colours = viridis::inferno(100))
p1
dev.off()

Idents(Seurat_object_cropped) <- Seurat_object_cropped$cell_type
p2 <- VlnPlot(Seurat_object_cropped, split.by = "cell_type", features = c("Gm2694", "Cbln1", "Calb1"), flip = T, stack = T)
p2

p3 <- DotPlot(Seurat_object_cropped, features = rev(c("Gm2694", "Cbln1", "Calb1")), cols = c("grey", "red"), group.by = 'cell_type') + coord_flip() + theme(axis.text.x = element_text(hjust = 1, angle = 45), axis.title = element_blank())
p3

pdf("/media/maoni/data/Gm2694/Cere_vlnplot_dotplot.pdf", width = 14, height = 5)
p2 + p3
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

pdf(paste0("Plot/Cere_cell_cycle_and_DotPlot.pdf"), height = 5, width = 15)
p1 <- DimPlot(Seurat_object_cropped_cc) + theme(aspect.ratio = 1)
p2 <- SpatialDimPlot(Seurat_object_cropped_cc, label = TRUE, label.size = 3, pt.size.factor = 0.9, stroke = NA) + theme(aspect.ratio = 1) # the spatial plot
p3 <- DotPlot(Seurat_object_cropped_cc, features = c("S.Score", "G2M.Score"), cols = c("grey", "red"), group.by = 'cell_type') + 
  coord_flip() + 
  theme(axis.text.x = element_text(hjust = 1, angle = 45), axis.title = element_blank()) +
  ggtitle(label = 'G2M.S.Score_DotPlot')
print(p1 + p2 + p3)
dev.off()


############################################### plot
setwd("/media/maoni/data/CZP/spatial_transcriptome/Cere-RNA-6/6_seurat_bigSpotShape_new/")
PC=20
res=2.5
sample="Cere-RNA-6"
load(paste0("Out/", sample, "_cropped_PC",PC,"_res",res,"_anno.RData"))
top_anno <- de_markers %>% group_by(cell_type) %>% top_n(n = 20, wt = avg_log2FC)  %>% arrange(cell_type, desc(avg_log2FC))

top_anno_GL <- top_anno[grep("Granule_1|Granule_2|Granule_3|Granule_4", top_anno$cell_type), ]
top_anno_GL_unique <- unique(top_anno_GL$gene)

pdf(paste0("Plot/3_findmarker_Spatial_data_PC",PC,"_res",res,"_anno_top20_GL_genes.pdf"), height = ceiling(length(top_anno_GL_unique)/5)*7, width = 5*7)
plot <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = top_anno_GL_unique, alpha = c(1, 1), ncol = 5, pt.size.factor = 0.9, stroke = NA) + theme(aspect.ratio = 1) + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) + theme(aspect.ratio = 1) # + scale_fill_gradientn(colours = viridis::inferno(100))
plot$layers[[1]]$aes_params=c(plot$layers[[1]]$aes_params, shape=22)
plot
dev.off()

### heatmap plot
library(scRNAtoolVis)
Idents(Seurat_object_cropped) <- Seurat_object_cropped$cell_type
pdf(file = paste0("Plot/3_findmarker_heatmap_PC",PC,"_res",res,"_anno_top20.pdf"), width = 10, height = 20)
AverageHeatmap(object = Seurat_object_cropped, markerGene = top_anno$gene, assays = "SCT", annoCol = TRUE, myanCol = my_color_palette, clusterAnnoName = F)
dev.off()

### volcano plot
Idents(Seurat_object_cropped) <- Seurat_object_cropped$cell_type
de_markers_cell_type <- FindAllMarkers(Seurat_object_cropped, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25) %>% arrange(cluster, desc(avg_log2FC))  #8434
top10 <- de_markers_cell_type %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)  %>% arrange(cluster, desc(avg_log2FC))

Seurat_object_cropped_subset <- Seurat_object_cropped[, grep("Granule_1|Granule_2|Granule_3", Seurat_object_cropped$cell_type)]
de_markers_GL <- de_markers_cell_type[grep("Granule_1|Granule_2|Granule_3", de_markers_cell_type$cluster), ]

markerVocalno(markers = de_markers_GL, topn = 5,labelCol = ggsci::pal_npg()(9))

p1 <- jjVolcano(diffData = de_markers_GL, topGeneN = 5)

mygene <- c('Kcnip4', 'Kcnd2', 'Dpp6', "Dpp10",'Grik2', "Gprin3", "Galntl6")

p2 <- jjVolcano(diffData = de_markers_GL, myMarkers = mygene)

pdf(file = paste0("Plot/3_findmarker_volcano_PC",PC,"_res",res,"_anno_selectgene.pdf"), width = 12, height = 8)
p1 + p2
dev.off()


### DotPlot
pdf(paste0("Plot/3_findmarker_dotplot_data_PC",PC,"_res",res,"_anno_top20.pdf"), width = 10, height = 20)
p <- DotPlot(Seurat_object_cropped, features = unique(top_anno$gene), group.by = 'cell_type') + 
  coord_flip() + 
  scale_color_viridis() +
  ggtitle(label = 'Top20_cellmarkers')
print(p)
dev.off()

### Vlnplot
pdf(paste0("Plot/3_findmarker_vlnPlot_data_PC",PC,"_res",res,"_anno_top20.pdf"), width = 10, height = 15)
p <- VlnPlot(Seurat_object_cropped, split.by = "cell_type", features = unique(top_anno$gene), flip = T, stack = T)
print(p)
dev.off()

### SpatialFeaturePlot hihglight cell_type
Seurat_object_cropped_magic$cell_type <- Seurat_object_cropped$cell_type
cell_type <- names(table(Seurat_object_cropped_magic$cell_type))
for(i in 1:length(table(Seurat_object_cropped_magic$cell_type))){
features <- top_anno$gene[top_anno$cell_type==cell_type[i]]
png(paste0("Plot/3_findmarker_Spatial_data_PC",PC,"_res",res,"_", cell_type[i], ".png"), height = ceiling(length(features)/5)*480, width = 5*480)
plot <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = features, alpha = c(0.1, 1), ncol = 5, pt.size.factor = 0.9) # + scale_fill_gradientn(colours = viridis::inferno(100))
plot$layers[[1]]$aes_params=c(plot$layers[[1]]$aes_params, shape=22)
print(plot)
dev.off()
}

### FeaturePlot
features_MLI <- c("Lypd6", "Sorcs3", "Ptprk", "Nxph1", "Cdh22")
png(paste0("Plot/3_findmarker_UMAP_data_PC",PC,"_res",res,"_MLI.png"), height = ceiling(length(features_MLI))/5*480, width = 5*480)
p <- FeaturePlot(Seurat_object_cropped, features = features_MLI, ncol = 5, pt.size = 0.1, order = T)
print(p)
dev.off()

############################################### plot（Final)


features_highlight <- list("Granular" = c("Rbfox1", "Rbfox3"),
                           "Purkinje" = c("Ppp1r17","Car8"),
                           "Bergmann" = c("Gdf10", "Hopx"),
                           "ODC" = c("Plp1", "Mbp", "Mobp"),
                           "VLMC" = c("Ptgds"),
                           "Endothelial" = c("Flt1"),
                           "MLI" = c("Sorcs3"))

color_exp <- list(Rbfox1 = c(1:7, 9:11),
                  Rbfox3 = c(1:6, 9:10),
                  Ppp1r17 = c(1:11),
                  Car8 = c(1:11),
                  Gdf10 = c(1:6, 9:11),
                  Hopx = c(1:6, 9:11),
                  Plp1 = c(1:6, 9:11),
                  Mbp = c(1:6, 9:11),
                  Mobp = c(1:6, 9:11),
                  Ptgds = c(1:6, 9:10),
                  Flt1 = c(1:6, 9:10),
                  Sorcs3 = c(1:7, 10:11))

pos <- which(rownames(Seurat_object_cropped_magic)=="Ptgds")
pos_exp <- which(Seurat_object_cropped_magic@assays$MAGIC_SCT@data[pos, ] > 3)
Seurat_object_cropped_magic@assays$MAGIC_SCT@data[pos, pos_exp] <- 3

pos <- which(rownames(Seurat_object_cropped_magic)=="Flt1")
pos_exp <- which(Seurat_object_cropped_magic@assays$MAGIC_SCT@data[pos, ] > 0.45)
Seurat_object_cropped_magic@assays$MAGIC_SCT@data[pos, pos_exp] <- 0.45

pos <- which(rownames(Seurat_object_cropped_magic)=="Sorcs3")
pos_exp <- which(Seurat_object_cropped_magic@assays$MAGIC_SCT@data[pos, ] > 0.4)
Seurat_object_cropped_magic@assays$MAGIC_SCT@data[pos, pos_exp] <- 0.4

p_list <- list()
for(i in 1:length(unlist(features_highlight))){
  # i <- 12
  p_list[[i]] <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = unlist(features_highlight)[i], alpha = c(1, 1), ncol = 1, pt.size.factor = 0.9, stroke = NA) + theme(aspect.ratio = 1) + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))[color_exp[[i]]]) + theme(aspect.ratio = 1) # + scale_fill_gradientn(colours = viridis::inferno(100))
  p_list[[i]]$layers[[1]]$aes_params=c(p_list[[i]]$layers[[1]]$aes_params, shape=22)
  p_list[[i]]
}

pdf(paste0("Plot/3_findmarker_Spatial_data_PC",PC,"_res",res,"_anno_top20_selectgene_Final.pdf"), height = ceiling(length(unique(unlist(features_highlight)))/5)*7, width = 5*7)
ggarrange(plotlist = p_list, ncol = 5, nrow = ceiling(length(p_list)/5))
dev.off()

features_highlight <- list("Granule_1" = c("Kcnip4", "Dpp6", "Kcnd2", "Ptprt"),
                           "Granule_2" = c("Grik2", "Wincr1", "Nav2", "Ntrk3", "Cntn6"),
                           "Granule_3" = c("Patj", "Slc9a9", "Cacna2d1"),
                           "Granule_4" = c("Cacna2d1"),
                           "Purkinje" = c("Ppp1r17", "Prkg1", "Kcnc3", "Car8"),
                           "Bergmann" = c("Sparcl1", "Slc1a2", "Mt3", "Hopx"),
                           "ODC_1" = c("Prkch", "Bc1", "Tlk1"),
                           "ODC_2" = c("Mbp", "Itga9", "Plp1"),
                           "MLI_1" = c("Adarb2", "Kit", "Slc24a3", "Asic2", "Cacna2d3", "Cacng2"),
                           "MLI_2" = c("S100b", "Shank1", "Pcsk1n", "Kif5b", "Camk2n1", "Rflnb"),
                           "Fibroblast" = c("Ptgds", "Apod", "Nnat", "Slc6a13", "Ranbp3l", "Slc7a11"),
                           "Endothelial_stalk" = c("Flt1", "Slco1a4", "Pecam1"))

png(paste0("Plot/3_findmarker_Spatial_data_PC",PC,"_res",res,"_anno_top20_selectgene.png"), height = ceiling(length(unique(unlist(features_highlight)))/5)*480, width = 5*480)
plot <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = unique(unlist(features_highlight)), alpha = c(1, 1), ncol = 5, pt.size.factor = 0.9, stroke = NA) + theme(aspect.ratio = 1) + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) + theme(aspect.ratio = 1) # + scale_fill_gradientn(colours = viridis::inferno(100))
plot$layers[[1]]$aes_params=c(plot$layers[[1]]$aes_params, shape=22)
plot
dev.off()

pdf(paste0("Plot/3_findmarker_Spatial_data_PC",PC,"_res",res,"_anno_top20_selectgene.pdf"), height = ceiling(length(unique(unlist(features_highlight)))/5)*7, width = 5*7)
plot <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = unique(unlist(features_highlight)), alpha = c(1, 1), ncol = 5, pt.size.factor = 0.9, stroke = NA) + theme(aspect.ratio = 1) + scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) + theme(aspect.ratio = 1) # + scale_fill_gradientn(colours = viridis::inferno(100))
plot$layers[[1]]$aes_params=c(plot$layers[[1]]$aes_params, shape=22)
plot
dev.off()


### Heatmap plot
library(scRNAtoolVis)
Idents(Seurat_object_cropped) <- factor(as.character(Seurat_object_cropped$cell_type), levels = c("Granule_1", "Granule_2", "Granule_3", "Granule_4", "Purkinje", "Bergmann", "ODC_1", "ODC_2", "MLI_1", "MLI_2", "Fibroblast", "Endothelial_stalk", "undefined"))
pdf(file = paste0("Plot/3_findmarker_heatmap_PC",PC,"_res",res,"_anno_top20_selectgene.pdf"), width = 6, height = 10)
averageHeatmap(object = Seurat_object_cropped, markerGene = rev(unlist(features_highlight)), assays = "SCT", annoCol = TRUE, myanCol = my_color_palette, clusterAnnoName = F)
dev.off()

### DotPlot
pdf(paste0("Plot/3_findmarker_dotplot_data_PC",PC,"_res",res,"_anno_top20_selectgene.pdf"), width = 10, height = 10)
p <- DotPlot(Seurat_object_cropped, features = unique(unlist(features_highlight)), group.by = 'cell_type') + 
  coord_flip() + 
  scale_color_viridis() +
  theme(axis.text.x = element_text(hjust = 1, angle = 45), axis.title = element_blank()) +
  ggtitle(label = 'Top20_cellmarkers')
print(p)
dev.off()

### Vlnplot
pdf(paste0("Plot/3_findmarker_vlnPlot_data_PC",PC,"_res",res,"_anno_top20_selectgene.pdf"), width = 10, height = 10)
p <- VlnPlot(Seurat_object_cropped, split.by = "cell_type", features = unique(unlist(features_highlight)), flip = T, stack = T)
print(p)
dev.off()

