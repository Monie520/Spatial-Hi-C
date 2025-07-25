source("/media/maoni/data/R_functions_ST_lung.R")
work_path <- "/media/maoni/data/CZP/spatial_transcriptome/Adult-RNA-23"
sample <- "Adult-RNA-23"

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
removed.genes3 <- gene.names$V2[gene.names$V3!="protein_coding"]
removed.genes <- unique(c(removed.genes1, removed.genes2, removed.genes3))  # 38300

setwd(work_path)
allCounts <- readRDS(paste0("4_zUMIs/", sample, ".dgecounts.rds"))
count_matrix <- allCounts$umicount$inex$all
count_matrix <- as.matrix(count_matrix)
dim(count_matrix) # 35193  9216

raw_stat = read.table(paste0("3_umi_tools/", sample, ".debarcoded_passed_reads_stat.csv"),sep =",", header = TRUE, dec =".", stringsAsFactors = F)
raw_stat_cell = paste(raw_stat$iB, 97-raw_stat$iA, sep="x")
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
dim(count_matrix)

# 1 prepare cropped matrix
stat_cropped = read.table(paste0("4_zUMIs/", sample,".debarcoded_passed_reads_stat_cropped.csv"),sep =",", header = TRUE, dec =".", stringsAsFactors = F)
iA_iB <- paste(stat_cropped$iB, 97-stat_cropped$iA, sep="x")
count_matrix_cropped <- count_matrix[,na.omit(match(iA_iB,colnames(count_matrix)))]
dim(count_matrix_cropped) # 18731

# 2 Seurat object
# See Load10X_Spatial; Read10X_Image
# See satijalab/seurat/issue/3539 4993
assay = "Spatial"
slice = sample
Seurat_object_cropped = CreateSeuratObject(counts = count_matrix_cropped, project = sample, assay = assay, min.cells = 3, min.features = 1) # min.cells min.features
image.dir = "/media/maoni/data/CZP/spatial_transcriptome/seurat_yuhao/core/96*96_blank_img" # "./Img"
image.nam = paste0("white_background_1080p.png") # "grey_pixel_1080p.png"
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

plot0 <- VlnPlot(Seurat_object_cropped, features = c("nCount_Spatial", "nFeature_Spatial"), ncol = 2, pt.size = 0)
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
# Seurat_object_cropped <- subset(Seurat_object_cropped, subset = nFeature_Spatial > 100)
Seurat_object_cropped <- SCTransform(Seurat_object_cropped, assay = "Spatial", verbose = FALSE)
Seurat_object_cropped <- RunPCA(Seurat_object_cropped, assay = "SCT", verbose = FALSE)
ElbowPlot(Seurat_object_cropped, ndims=50)

#### PCA拐点定量识别
pct <- Seurat_object_cropped[["pca"]]@stdev / sum(Seurat_object_cropped[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
ElbowPlot(Seurat_object_cropped, ndims=50)$data %>% ggplot() +
  geom_point(aes(x = dims,y = stdev)) +
  geom_vline(xintercept = pc.use, color = "darkred") +
  theme_bw() + labs(title = "Elbow plot: quantitative approach")

# 5. Cluster
dim_p <- list()
spatial_p <- list()
pcNO <- 30
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
PC=30
res=0.8
Seurat_object_cropped <- FindNeighbors(Seurat_object_cropped, reduction = "pca", dims = 1:PC)
Seurat_object_cropped <- FindClusters(Seurat_object_cropped, verbose = FALSE, resolution = res)
Seurat_object_cropped <- RunUMAP(Seurat_object_cropped, reduction = "pca", dims = 1:PC, min.dist = 0.1)

source("/media/maoni/data/Zhanglin_Junjie/ColorPalettes.R")
my_color_palette <- c(ArchRPalettes$stallion, ArchRPalettes$bear[1:(length(table(Seurat_object_cropped$seurat_clusters))-20)])
scales::show_col(my_color_palette)
names(my_color_palette) <- 0:(length(table(Seurat_object_cropped$seurat_clusters))-1)

umap_Spa <- list()
umap_Spa[["p1"]] <- DimPlot(Seurat_object_cropped, reduction = "umap", group.by = "seurat_clusters", label = FALSE, pt.size = 0.5, cols = my_color_palette) + theme(aspect.ratio = 1)
umap_Spa[["p2"]] <- SpatialDimPlot(Seurat_object_cropped, group.by = "seurat_clusters", label = FALSE, label.size = 3, pt.size.factor = 0.9, stroke = NA, cols = my_color_palette) + theme(aspect.ratio = 1)
umap_Spa[["p2"]]$layers[[1]]$aes_params=c(umap_Spa[["p2"]]$layers[[1]]$aes_params, shape=22)
umap_Spa[["p3"]] <- DimPlot(Seurat_object_cropped, reduction = "umap", group.by = "seurat_clusters", label = FALSE, pt.size = 0.5, cols = my_color_palette)+ theme(aspect.ratio = 1)
umap_Spa[["p4"]] <- SpatialDimPlot(Seurat_object_cropped, group.by = "seurat_clusters", label = FALSE, label.size = 3, pt.size.factor = 0.9, stroke = NA, cols = my_color_palette)+ theme(aspect.ratio = 1)
umap_Spa[["p4"]]$layers[[1]]$aes_params=c(umap_Spa[["p4"]]$layers[[1]]$aes_params, shape=22)

pdf(file = paste0("Plot/2_UMAP_clusters_PC",PC,"_res",res,"_final.pdf"), width = 28, height = ceiling(length(umap_Spa)/4)*7)
p <- ggarrange(plotlist = umap_Spa, ncol = 4, nrow = ceiling(length(umap_Spa)/4))
print(p)
dev.off()

spatial_p_cluster <- list()
for(i in names(table(Seurat_object_cropped$seurat_clusters))){
  Seurat_object_cropped_subset <- Seurat_object_cropped[, Seurat_object_cropped$seurat_clusters==i]
  spatial_p_cluster[[i]] <- SpatialDimPlot(Seurat_object_cropped_subset, group.by = "seurat_clusters", crop = FALSE, label = FALSE, label.size = 3, pt.size.factor = 0.9, stroke = NA,  label.box = FALSE, cols = my_color_palette) + theme(aspect.ratio = 1) + ggtitle(i) + theme(legend.position="none") # the spatial plot
  spatial_p_cluster[[i]]$layers[[1]]$aes_params=c(spatial_p_cluster[[i]]$layers[[1]]$aes_params, shape=22)
}

pdf(file = paste0("Plot/2_Spatial_clusters_PC",PC,"_res",res,"_final.pdf"), width = 28, height = ceiling(length(spatial_p_cluster)/4)*7)
p <- ggarrange(plotlist = spatial_p_cluster, ncol = 4, nrow = ceiling(length(spatial_p_cluster)/4))
print(p)
dev.off()

pdf_combine(input = c(paste0("Plot/2_UMAP_clusters_PC",PC,"_res",res,"_final.pdf"), paste0("Plot/2_Spatial_clusters_PC",PC,"_res",res,"_final.pdf")), output = paste0("Plot/2_UMAP_Spatial_clusters_PC",PC,"_res",res,"_final.pdf"))

saveRDS(Seurat_object_cropped, file = paste0("Out/", sample, "_cropped_PC", PC, "_res", res, ".RDS"))


### 6 annotation
Seurat_object_cropped <- readRDS(paste0("Out/", sample, "_cropped_PC", PC, "_res", res, ".RDS"))

anno <- rep(NA,ncol(Seurat_object_cropped))
anno[Seurat_object_cropped$seurat_clusters %in% c(0)] <- "Other"
anno[Seurat_object_cropped$seurat_clusters %in% c(1)] <- "HPF-CA1-so-ExN"
anno[Seurat_object_cropped$seurat_clusters %in% c(2,6)] <- "ODC"
anno[Seurat_object_cropped$seurat_clusters %in% c(3)] <- "HPF-DG-mo-ExN"
anno[Seurat_object_cropped$seurat_clusters %in% c(4)] <- "PVT-ExN"
anno[Seurat_object_cropped$seurat_clusters %in% c(5)] <- "TH-ExN"
anno[Seurat_object_cropped$seurat_clusters %in% c(7)] <- "Endo"
anno[Seurat_object_cropped$seurat_clusters %in% c(8)] <- "VLMC"
anno[Seurat_object_cropped$seurat_clusters %in% c(9)] <- "MH-Cho"
anno[Seurat_object_cropped$seurat_clusters %in% c(10)] <- "HPF-DG-sg-Granule"
anno[Seurat_object_cropped$seurat_clusters %in% c(11)] <- "TH-Astro"
anno[Seurat_object_cropped$seurat_clusters %in% c(12)] <- "LH-Astro"
anno[Seurat_object_cropped$seurat_clusters %in% c(13)] <- "HPF-CA1-sp-ExN"
anno[Seurat_object_cropped$seurat_clusters %in% c(14)] <- "ChP"
anno[Seurat_object_cropped$seurat_clusters %in% c(15)] <- "HPF-PN"
anno[Seurat_object_cropped$seurat_clusters %in% c(16)] <- "HPF-InN"
anno[Seurat_object_cropped$seurat_clusters %in% c(17)] <- "MG"
anno[Seurat_object_cropped$seurat_clusters %in% c(18)] <- "HPF-DG-po-ExN"
anno[Seurat_object_cropped$seurat_clusters %in% c(19)] <- "OPC"
anno[Seurat_object_cropped$seurat_clusters %in% c(20)] <- "Epen"
anno[Seurat_object_cropped$seurat_clusters %in% c(21)] <- "SM-Cho"

Seurat_object_cropped$cell_type <- factor(anno, levels = names(table(anno))[c(1:15, 17:21, 16)])
color_order <- c("Other", "HPF-CA1-so-ExN", "ODC", "HPF-DG-mo-ExN", "PVT-ExN", "TH-ExN", "Endo", "VLMC", "MH-Cho", 
                 "HPF-DG-sg-Granule", "TH-Astro", "LH-Astro", "HPF-CA1-sp-ExN", "ChP", "HPF-PN", "HPF-InN", "MG", 
                 "HPF-DG-po-ExN", "OPC", "Epen", "SM-Cho")

source("/media/maoni/data/Zhanglin_Junjie/ColorPalettes.R")
my_color_palette <- c(ArchRPalettes$stallion, ArchRPalettes$bear[1:(length(table(Seurat_object_cropped$seurat_clusters))-20)])
my_color_palette_2 <- c("grey80", my_color_palette[c(2:6,8:22)])
scales::show_col(my_color_palette_2)
names(my_color_palette_2) <- color_order

umap_Spa <- list()
umap_Spa[["p1"]] <- DimPlot(Seurat_object_cropped, reduction = "umap", group.by = "cell_type", label = FALSE, pt.size = 0.1, cols = my_color_palette_2) + theme(aspect.ratio = 1)
umap_Spa[["p2"]] <- SpatialDimPlot(Seurat_object_cropped, group.by = "cell_type", label = FALSE, label.size = 3, pt.size.factor = 0.9, stroke = NA, cols = my_color_palette_2) + theme(aspect.ratio = 1)
umap_Spa[["p2"]]$layers[[1]]$aes_params=c(umap_Spa[["p2"]]$layers[[1]]$aes_params, shape=22)
umap_Spa[["p3"]] <- DimPlot(Seurat_object_cropped, reduction = "umap", group.by = "cell_type", label = FALSE, pt.size = 0.1, cols = my_color_palette_2) + theme(aspect.ratio = 1)
umap_Spa[["p4"]] <- SpatialDimPlot(Seurat_object_cropped, group.by = "cell_type", label = FALSE, label.size = 3, pt.size.factor = 0.9, stroke = NA, cols = my_color_palette_2) + theme(aspect.ratio = 1)
umap_Spa[["p4"]]$layers[[1]]$aes_params=c(umap_Spa[["p4"]]$layers[[1]]$aes_params, shape=22)

pdf(file = paste0("Plot/2_UMAP_clusters_PC",PC,"_res",res,"_final_anno.pdf"), width = 28, height = ceiling(length(umap_Spa)/4)*7)
p <- ggarrange(plotlist = umap_Spa, ncol = 4, nrow = ceiling(length(umap_Spa)/4))
print(p)
dev.off()

spatial_p_cluster <- list()
for(i in names(table(Seurat_object_cropped$cell_type))){
  Seurat_object_cropped_subset <- Seurat_object_cropped[, Seurat_object_cropped$cell_type==i]
  spatial_p_cluster[[i]] <- SpatialDimPlot(Seurat_object_cropped_subset, group.by = "cell_type", crop = FALSE, label = FALSE, label.size = 3, pt.size.factor = 0.9, stroke = NA,  label.box = FALSE, cols = my_color_palette_2) + theme(aspect.ratio = 1) + ggtitle(i) + theme(legend.position="none") # the spatial plot
  spatial_p_cluster[[i]]$layers[[1]]$aes_params=c(spatial_p_cluster[[i]]$layers[[1]]$aes_params, shape=22)
}

pdf(file = paste0("Plot/2_Spatial_clusters_PC",PC,"_res",res,"_final_anno.pdf"), width = 28, height = ceiling(length(spatial_p_cluster)/4)*7)
p <- ggarrange(plotlist = spatial_p_cluster, ncol = 4, nrow = ceiling(length(spatial_p_cluster)/4))
print(p)
dev.off()

pdf_combine(input = c(paste0("Plot/2_UMAP_clusters_PC",PC,"_res",res,"_final_anno.pdf"), paste0("Plot/2_Spatial_clusters_PC",PC,"_res",res,"_final_anno.pdf")), output = paste0("Plot/2_UMAP_Spatial_clusters_PC",PC,"_res",res,"_final_anno.pdf"))


saveRDS(Seurat_object_cropped, file = paste0("Out/", sample, "_cropped_PC", PC, "_res", res, ".RDS"))
qsave(Seurat_object_cropped, paste0("Out/", sample, "_cropped_PC", PC, "_res", res, ".qs"))

Seurat_object_cropped <- readRDS(paste0("Out/", sample, "_cropped_PC", PC, "_res", res, ".RDS"))

library(reticulate)
use_python("/home/maoni/miniconda3/envs/r4.3/bin/python")
seurat2scanpy(Seurat_object_cropped, ann.X = "Spatial-data", ann.raw.X = NULL, h5ad_path = paste0("Out/", sample, "_Seurat_object_cropped.h5ad"))



#### 7 DE
de_markers <- FindAllMarkers(Seurat_object_cropped, only.pos = T, min.pct = 0.1, logfc.threshold = 0.25) %>% arrange(cluster, desc(avg_log2FC))  #8434
top10 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)  %>% arrange(cluster, desc(avg_log2FC))
top50 <- de_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)  %>% arrange(cluster, desc(avg_log2FC))
top100 <- de_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)  %>% arrange(cluster, desc(avg_log2FC))
top_FC <- de_markers[de_markers$avg_log2FC > 1, ]
top_FC <- top_FC %>% group_by(cluster) %>% arrange(cluster, desc(avg_log2FC))

Seurat_object_cropped <- ScaleData(Seurat_object_cropped, features = top10$gene)
pdf(paste0("Plot/3_findmarker_heatmap_PC",PC,"_res",res,".pdf"),height = 14,width = 14)
p <- DoHeatmap(Seurat_object_cropped, features = top10$gene,size = 5.5)+
  theme(plot.title = element_text(hjust = 0.5, size = 5, face = "bold"),
        axis.text=element_text(size=5,face="bold"),
        axis.title=element_text(size=5,face="bold"),
        legend.text=element_text(size=10),
        legend.title = element_blank())
print(p)
dev.off()

# VlnPlot
pdf(paste0("Plot/3_findmarker_vlnPlot_data_PC",PC,"_res",res,".pdf"), height = 15, width = 10)
p <- VlnPlot(Seurat_object_cropped, split.by = "seurat_clusters", features = top10$gene, flip = T, stack = T)
print(p)
dev.off()

# DotPlot
pdf(paste0("Plot/3_findmarker_dotplot_data_PC",PC,"_res",res,".pdf"), height = 30, width = 15)
p <- DotPlot(Seurat_object_cropped, features = unique(top10$gene), group.by = 'seurat_clusters') +
  coord_flip() +
  scale_color_viridis() +
  ggtitle(label = 'Top10_cellmarkers')
print(p)
dev.off()

# # FeaturePlot
# png(paste0("Plot/3_findmarker_UMAP_data_PC",PC,"_res",res,".png"), height = ceiling(length(table(top10$cluster)))*2*480, width = 5*480)
# p <- FeaturePlot(Seurat_object_cropped, features = unique(top10$gene),ncol = 5,pt.size = 0.1,order = T)
# print(p)
# dev.off()

write.csv(top10,row.names = F,file = paste0("Out/1_de_genes_top10_PC",PC,"_res",res,".csv"))
write.csv(top100,row.names = F,file = paste0("Out/1_de_genes_top100_PC",PC,"_res",res,".csv"))
write.csv(top_FC,row.names = F,file = paste0("Out/1_de_genes_top_FC_PC",PC,"_res",res,".csv"))
write.csv(de_markers,row.names = F,file = paste0("Out/1_de_genes_PC",PC,"_res",res,".csv"))
write.csv(Seurat_object_cropped@meta.data,row.names = T,file = paste0("Out/1_cell_metadata_PC",PC,"_res",res,".csv"))
save(top10, de_markers, Seurat_object_cropped, file = paste0("Out/", sample, "_cropped_PC",PC,"_res",res,".RData"))
load(paste0("Out/", sample, "_cropped_PC",PC,"_res",res,".RData"))

# # markers in paper
# marker_2021_nature <- read.csv("/media/maoni/data/CZP/spatial_transcriptome/Brain_marker_list_2021_nature_all.txt", sep="\t", header = TRUE)
# cluster <- names(table(top_FC$cluster))
# for(j in 1:length(cluster)){
#   top_FC_cluster <- top_FC[top_FC$cluster==cluster[j], ]
#   top_FC_cluster_anno <- c()
#   for(i in 1:length(top_FC_cluster$gene)){
#     anno <- cbind(top_FC_cluster[i, ], marker_2021_nature[grep(top_FC_cluster$gene[i], marker_2021_nature$MarkerGenes), ])
#     top_FC_cluster_anno <- rbind(top_FC_cluster_anno, anno)
#   }
#   write.csv(top_FC_cluster_anno,row.names = F,file = paste0("Out/1_de_genes_top_FC_PC",PC,"_res",res,"_add_anno_2021_nature_C", cluster[j], ".csv"))
# }


# 8 magic
source("/media/maoni/data/R_functions_ST_lung.R")
Seurat_object_cropped_magic <- magic.Seurat(Seurat_object_cropped)
DefaultAssay(object = Seurat_object_cropped_magic) <- "MAGIC_SCT"
saveRDS(Seurat_object_cropped_magic, file = paste0("Out/", sample, "_cropped_magic_PC", PC, "_res", res, ".RDS"))

# png(paste0("Plot/3_findmarker_Spatial_data_PC",PC,"_res",res,".png"), height = ceiling(length(top10$cluster)/5)*480, width = 5*480)
# plot <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = unique(top10$gene), alpha = c(0.1, 1), ncol = 5, pt.size.factor = 0.9) # + scale_fill_gradientn(colours = viridis::inferno(100))
# plot$layers[[1]]$aes_params=c(plot$layers[[1]]$aes_params, shape=22)
# plot
# dev.off()
# 
# 
# cluster <- names(table(top100$cluster))
# for(i in 1:length(cluster)){
#   features <- top100$gene[top100$cluster==cluster[i]]
#   png(paste0("Plot/3_findmarker_Spatial_data_PC",PC,"_res",res,"_C", cluster[i], ".png"), height = ceiling(length(features)/5)*480, width = 5*480)
#   plot <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = features, alpha = c(0.1, 1), ncol = 5, pt.size.factor = 0.9) # + scale_fill_gradientn(colours = viridis::inferno(100))
#   plot$layers[[1]]$aes_params=c(plot$layers[[1]]$aes_params, shape=22)
#   print(plot)
#   dev.off()
# }

# # 9 highlight genes in papers
# Seurat_object_cropped <- readRDS(paste0("Out/", sample, "_cropped_PC", PC, "_res", res, ".RDS"))
# Seurat_object_cropped_magic <- readRDS(paste0("Out/", sample, "_cropped_magic_PC", PC, "_res", res, ".RDS"))

# SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = c("Dlx1", "Dlx2", "Gad1", "Gad2", "Emx1", "Npy", "Sst", "Lhx6", "Nxph1", "Htr3a", "Prox1", "Cxcl14", "Meis2", "Etv1", "Sp8"), alpha = c(0.1, 1), ncol = 4, pt.size.factor = 0.9) # interneuron
# SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = c("Mdk", "Sstr2", "Neurod1", "Pcp4", "Nrp1", "Rnd2", "Id2", "Ina"), alpha = c(0.1, 1), ncol = 4, pt.size.factor = 0.9) # migrating neuron
# SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = c("Nrn1", "S100b", "Gfap", "Aldh1l1", "Cnp", "Mbp", "Tmem119", "Nkx2-1", "Nkx6-2", "Hmx3"), alpha = c(0.1, 1), ncol = 4, pt.size.factor = 0.9) # excitatory
# SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = c("Aldh1a3", "Foxp1", "Calb1", "Cdh8", "Gbx1", "Lhx6", "Lhx7", "Lhx8", "Shh", "Npy"), alpha = c(0.1, 1), ncol = 4, pt.size.factor = 0.9) # LGE, CP

pdf("Plot/Markers_highlight.pdf", width = 4*7, height = 3*7)
SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = c("Sox2", "Pax6", "Eomes", "Tbr1", "Satb2", "Hes1", "Zic1", "Zic4", "Epha3", "Nrn1", "Lrrtm3"), alpha = c(0.1, 1), ncol = 4, pt.size.factor = 0.9)
dev.off()

pdf("Plot/Markers_highlight.pdf", width = 4*7, height = 3*7)
SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = c("Tmem212", "Wdr63", "Foxj1"), alpha = c(1, 1), ncol = 4, pt.size.factor = 0.9)
dev.off()

features_highlight <- list("Global" = c("Sox2", "Pax6", "Hes5", "Eomes", "Neurog2", "Btg2", "Neurod2", "Neurod6", "Tubb3", "Nrp1", "Neurod1", "Reln", "Lhx5", "Dlx2", "Gad1", "Foxj1"),
                           "Glial_cells" = c("Aldh1l1", "Apoe", "Slc1a3", "Cst3", "Olig1", "Olig2", "Pdgfra", "Aif1", "P2ry12", "C1qb", "Pcna", "Top2a"),
                           "Cortico_Fugal_neuron" = c("Fezf2", "Bcl11b", "Tle4", "Foxp2", "Sox5", "Ldb2", "Crym", "Pcp4", "Nr4a2", "Tshz2"),
                           "Intra_Telencephalic_neuron" = c("Satb2", "Ptn", "Cux1", "Cux2", "Rorb", "Lmo4", "Lpl", "Ldb2", "Pantr1"),
                           "Vascular_cells" = c("Cldn5", "Igfbp7", "Cspg4", "Pdgfrb", "Rgs5", "Lgasl1", "Col3a1", "Lum", "Car2", "Hemgn"),
                           "Ependymal" = c("3300002A11Rik", "Ccdc153", "Tmem212", "1700001C02Rik", "1700012B09Rik", "Dmkn", "Cdc20b", "Foxn4", "Deup1", "Meig1", "Tekt4", "Lrrc23", "Cybrd1", "Dynlrb2", "Meig1", "1700007G11Rik", "Wdr63", "Acta2", "Ctnnb1", "Prom1", "Foxj1"),
                           "Inhibitory" = c("Dlx1", "Dlx2", "Gad1", "Gad2", "Emx1", "Npy", "Sst", "Lhx6", "Nxph1", "Htr3a", "Prox1", "Cxcl14", "Meis2", "Etv1", "Sp8"),
                           "Excitatory" = c("Nrn1", "S100b", "Gfap", "Aldh1l1", "Cnp", "Mbp", "Tmem119", "Nkx2-1", "Nkx6-2", "Hmx3"),
                           "Migrating_neuron" = c("Mdk", "Sstr2", "Neurod1", "Pcp4", "Nrp1", "Rnd2", "Id2", "Ina"),
                           "Other" = c("Aldh1a3", "Foxp1", "Calb1", "Cdh8", "Gbx1", "Lhx6", "Lhx7", "Lhx8", "Shh", "Npy", "Pde10a"))
for(i in 6){
  png(paste0("Plot/3_findmarker_Spatial_data_PC",PC,"_res",res,"_", names(features_highlight)[i], ".png"), height = ceiling(length(unique(unlist(features_highlight[[i]])))/4)*300, width = 4*300)
  plot <- SpatialFeaturePlot(object = Seurat_object_cropped_magic, features = unique(features_highlight[[i]]), alpha = c(0.1, 1), ncol = 4, pt.size.factor = 0.9) # + scale_fill_gradientn(colours = viridis::inferno(100))
  plot$layers[[1]]$aes_params=c(plot$layers[[1]]$aes_params, shape=22)
  print(plot)
  dev.off()
}

################################################################################################
###########################################  revision  #########################################
################################################################################################

############################### Seurat_object_cropped_crop (dapi 1080*1080)
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
removed.genes3 <- gene.names$V2[gene.names$V3!="protein_coding"]
removed.genes <- unique(c(removed.genes1, removed.genes2, removed.genes3))  # 38300

setwd(work_path)
allCounts <- readRDS(paste0("4_zUMIs/", sample, ".dgecounts.rds"))
count_matrix <- allCounts$umicount$inex$all
count_matrix <- as.matrix(count_matrix)
dim(count_matrix) # 35193  9216

raw_stat = read.table(paste0("3_umi_tools/", sample, ".debarcoded_passed_reads_stat.csv"),sep =",", header = TRUE, dec =".", stringsAsFactors = F)
raw_stat_cell = paste(raw_stat$iB, 97-raw_stat$iA, sep="x")
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
dim(count_matrix)

# 1 prepare cropped matrix
stat_cropped = read.table(paste0("4_zUMIs/", sample,".debarcoded_passed_reads_stat_cropped.csv"),sep =",", header = TRUE, dec =".", stringsAsFactors = F)
iA_iB <- paste(stat_cropped$iB, 97-stat_cropped$iA, sep="x")
count_matrix_cropped <- count_matrix[,na.omit(match(iA_iB,colnames(count_matrix)))]
dim(count_matrix_cropped) # 18731

# 2 Seurat object
# See Load10X_Spatial; Read10X_Image
# See satijalab/seurat/issue/3539 4993
assay = "Spatial"
slice = sample
Seurat_object_cropped = CreateSeuratObject(counts = count_matrix_cropped, project = sample, assay = assay, min.cells = 3, min.features = 1) # min.cells min.features
image.dir = "/media/maoni/data/CZP/spatial_transcriptome/seurat_yuhao/core/96*96_blank_img" # "./Img"
image.nam = paste0(sample, "_dapi.png") # "grey_pixel_1080p.png"
coord.nam = "combine_barcode.round2round1_index1_index2.Seurat_3.txt"
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

plot0 <- VlnPlot(Seurat_object_cropped, features = c("nCount_Spatial", "nFeature_Spatial"), ncol = 2, pt.size = 0)
plot0

plot1 <- FeatureScatter(Seurat_object_cropped, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plot2 <- FeatureScatter(Seurat_object_cropped, feature1 = "percent.MT", feature2 = "nFeature_Spatial")
plot1 + plot2

plot3 <- SpatialFeaturePlot(Seurat_object_cropped, features = "nCount_Spatial", pt.size.factor = 0.9, stroke = NA) + theme(legend.position = "top")
plot3$layers[[1]]$aes_params=c(plot3$layers[[1]]$aes_params, shape=22)
plot4 <- SpatialFeaturePlot(Seurat_object_cropped, features = "nFeature_Spatial", pt.size.factor = 0.9, stroke = NA) + theme(legend.position = "top")
plot4$layers[[1]]$aes_params=c(plot4$layers[[1]]$aes_params, shape=22)
plot3 + plot4

pdf("Plot/1_QC_check_vlnplot_dapi_2.pdf", width = 10, height = 4)
print(plot0)
dev.off()
pdf("Plot/1_QC_check_spatial_dapi_2.pdf", width = 10, height = 8)
print(plot1 + plot2 + plot3 + plot4)
dev.off()

save(Seurat_object_cropped, file = paste0("Matrix/", sample, ".cropped_dapi_2.RData"))
load(paste0("Matrix/", sample, ".cropped_dapi_2.RData"))


### crop some regions
# Seurat_object_cropped$crop <- "original"
# Seurat_object_cropped_crop = shiny_st(Seurat_object_cropped_crop, assay = "Spatial", slot = "data", image = Images(Seurat_object_cropped_crop), python_env = "~/anaconda3/envs/python37/bin/python", script = "/media/maoni/data/Software/DBIT_tools/filter_pixel_AI.py")
# table(Seurat_object_cropped_crop$crop)
save(Seurat_object_cropped_crop, file = paste0("Matrix/", sample, ".cropped_dapi_2_cropregions.RData"))
load(paste0("Matrix/", sample, ".cropped_dapi_2_cropregions.RData"))


############################### Seurat_object_cropped_segments (4800*4800)
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
removed.genes3 <- gene.names$V2[gene.names$V3!="protein_coding"]
removed.genes <- unique(c(removed.genes1, removed.genes2, removed.genes3))  # 38300

setwd(work_path)
allCounts <- readRDS(paste0("4_zUMIs/", sample, ".dgecounts.rds"))
count_matrix <- allCounts$umicount$inex$all
count_matrix <- as.matrix(count_matrix)
dim(count_matrix) # 35193  9216

raw_stat = read.table(paste0("3_umi_tools/", sample, ".debarcoded_passed_reads_stat.csv"),sep =",", header = TRUE, dec =".", stringsAsFactors = F)
raw_stat_cell = paste(raw_stat$iB, 97-raw_stat$iA, sep="x")
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
dim(count_matrix)

# 1 prepare cropped matrix
stat_cropped = read.table(paste0("4_zUMIs/", sample,".debarcoded_passed_reads_stat_cropped.csv"),sep =",", header = TRUE, dec =".", stringsAsFactors = F)
iA_iB <- paste(stat_cropped$iB, 97-stat_cropped$iA, sep="x")
count_matrix_cropped <- count_matrix[,na.omit(match(iA_iB,colnames(count_matrix)))]
dim(count_matrix_cropped) # 18731

# 2 Seurat object
# See Load10X_Spatial; Read10X_Image
# See satijalab/seurat/issue/3539 4993
assay = "Spatial"
slice = sample
Seurat_object_cropped = CreateSeuratObject(counts = count_matrix_cropped, project = sample, assay = assay, min.cells = 3, min.features = 1) # min.cells min.features
image.dir = "/media/maoni/data/CZP/spatial_transcriptome/seurat_yuhao/core/96*96_blank_img" # "./Img"
image.nam = paste0(sample, "_segments.png") # "grey_pixel_1080p.png"
coord.nam = "combine_barcode.round2round1_index1_index2.Seurat_4.txt"
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

plot0 <- VlnPlot(Seurat_object_cropped, features = c("nCount_Spatial", "nFeature_Spatial"), ncol = 2, pt.size = 0)
plot0

plot1 <- FeatureScatter(Seurat_object_cropped, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
plot2 <- FeatureScatter(Seurat_object_cropped, feature1 = "percent.MT", feature2 = "nFeature_Spatial")
plot1 + plot2

plot3 <- SpatialFeaturePlot(Seurat_object_cropped, features = "nCount_Spatial", pt.size.factor = 0, stroke = NA) + theme(legend.position = "top")
plot3$layers[[1]]$aes_params=c(plot3$layers[[1]]$aes_params, shape=22)
plot4 <- SpatialFeaturePlot(Seurat_object_cropped, features = "nFeature_Spatial", pt.size.factor = 0.9, stroke = NA) + theme(legend.position = "top")
plot4$layers[[1]]$aes_params=c(plot4$layers[[1]]$aes_params, shape=22)
plot3 + plot4

pdf("Plot/1_QC_check_vlnplot_segments.pdf", width = 10, height = 4)
print(plot0)
dev.off()
pdf("Plot/1_QC_check_spatial_segments.pdf", width = 10, height = 8)
print(plot1 + plot2 + plot3 + plot4)
dev.off()

Seurat_object_cropped_segments <- Seurat_object_cropped
save(Seurat_object_cropped_segments, file = paste0("Matrix/", sample, ".cropped_segments.RData"))
load(paste0("Matrix/", sample, ".cropped_segments.RData"))

############################################################################### start plot ##############################################################
setwd("/media/maoni/data/CZP/spatial_transcriptome/Adult-RNA-23/6_seurat_bigSpotShape_new")
load(paste0("Matrix/", sample, ".cropped_dapi_2_cropregions.RData"))   # Seurat_object_cropped_crop
load(paste0("Matrix/", sample, ".cropped_segments.RData")) # Seurat_object_cropped_segments

# calculate segmentations in each spot
spot_size <- 1080/(95+96)
spot_meta <- Seurat_object_cropped_crop@images$Adult.RNA.23@coordinates
spot_meta$imagerow_up <- round(spot_meta$imagerow+0.5*spot_size)
spot_meta$imagerow_down <- round(spot_meta$imagerow-0.5*spot_size)
spot_meta$imagecol_up <- round(spot_meta$imagecol+0.5*spot_size)
spot_meta$imagecol_down <- round(spot_meta$imagecol-0.5*spot_size)
squidpy_segments_file <- read.table(paste0("/media/maoni/data/CZP/squidpy/", sample, "_segmentation_boundaries.csv"),sep =",", header = FALSE, dec =".", stringsAsFactors = F)
squidpy_segments_table <- table(unlist(squidpy_segments_file))[-1]

### barplot for all classify strategies
for(seq in c(0, 25, 50, 75)){
  
  spot_segments_NOs <- c()
  
  for(i in 1:length(spot_meta[,1])){
    # i <- 100
    spot_segments <- table(unlist(squidpy_segments_file[spot_meta$imagerow_down[i]:spot_meta$imagerow_up[i], spot_meta$imagecol_down[i]:spot_meta$imagecol_up[i]]))[-1]
    spot_segments_percent <- spot_segments/squidpy_segments_table[names(spot_segments)]
    spot_segments_NO <- length(which(spot_segments_percent > seq/100))
    spot_segments_NOs <- c(spot_segments_NOs, spot_segments_NO)
  }
  Seurat_object_cropped_crop$spot_segments_NOs <- spot_segments_NOs
  Seurat_object_cropped_crop[[paste0("spot_segments_NOs_", seq)]] <- spot_segments_NOs
  print(table(Seurat_object_cropped_crop[[paste0("spot_segments_NOs_", seq)]]))
  
}

meta_data <- Seurat_object_cropped_crop@meta.data[, c("spot_segments_NOs_0", "spot_segments_NOs_25", "spot_segments_NOs_50", "spot_segments_NOs_75")]

# Create frequency tables
df0 <- as.data.frame(table(meta_data$spot_segments_NOs_0))
df0$Threshold <- "with_area_largethan0"

df25 <- as.data.frame(table(meta_data$spot_segments_NOs_25))
df25$Threshold <- "with_area_largethan25%"

df50 <- as.data.frame(table(meta_data$spot_segments_NOs_50))
df50$Threshold <- "with_area_largethan50%"

df75 <- as.data.frame(table(meta_data$spot_segments_NOs_75))
df75$Threshold <- "with_area_largethan75%"

combined_df <- bind_rows(df0, df25, df50, df75)
colnames(combined_df) <- c("Segment", "Count", "Threshold")

combined_df$Threshold <- factor(combined_df$Threshold, 
                                levels = c("with_area_largethan0", 
                                           "with_area_largethan25%", 
                                           "with_area_largethan50%", 
                                           "with_area_largethan75%"))

combined_df$Segment <- factor(combined_df$Segment, 
                              levels = sort(unique(as.numeric(as.character(combined_df$Segment)))))


p_list <- list()

p0 <- ggplot(combined_df, aes(x = Threshold, y = Count, fill = Segment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "Segment Distribution by Threshold",
       x = "Threshold",
       y = "Spot Count") +
  theme_classic(base_size = 14) +
  scale_fill_brewer(palette = "Set2") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p0
p_list[[1]] <- p0



### present example for seq=50
for(seq in c(0, 25, 50, 75)){
  
  seq <- 50
  
  spot_segments_NOs <- c()
  for(i in 1:length(spot_meta[,1])){
    # i <- 100
    spot_segments <- table(unlist(squidpy_segments_file[spot_meta$imagerow_down[i]:spot_meta$imagerow_up[i], spot_meta$imagecol_down[i]:spot_meta$imagecol_up[i]]))[-1]
    spot_segments_percent <- spot_segments/squidpy_segments_table[names(spot_segments)]
    spot_segments_NO <- length(which(spot_segments_percent > seq/100))
    spot_segments_NOs <- c(spot_segments_NOs, spot_segments_NO)
  }
  Seurat_object_cropped_crop$spot_segments_NOs <- spot_segments_NOs
  Seurat_object_cropped_crop[[paste0("spot_segments_NOs_", seq)]] <- spot_segments_NOs
  print(table(Seurat_object_cropped_crop[[paste0("spot_segments_NOs_", seq)]]))
  
  
  # Plot spatial spot_segments_NOs
  p1 <- SpatialFeaturePlot(Seurat_object_cropped_crop, features = "spot_segments_NOs", pt.size.factor = 0.5, stroke = NA, alpha = 1) + theme(legend.position = "top") + theme(aspect.ratio = 1) # + scale_fill_viridis_c()
  p1$layers[[1]]$aes_params=c(p1$layers[[1]]$aes_params, shape=22)
  p1
  p_list[[2]] <- p1

  
  # Plot regions
  Seurat_object_cropped_crop_subset <- Seurat_object_cropped_crop[, grep("region0|region2|region3", Seurat_object_cropped_crop$crop)]
  p2 <- SpatialDimPlot(Seurat_object_cropped_crop_subset, group.by = "crop", label = TRUE, label.size = 3, pt.size.factor = 0.5, stroke = NA, crop = FALSE) + theme(legend.position = "top") + theme(aspect.ratio = 1)
  p2$layers[[1]]$aes_params=c(p2$layers[[1]]$aes_params, shape=22)
  p2
  p_list[[3]] <- p2
  
  p0 + p1 + p2
  
  
  ### Seurat_object_cropped_anno(1080*1080)
  
  # Seurat_object_cropped_anno <- readRDS(paste0("Out/", sample, "_cropped_PC", PC, "_res", res, ".RDS"))
  # Seurat_object_cropped_anno$spot_segments_NOs <- spot_segments_NOs
  # Seurat_object_cropped_anno$crop <- Seurat_object_cropped_crop$crop
  # Idents(Seurat_object_cropped_anno) <- Seurat_object_cropped_anno$cell_type
  # 
  # p3 <- VlnPlot(Seurat_object_cropped_anno, features = c("spot_segments_NOs"), ncol = 1, pt.size = 0)
  # p3
  
  color <- "red"
  names(color) <- "region0"
  Seurat_object_cropped_crop_regions1 <- Seurat_object_cropped_crop[, Seurat_object_cropped_crop$crop=="region0"]
  p4 <- SpatialDimPlot(Seurat_object_cropped_crop_regions1, group.by = "crop", label = FALSE, label.size = 3, pt.size.factor = 5, stroke = NA, alpha = 0.5, cols = color) + theme(aspect.ratio = 1)
  p4$layers[[1]]$aes_params=c(p4$layers[[1]]$aes_params, shape=22)
  p4
  p_list[[4]] <- p4
  
  color <- "red"
  names(color) <- "region2"
  Seurat_object_cropped_crop_regions2 <- Seurat_object_cropped_crop[, Seurat_object_cropped_crop$crop=="region2"]
  p5 <- SpatialDimPlot(Seurat_object_cropped_crop_regions2, group.by = "crop", label = FALSE, label.size = 3, pt.size.factor = 5, stroke = NA, alpha = 0.5, cols = color) + theme(aspect.ratio = 1)
  p5$layers[[1]]$aes_params=c(p5$layers[[1]]$aes_params, shape=22)
  p5
  p_list[[5]] <- p5
  
  color <- "red"
  names(color) <- "region3"
  Seurat_object_cropped_crop_regions3 <- Seurat_object_cropped_crop[, Seurat_object_cropped_crop$crop=="region3"]
  p6 <- SpatialDimPlot(Seurat_object_cropped_crop_regions3, group.by = "crop", label = FALSE, label.size = 3, pt.size.factor = 5, stroke = NA, alpha = 0.5, cols = color) + theme(aspect.ratio = 1)
  p6$layers[[1]]$aes_params=c(p6$layers[[1]]$aes_params, shape=22)
  p6
  p_list[[6]] <- p6
  
  
  Seurat_object_cropped_segments$crop <- Seurat_object_cropped_crop$crop
  Seurat_object_cropped_segments$spot_segments_NOs <- Seurat_object_cropped_crop$spot_segments_NOs
  
  color <- "red"
  names(color) <- "region0"
  Seurat_object_cropped_segments_regions1 <- Seurat_object_cropped_segments[, Seurat_object_cropped_segments$crop=="region0"]
  p7 <- SpatialDimPlot(Seurat_object_cropped_segments_regions1, group.by = "crop", label = FALSE, label.size = 3, pt.size.factor = 5, stroke = NA, alpha = 0.5, cols = color) + theme(aspect.ratio = 1)
  p7$layers[[1]]$aes_params=c(p7$layers[[1]]$aes_params, shape=22)
  p7
  p_list[[7]] <- p7
  
  color <- "red"
  names(color) <- "region2"
  Seurat_object_cropped_segments_regions2 <- Seurat_object_cropped_segments[, Seurat_object_cropped_segments$crop=="region2"]
  p8 <- SpatialDimPlot(Seurat_object_cropped_segments_regions2, group.by = "crop", label = FALSE, label.size = 3, pt.size.factor = 5, stroke = NA, alpha = 0.5, cols = color) + theme(aspect.ratio = 1)
  p8$layers[[1]]$aes_params=c(p8$layers[[1]]$aes_params, shape=22)
  p8
  p_list[[8]] <- p8
  
  color <- "red"
  names(color) <- "region3"
  Seurat_object_cropped_segments_regions3 <- Seurat_object_cropped_segments[, Seurat_object_cropped_segments$crop=="region3"]
  p9 <- SpatialDimPlot(Seurat_object_cropped_segments_regions3, group.by = "crop", label = FALSE, label.size = 3, pt.size.factor = 5, stroke = NA, alpha = 0.5, cols = color) + theme(aspect.ratio = 1)
  p9$layers[[1]]$aes_params=c(p9$layers[[1]]$aes_params, shape=22)
  p9
  p_list[[9]] <- p9
  
  # 使用 ggplot 绘制 segments_NOs 直方图
  data <- data.frame(table(Seurat_object_cropped_segments_regions1$spot_segments_NOs))
  colnames(data) <- c("Segment_number", "Spot_number")
  p_list[[10]] <- ggplot(data, aes(x = Segment_number, y = Spot_number)) +
    geom_bar(stat = "identity", fill = "#A4A4A4", color = "#A4A4A4", width = 0.6) +
    labs(title = paste0("region0"), x = "Segment_number", y = "Spot_number") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +   # 旋转 X 轴标签为 45 度 
    theme(aspect.ratio = 1) + 
    scale_y_continuous(expand = c(0, 0))  # ✅ y轴从0开始
  
  > data
  Segment_number Spot_number
  1              0          68
  2              1          29
  3              2           3
    
  data <- data.frame(table(Seurat_object_cropped_segments_regions2$spot_segments_NOs))
  colnames(data) <- c("Segment_number", "Spot_number")
  p_list[[11]] <- ggplot(data, aes(x = Segment_number, y = Spot_number)) +
    geom_bar(stat = "identity", fill = "#A4A4A4", color = "#A4A4A4", width = 0.6) +
    labs(title = paste0("region2"), x = "Segment_number", y = "Spot_number") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +   # 旋转 X 轴标签为 45 度 
    theme(aspect.ratio = 1) + 
    scale_y_continuous(expand = c(0, 0))  # ✅ y轴从0开始
  
  > data
  Segment_number Spot_number
  1              0          59
  2              1          29
  3              2          11
  4              3           1
  
  data <- data.frame(table(Seurat_object_cropped_segments_regions3$spot_segments_NOs))
  colnames(data) <- c("Segment_number", "Spot_number")
  p_list[[12]] <- ggplot(data, aes(x = Segment_number, y = Spot_number)) +
    geom_bar(stat = "identity", fill = "#A4A4A4", color = "#A4A4A4", width = 0.6) +
    labs(title = paste0("region3"), x = "Segment_number", y = "Spot_number") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +   # 旋转 X 轴标签为 45 度 
    theme(aspect.ratio = 1) + 
    scale_y_continuous(expand = c(0, 0))  # ✅ y轴从0开始
  
  > data  
  Segment_number Spot_number
  1              0          80
  2              1          19
  3              2           1
  
  p_list[[10]] + p_list[[11]] + p_list[[12]]
  
  pdf(paste0("/media/maoni/data/CZP/squidpy/", sample, "_spot_segments_NOs_seq_", seq, "_and_regions_highlight_alpha_0.5.pdf"), height = 20, width = 15)
  p <- ggarrange(plotlist = p_list, ncol = 3, nrow = 4)
  print(p)
  dev.off()
  
  
  save(combined_df, Seurat_object_cropped_segments, Seurat_object_cropped_crop, Seurat_object_cropped_segments_regions1, Seurat_object_cropped_segments_regions2, Seurat_object_cropped_segments_regions3, file = paste0("Matrix/", sample, ".cropregions_segmentation_data_for_plot.RData"))
  
}



