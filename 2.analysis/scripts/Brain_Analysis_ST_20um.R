source("/home/yiyelinfeng/scripts/Rscripts/scRNASeqSeuratAnalysisFunction.R")
source("/home/yiyelinfeng/scripts/Rscripts/R.analysis.Functions.R")
source("/home/yiyelinfeng/scripts/Rscripts/lung_project/IPF/spatial-lung-fibrosis/scripts/custom_colors.R")
outPath = "/media/yiyelinfeng/data/Projects/brain_project/Spatial_analysis"
dir.create(outPath)
setwd(outPath)
library(Seurat)
library(qs)
library(ggplot2)
library(cowplot)
library(patchwork)
library(limma)
library(dplyr)
#library(Nebulosa)
require(stringr)
library(SeuratWrappers)
options(stringsAsFactors = FALSE)

library(scCustomize)
selected_colors <- DiscretePalette_scCustomize(num_colors = 40, palette = "ditto_seq")
selected_colors <- selected_colors[-c(8, 16, 24)]
names(selected_colors) <- 0:36

library(pals)
mypal <- kelly()[-1]

sample.name <- "brain_mouse_ST"

outPathT <- paste0(outPath, "/", sample.name)
if(!file.exists(outPathT)){
	dir.create(outPathT)
}
setwd(outPathT)

cluster_colors <- c("E14_5" = "#34D916", "E16_5" = "#00D4E6", "E18_5" = "#1E90FF", "GZ" = "#66C2A5", "CP" = "#FC8D62",
"Cortex" = "#59B375", "MGE" = "#9E7BFF", "LGE" = "#89288F", "LCS" = "#C1FF73", "CLA" = "#33FF00", "LS" = "#CC79A7", "Epd" = "#94FFB5", "PIR" = "#B3823E", "others" = "#FEE52C", "Others" = "#dde2e6",
"Cortex_AP" = "#03FFF4", "Cortex_IPC_MigN" = "#036DF4", "Cortex_IPC" = "#0BD3B1", "Cortex_MigN" = "#62CFE8",
"Cortex_Layer_6b" = "#99FFBC", "Cortex_Layer_6b_CThPN" = "#009E73", "Cortex_CThPN_SCPN" = "#34A047", "Cortex_CThPN" = "#7EC136", "Cortex_SCPN_UL_CPN" = "#007756", "Cortex_SCPN" = "#00991F", "Cortex_UL_CPN" = "#01545a", "Cortex_CR" = "#F2F318",
"MGE_LGE_LS_AP" = "#9E7BFF", "MGE_AP" = "#8D73DF", "LGE_AP" = "#9E7BFF", "MGE_InN" = "#BE50FF", "LGE_IPC" = "#D85FF7", "LGE_InN" = "#AA0DFE", "CPU_pre_MSN" = "#F07F92", "CPU_MSN" = "#D38B5C",
"LCS_IMN" = "#C1FF73", "CLA_ExN" = "#33FF00", "PIR_ExN" = "#B5EFB5", "PAL" = "#96C885", "LS_AP" = "#7959FA", "LS_IPC" = "#E0AFCA", "LS" = "#CC79A7",
"MS" = "#07F4B3", "DB" = "#2271A9", "AVPV" = "#FF9B4D", "POA" = "#B14380", "IIIV" = "#F4CE70",
"Endo" = "#994567", "Stromal" = "#886C00", "unknown" = "#dde2e6")

ClusterNames <- c("Cortex_AP", "Cortex_IPC_MigN", "Cortex_IPC", "Cortex_MigN", "Cortex_Layer_6b", "Cortex_Layer_6b_CThPN", "Cortex_CThPN_SCPN", "Cortex_CThPN", "Cortex_SCPN_UL_CPN", "Cortex_SCPN", "Cortex_UL_CPN", "Cortex_CR", "MGE_LGE_LS_AP", "MGE_AP", "LGE_AP", "MGE_InN", "LGE_IPC", "LGE_InN", "CPU_pre_MSN", "CPU_MSN", "LCS_IMN", "PIR_ExN", "CLA_ExN", "PAL", "LS_AP", "LS_IPC", "LS", "MS", "DB", "AVPV", "POA", "IIIV", "Endo", "Stromal")
ClusterID <- paste0("C", 1:length(ClusterNames))
ClusterIDColors <- cluster_colors[match(ClusterNames, names(cluster_colors))]
names(ClusterIDColors) <- ClusterID
names(ClusterNames) <- ClusterID
#!-----------------------------------------------------------------------------------------------------------
gene.names <- read.table("../../data/Mus_musculus.GRCm38.102.gene.id_name.tsv", sep = "\t", header = T)

# remove unnecessary pseudogene, rRNA, snoRNA, tRNA, ribozyme
removed.genes <- gene.names$gene.name[-grep("protein_coding", gene.names$gene.bio_type)]

min.features <- 200
min.cells <- 3

sample_names <- list.files("../../data/rds/")
sample_names <- sample_names[grep(".dgecounts.rds", sample_names)]
sample_names <- sub(pattern = '.dgecounts.rds', replacement = '', x = sample_names)

sample_names <- sample_names[c(3:5, 2, 1)]

dup_rates <- c()
seurat_list <- list()
count = 1

dir.create(paste0(outPathT, "/output"))

for(sample_name in sample_names){
	print(paste0("now start deal sample with ", sample_name, "..."))
	barcodeA <- read.table(paste0("../../data/barcode/barcodesA.txt"))
	barcodeB <- read.table(paste0("../../data/barcode/barcodesB.txt"))

	ST_barcodes <- c()
	index <- c()
	for(i in 1:length(barcodeB$V1)){
		for(j in 1:length(barcodeA$V1)){
			ST_barcodes <- c(ST_barcodes, paste0(barcodeB$V1[i], barcodeA$V1[j]))
			index <- c(index, paste0(i, 'x', j))
		}
	}

	ST_barcodes <- data.frame(ST_barcodes, index)

	allCounts <- readRDS(paste0("../../data/rds/", sample_name, ".dgecounts.rds"))

	count_matrix <- allCounts$umicount$inex$all
	count_matrixT <- allCounts$readcount$inex$all

	dup_rate <- 1 - sum(count_matrix) / sum(count_matrixT)
	dup_rates <- c(dup_rates, dup_rate)

	temp <- intersect(ST_barcodes$ST_barcodes, colnames(count_matrix))
	ST_barcode <- ST_barcodes[match(temp, ST_barcodes$ST_barcodes),]
	ST_barcode <- ST_barcode[match(colnames(count_matrix), ST_barcode$ST_barcodes),]

	colnames(count_matrix) <- ST_barcode$index
	tmp <- as.data.frame(strsplit(colnames(count_matrix), split = 'x'))
	colnames(count_matrix) <- paste0(tmp[2,], 'x', tmp[1,])
	
	print(dim(count_matrix))
	rownames(count_matrix) <- gene.names[,2][match(rownames(count_matrix), gene.names[,1])]
	print(dim(count_matrix))

	UMI <- sort(rowSums(count_matrix), decreasing = T)
	Reads <- sort(rowSums(count_matrixT), decreasing = T)
	all_counts <- as.data.frame(cbind(UMI, Reads))
	all_counts$UMI_rate <- all_counts$UMI / sum(all_counts$UMI)
	all_counts$read_rate <- all_counts$Reads / sum(all_counts$Reads)
	write.csv(all_counts[,c(1,3, 2, 4)], paste0(sample_name, "_mun_UMI_reads_per_genes.csv"), quote = F)
	rm(count_matrixT)

	ID <- as.character(rownames(count_matrix))
	d <- duplicated(ID)
	ID <- factor(ID, levels = unique(ID))
	count_matrix <- rowsum(as.matrix(count_matrix), ID, reorder = FALSE, na.rm = TRUE)
	print(dim(count_matrix))

	removed.gene <- intersect(rownames(count_matrix), removed.genes)
	count_matrix <- count_matrix[-match(removed.gene, rownames(count_matrix)),]
	print(dim(count_matrix))

	temp_object <- CreateSeuratObject(counts = count_matrix, project = sample_name, min.cells = ifelse(floor(ncol(count_matrix) * 0.001) > min.cells, floor(ncol(count_matrix) * 0.001), min.cells), min.features = min.features)
	temp_object[["percent.mt"]] <- PercentageFeatureSet(object = temp_object, pattern = "^mt")
	temp_object[["percent.ribo"]] <- PercentageFeatureSet(object = temp_object, pattern = "^Rp[sl]")
	
	temp_object <- create_image_object(temp_object, png_path = paste0("../../data/fix/", sample_name, "_fix.png"), new_assay = sample_name)
	tmp <- as.data.frame(strsplit(colnames(temp_object), split = 'x'))
	temp_object$row <- as.numeric(tmp[2,])
	temp_object$col <- as.numeric(tmp[1,])
	seurat_list[count] <- temp_object
	count <- count + 1
}
names(seurat_list) <- sample_names

source("~/scripts/seurat2scanpy/shiny_st.R")
options(browser = "/usr/bin/firefox")

i = 1
seurat_list[[i]] <- shiny_st(seurat = seurat_list[[i]], isVisium = F, assay = "RNA", image = names(seurat_list)[i])

for(i in 1:length(seurat_list)){
	seurat_list[[i]] <- subset(seurat_list[[i]], subset = seurat_clusters != "removed_cells")
}

seurat_list[[1]] <- subset(seurat_list[[1]], subset = nCount_RNA < 25000 & nCount_RNA > 500)
seurat_list[[2]] <- subset(seurat_list[[2]], subset = nCount_RNA < 25000 & nCount_RNA > 500)
seurat_list[[3]] <- subset(seurat_list[[3]], subset = nCount_RNA < 30000 & nCount_RNA > 500)

saveRDS(seurat_list, "output/seurat_list.rds")


cc.genes <- cc.genes.updated.2019
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
gene.mouse <- read.csv("/media/yiyelinfeng/data/Repository/baseFiles/ortholog_mouse.csv", header = T)
#!-------------------------------------------------------------------------------------------------------------------------------------
vf.nfeatures = 3000
HVGs <- SelectIntegrationFeatures(seurat_list, nfeatures = 3000)

library(reticulate)
use_python("/home/yiyelinfeng/softwares/miniconda3/envs/R4/bin")
SPARC <- import("SPARC", convert = FALSE)
assay_name <- paste0("spARC_SCT")

for(i in 1:length(seurat_list)){
	DefaultAssay(seurat_list[[i]]) <- "RNA"
	seurat_list[[i]] = SCTransform(seurat_list[[i]], vst.flavor = "v2", method = "glmGamPoi", verbose = TRUE)
	s.gene <- intersect(rownames(seurat_list[[i]]), s.genes)
	g2m.gene <- intersect(rownames(seurat_list[[i]]), g2m.genes)
	s.gene <- intersect(rownames(seurat_list[[i]]), gene.mouse$Mouse.gene.name[match(s.genes, gene.mouse$Gene.name)])
	g2m.gene <- intersect(rownames(seurat_list[[i]]), gene.mouse$Mouse.gene.name[match(g2m.genes, gene.mouse$Gene.name)])
	seurat_list[[i]] <- CellCycleScoring(seurat_list[[i]], s.features = s.gene, g2m.features = g2m.gene, set.ident = FALSE, nbin = 12)
	seurat_list[[i]]$CC.Difference <- seurat_list[[i]]$S.Score - seurat_list[[i]]$G2M.Score
	DefaultAssay(seurat_list[[i]]) <- "RNA"
	seurat_list[[i]] = SCTransform(seurat_list[[i]], vst.flavor = "v2", method = "glmGamPoi", vars.to.regress = c("S.Score", "G2M.Score"), verbose = TRUE)
	spatial_X <- Seurat::GetTissueCoordinates(seurat_list[[i]])[,seq_len(2)]
	data_sparc <- SPARC$spARC()$fit_transform(expression_X = t(x = as.matrix(LayerData(seurat_list[[i]], assay = "SCT", layer = "data"))), spatial_X = spatial_X)
	data_sparc <- t(x = as.matrix(x = data_sparc))
	colnames(data_sparc) <- colnames(seurat_list[[i]])
    rownames(data_sparc) <- rownames(seurat_list[[i]])
    seurat_list[[i]][[assay_name]] <- CreateAssay5Object(data = data_sparc)
	seurat_list[[i]] <- RunPCA(seurat_list[[i]], assay = "SCT", reduction.name = "pca", npcs = 30)
	seurat_list[[i]] <- RunUMAP(seurat_list[[i]], reduction = "pca", dims = 1:30, min.dist = 0.3, reduction.name = "SCT_umap", return.model = T)
	seurat_list[[i]] <- FindNeighbors(seurat_list[[i]], reduction = "pca", dims = 1:30, annoy.metric = "cosine", graph.name = c('SCT_nn', 'SCT_snn'))
}
saveRDS(seurat_list, "output/seurat_list_SCT.rds")
#!---------------------------------------------------------------------------------------------------------------
# deal E14_5 spatial RNA
tmp_object <- seurat_list$E14_5
tmp_object <- FindClusters(tmp_object, resolution = seq(0.1, 2, 0.1), graph.name = "SCT_snn")
wrap_plots(map(seq(0.6, 2, 0.1), function(x) DimPlot(tmp_object, reduction = "SCT_umap", group.by = paste0("SCT_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)
Seurat::SpatialPlot(tmp_object, shape = 22, cols = selected_colors, group.by = paste0("SCT_snn_res.", seq(1.1, 2, 0.1)), ncol = 5) & NoLegend()

tmp_object <- FindClusters(tmp_object, cluster.name = "SCT_cluster", resolution = 1.6, graph.name = "SCT_snn")
tmp_object$SCT_cluster[which(tmp_object$SCT_snn_res.0.6 == 4)] <- 2
Idents(tmp_object) <- "SCT_cluster"
selected_cells <- CellsByIdentities(tmp_object)
Seurat::SpatialPlot(tmp_object, cells.highlight = selected_cells[setdiff(names(selected_cells), "NA")], cols.highlight = c("#FFFF00", alpha("lightgrey", 0)), crop = FALSE, shape = 22, facet.highlight = T, combine = T, ncol = 7, alpha = NULL) & NoLegend()

all_markers <- FindAllMarkers(tmp_object, only.pos = T, min.pct = 0.1)
E14_de_markers <- all_markers[which(all_markers$avg_log2FC > 1 & all_markers$p_val_adj < 0.05),]
Seurat::SpatialPlot(tmp_object, shape = 22, crop = F, ncol = 5, features = E14_de_markers$gene[which(E14_de_markers$cluster == 0)][1:15]) & scale_fill_gradientn(colours = paletteContinuous(set = "sambaNight"))

tmp_object$cell_identity <- as.character(tmp_object$SCT_cluster)
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(3))] <- "Cortex_AP" # "Eomes", "Pax6", "Celsr1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(7))] <- "Cortex_IPC" # "Slc17a6", "Unc5d", "Plcb1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(4))] <- "Cortex_MigN" # "Cntn2", "Satb2", "Nos1", "Robo2"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(5, 19))] <- "Cortex_CThPN_SCPN" # "Fezf2", "Dync1i1", "Nr4a3", "Nin", "Fat4", "Rgs6"

tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(9))] <- "MGE_AP" # "Vit", "Mt3", "Acss1", "Ildr2"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(18, 24))] <- "MGE_InN" # "Lhx6", "Kcnc2", "Lhx8", "Gbx1", "Nkx2-1"

tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(21))] <- "LGE_AP" # "Vit", "Mt3", "Acss1", "Ildr2"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(10))] <- "LGE_IPC" # "Esrrg", "Dlx6", "Slc18a2", "Meis1", "Cemip2"

tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(8))] <- "CPU_pre_MSN" # "Ebf1", "Rarb", "Foxp1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(13))] <- "CPU_MSN"  # "Drd1", "Kcnh1", "Rxrg", "Trpc7", "Myo3b"

tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(17))] <- "LCS_IMN" # "Tfap2d"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(14))] <- "CLA_ExN" # "Tafa1", "Ppp1r14c", "Unc13c", "Lmo3", "Nr4a2"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(12))] <- "PIR_ExN" # "Tafa1", "Ppp1r14c", "Unc13c", "Lmo3", "Nr4a2"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(0, 15))] <- "PAL" # "Adarb2", 15 "Cacna2d3", "Cntnap5b", "Asic2"

tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(25))] <- "LS_AP" # "Vit", "Mt3", "Acss1", "Ildr2"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(20))] <- "LS_IPC" # "St18", "Fgd5"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(11))] <- "LS" # "Trpc4", "Sntg1", "Fgd5"

tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(16))] <- "MS" # "Onecut1", "Onecut2", "Onecut3", "Ntng1", "Prox1", "Shisa6"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(1))] <- "AVPV" # "Snhg11", "Samd5", "Ahi1", "Gbx1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(6))] <- "POA" # "Peg3", "Pou2f2", "Nova1", "Hmx2", "Hmx3"

tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(22))] <- "IIIV" # "Lrp2", "Slit2", "Peg3"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(23))] <- "IIIV" # "Pax2", "Cp", "Clybl", "Wnt7b", "Tenm3"

tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(2, 26))] <- "Stromal" # "Col4a1", "Col4a2", "Cldn5"

tmp_object$cell_identity <- factor(tmp_object$cell_identity, levels = c("Cortex_AP", "Cortex_IPC", "Cortex_MigN", "Cortex_CThPN_SCPN", "MGE_AP", "MGE_InN", "LGE_AP", "LGE_IPC", "CPU_pre_MSN", "CPU_MSN", "LCS_IMN", "CLA_ExN", "PIR_ExN", "PAL", "LS_AP", "LS_IPC", "LS", "MS", "AVPV", "POA", "IIIV", "Stromal"))

Seurat::SpatialPlot(tmp_object, group.by = "cell_identity", shape = 22, cols = cluster_colors, crop = F)
#!------------------------------
Idents(tmp_object) <- "cell_identity"
all_markers <- FindAllMarkers(tmp_object, min.pct = 0.1, only.pos = T)
de_markers <- all_markers[which(all_markers$avg_log2FC > 1 & all_markers$p_val_adj < 0.01 & all_markers$pct.1 > 0.3),]
write.csv(de_markers, "brain_20um_analysis/E14_5_ST_cell_identity_DE_markers.csv", quote = F)

key_markers <- c("Satb2", "Top2a", "Mki67", "Neurog2", "Eomes", "Hmga2", "Pax6", "Dct", "Adamts19", "Emid1", "Slc17a6", "Nrn1", "Neurod1", "Adamtsl3", "Cdh12", "Neurod2", "Neurod6", "Cntn2", "Abcc8", "Nfe2l3", "Trpc3", "Rgs6", "Dync1i1", "Opcml", "Nr4a3", "Fezf2", "Reln", "Dpp10", "Trp73",  "Zhx2","Vit", "Ednrb", "Ildr2", "Sp9", "Dlx1", "Dlx6", "Stk33", "Nkx2-1", "Lhx6", "Lhx8", "Rai2", "Esrrg", "Rffl", "Slc18a2", "Dlx2", "Gucy1a1", "Zfp503", "Isl1", "Rarb", "Oprm1", "Frmd6", "Ikzf1", "Grid1", "Gad2", "Gad1", "Tfap2d", "Nr2f1", "Tafa2", "Tafa1", "Lhfpl3", "Unc13c", "Fstl4", "Sphkap", "Plppr5", "Nr4a2", "Galnt14", "Zic4", "Otx2", "Fgf17", "Lrp2", "Bcan", "Slc6a11", "Kdr", "Flt1", "Ptprb", "Ahnak", "Col5a1", "Dcn", "Bmp6", "Ranbp3l", "Itih5")

DotPlot(tmp_object, features = key_markers, col.min = 0.1, dot.min = 0.1) + RotatedAxis() + labs(y = "cell identity", x = "Features") + scale_color_gradientn(colours = paletteContinuous(set = "whitePurple"))
ggsave("brain_20um_analysis/E14_5_cell_identity_key_markers.pdf")
#!------------------------------
tmp_object$Cluster <- names(ClusterNames)[match(as.character(tmp_object$cell_identity), ClusterNames)]
tmp_object$Cluster <- factor(tmp_object$Cluster, levels = intersect(ClusterID, unique(tmp_object$Cluster)))
Seurat::SpatialPlot(tmp_object, group.by = "Cluster", shape = 22, cols = ClusterIDColors)

Idents(tmp_object) <- "cell_identity"
seurat_list$E14_5 <- tmp_object
#!---------------------------------------------------------------------------------------------------------------
# deal E16_5 spatial RNA
tmp_object <- seurat_list$E16_5
tmp_object <- FindClusters(tmp_object, resolution = seq(0.1, 2, 0.1), graph.name = "SCT_snn")
wrap_plots(map(seq(0.6, 2, 0.1), function(x) DimPlot(tmp_object, reduction = "SCT_umap", group.by = paste0("SCT_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)
Seurat::SpatialPlot(tmp_object, shape = 22, cols = selected_colors, group.by = paste0("SCT_snn_res.", seq(1.1, 2, 0.1)), ncol = 5) & NoLegend()

tmp_object <- FindClusters(tmp_object, cluster.name = "SCT_cluster", resolution = 1.8, graph.name = "SCT_snn")
tmp_object <- FindSubCluster(tmp_object, graph.name = "SCT_snn", cluster = 6, subcluster.name = "SCT_sub_cluster", resolution = 0.2)
tmp_object$SCT_cluster <- as.character(tmp_object$SCT_cluster)
tmp_object$SCT_cluster[which(tmp_object$SCT_snn_res.0.6 == 14)] <- 27
tmp_object$SCT_cluster[which(tmp_object$SCT_sub_cluster == "6_2")] <- 29
tmp_object$SCT_cluster <- factor(tmp_object$SCT_cluster, levels = 0:29)
Idents(tmp_object) <- "SCT_cluster"
selected_cells <- CellsByIdentities(tmp_object)
Seurat::SpatialPlot(tmp_object, cells.highlight = selected_cells[setdiff(names(selected_cells), "NA")], cols.highlight = c("#FFFF00", alpha("lightgrey", 0)), crop = FALSE, shape = 22, facet.highlight = T, combine = T, ncol = 8, alpha = NULL) & NoLegend()

all_markers <- FindAllMarkers(tmp_object, only.pos = T, min.pct = 0.1)
E16_de_markers <- all_markers[which(all_markers$avg_log2FC > 1 & all_markers$p_val_adj < 0.05),]
Seurat::SpatialPlot(tmp_object, shape = 22, crop = F, ncol = 5, features = E16_de_markers$gene[which(E16_de_markers$cluster == 0)][1:15]) & scale_fill_gradientn(colours = paletteContinuous(set = "sambaNight"))

tmp_object$cell_identity <- as.character(tmp_object$SCT_cluster)

tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(14))] <- "Cortex_AP" # "Eomes", "Celsr2", "Pax6"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(4))] <- "Cortex_IPC" # "Slc17a6", "Unc5d", "Sema3c", "Cux2", "Gramd1b", "Mllt3"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(7, 26))] <- "Cortex_MigN" # "Satb2", "Cntn2", "Osbpl6", "Zfpm2"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(3))] <- "Cortex_Layer_6b_CThPN" # "Hs3st4", "Pdzd2", "Zfpm2", "Fezf2"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(1))] <- "Cortex_SCPN_UL_CPN" # "Ccbe1", "Kcnq5", "Rorb", "Mef2c", "Tmtc1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(23))] <- "Cortex_CR" # Reln

tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(17))] <- "MGE_LGE_LS_AP" # "Mt3", "Acsbg1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(13))] <- "MGE_InN" # "Lhx8", "Lhx6", "Pde1a", "Kcnc2"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(5))] <- "LGE_IPC" # "Dlx6", "Mpped2", "Dlx1", "Gad2"

tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(19))] <- "CPU_pre_MSN" # "Ikzf1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(8))] <- "CPU_MSN" # "Rxrg", "Myo3b", "Rarb", "Foxp1"

tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(10))] <- "CLA_ExN" # "Tafa1", "Ppp1r14c", "Lmo3", "Fstl4"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(24))] <- "CLA_ExN" # "Nr4a2", "Sv2b", "Galnt14"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(22))] <- "PIR_ExN" # "Tafa1", "Ppp1r14c", "Lmo3", "Fstl4"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(25))] <- "LCS_IMN" # Tfap2d

tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(2))] <- "PAL" # "Adarb2", "Tmeff2", "Samd5"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(12))] <- "PAL" # "Cacna2d3", "Fat3"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(9))] <- "PAL" # "Cacna2d3", "Cntnap5b", "Pde8b", "Mgat4c"

tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(27))] <- "LS_IPC"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(6))] <- "LS" # "Trpc4", "Gria1", "Ripor2", "Sntg1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(15))] <- "LS" # "Fgd5", "Ano1", "Gfra1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(29))] <- "MS"

tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(18))] <- "DB" # "Grik1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(21))] <- "DB" # "Casz1", "Gbx1", "Pde3a", "Nefm"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(11))] <- "AVPV"  # "Prlr", "Samd3", "Zdbf2", "Car10"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(16))] <- "POA" # "Hmx2", "Hmx3", "Skor2", "Esr1", "Hap1", "Scg2"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(28))] <- "IIIV" # "Traf1", "Ccdc60"

tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(0, 20))] <- "Stromal" # "Col4a1", "Col4a2", "Col3a1", "Myh9", "Sparc"

tmp_object$cell_identity <- factor(tmp_object$cell_identity, levels = c("Cortex_AP", "Cortex_IPC", "Cortex_MigN", "Cortex_Layer_6b_CThPN", "Cortex_SCPN_UL_CPN", "Cortex_CR", "MGE_LGE_LS_AP", "MGE_InN", "LGE_IPC", "CPU_pre_MSN", "CPU_MSN", "CLA_ExN", "PIR_ExN", "LCS_IMN", "PAL", "LS_IPC", "LS", "MS", "DB", "AVPV", "POA", "IIIV", "Stromal"))
Seurat::SpatialPlot(tmp_object, group.by = "cell_identity", shape = 22, cols = cluster_colors, crop = F)
#!------------------------------
Idents(tmp_object) <- "cell_identity"
all_markers <- FindAllMarkers(tmp_object, min.pct = 0.1, only.pos = T)
de_markers <- all_markers[which(all_markers$avg_log2FC > 1 & all_markers$p_val_adj < 0.01 & all_markers$pct.1 > 0.3),]
write.csv(de_markers, "brain_20um_analysis/E16_5_ST_cell_identity_DE_markers.csv", quote = F)

key_markers <- c("Satb2", "Top2a", "Mki67", "Eomes", "Pax6", "Veph1", "Ccdc80", "Aldoc", "Neurog2", "Hmga2", "Dct", "Pter", "Hopx", "Ltbp1", "Slc17a6", "Nrn1", "Adamts2","Calb2", "Neurod1", "Cemip", "Nr3c2", "Abca8a", "Tmem132d", "Pappa2", "Nfe2l3", "Klhl1", "Npy", "Npr3",  "Nr4a3","Igfbp3", "Nwd2", "Tafa1", "Tcerg1l", "Rorb", "Ntf3", "Trmt9b", "Chst8", "Usp43", "Tafa2", "Bcl6", "Inhba", "Frem2", "Mkx", "Mafb", "Lhx6", "Npas1", "Nxph1", "Lhcgr", "Pdzph1", "Vcam1", "Sp9", "Dlx6", "Dlx1", "Dlx2", "Crb1", "Six3", "Drd2", "Dchs2", "Zfp503", "Ikzf1", "Cntnap3", "Dclk3", "Oprm1", "Ano3", "Rxrg", "Myo3b", "Kdr", "Zfp366", "Fli1", "Cldn5", "Tmem72", "Folr1", "Clic6", "Gmnc", "Tfap2d", "Lypd6", "Nr4a2", "Prss12", "Runx2", "Daw1", "Hepacam", "Mlc1", "Glis3", "Zic4", "Zic1", "Zic2", "Epas1", "Igfbp7", "Ndnf", "Col3a1", "Col1a1")

DotPlot(tmp_object, features = key_markers, col.min = 0.1, dot.min = 0.1) + RotatedAxis() + labs(y = "cell identity", x = "Features") + scale_color_gradientn(colours = paletteContinuous(set = "whitePurple"))
ggsave("brain_20um_analysis/E16_5_cell_identity_key_markers.pdf")
#!----------------------------
tmp_object$Cluster <- ClusterID[match(as.character(tmp_object$cell_identity), ClusterNames)]
tmp_object$Cluster <- factor(tmp_object$Cluster, levels = intersect(ClusterID, unique(tmp_object$Cluster)))
Seurat::SpatialPlot(tmp_object, group.by = "Cluster", shape = 22, cols = ClusterIDColors, crop = F)

Idents(tmp_object) <- "cell_identity"
seurat_list$E16_5 <- tmp_object
#!---------------------------------------------------------------------------------------------------------------
# deal E18_5 spatial RNA
tmp_object <- seurat_list$E18_5
tmp_object <- FindClusters(tmp_object, resolution = seq(0.1, 2, 0.1), graph.name = "SCT_snn")
wrap_plots(map(seq(0.6, 2, 0.1), function(x) DimPlot(tmp_object, reduction = "SCT_umap", group.by = paste0("SCT_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)
Seurat::SpatialPlot(tmp_object, shape = 22, cols = selected_colors, group.by = paste0("SCT_snn_res.", seq(1.1, 2, 0.1)), ncol = 5) & NoLegend()

tmp_object <- FindClusters(tmp_object, cluster.name = "SCT_cluster", resolution = 1.7, graph.name = "SCT_snn")
tmp_object <- FindSubCluster(tmp_object, graph.name = "SCT_snn", cluster = 3, subcluster.name = "SCT_sub_cluster", resolution = 0.2)
tmp_object <- FindSubCluster(tmp_object, graph.name = "SCT_snn", cluster = 4, subcluster.name = "SCT_sub_cluster1", resolution = 0.2)
tmp_object$SCT_cluster <- as.character(tmp_object$SCT_cluster)
tmp_object$SCT_cluster[which(tmp_object$SCT_sub_cluster == "3_1")] <- 28
tmp_object$SCT_cluster[which(tmp_object$SCT_sub_cluster1 == "4_2")] <- 29
tmp_object$SCT_cluster <- factor(tmp_object$SCT_cluster, levels = 0:29)
Idents(tmp_object) <- "SCT_cluster"
selected_cells <- CellsByIdentities(tmp_object)
Seurat::SpatialPlot(tmp_object, cells.highlight = selected_cells[setdiff(names(selected_cells), "NA")], cols.highlight = c("#FFFF00", alpha("lightgrey", 0)), crop = FALSE, shape = 22, facet.highlight = T, combine = T, ncol = 8, alpha = NULL) & NoLegend()

all_markers <- FindAllMarkers(tmp_object, only.pos = T, min.pct = 0.1)
E18_de_markers <- all_markers[which(all_markers$avg_log2FC > 1 & all_markers$p_val_adj < 0.05),]
Seurat::SpatialPlot(tmp_object, shape = 22, crop = F, ncol = 5, features = E18_de_markers$gene[which(E18_de_markers$cluster == 26)][1:15]) & scale_fill_gradientn(colours = paletteContinuous(set = "sambaNight"))

tmp_object$cell_identity <- as.character(tmp_object$SCT_cluster)
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(19))] <- "Cortex_AP" # "Pax6", "Celsr1", "Ccdc80", "Lrp4"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(23))] <- "Cortex_IPC_MigN" # "Cemip", "Fgfr1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(13))] <- "Cortex_IPC" # "Unc5d", "Sema3c", "Slc17a6"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(2))] <- "Cortex_MigN" # "Adamts2", "Nrp1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(25))] <- "Cortex_Layer_6b" # "Tmem132d", "Cdh18"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(10, 21))] <- "Cortex_CThPN" # "Wnt7b", "Nef2l3", "Hs3st4", "Necab1", "Lmo7"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(14, 20))] <- "Cortex_SCPN" # "Inhba", "Mlip", "Adamts20"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(15))] <- "Cortex_UL_CPN" # "Tmtc1", "Ccbe1", "Kif26b"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(26))] <- "Cortex_CR" # "Reln"

tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(24))] <- "MGE_LGE_LS_AP"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(28))] <- "MGE_InN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(22))] <- "LGE_IPC" # "Dlx2", "Cdca7"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(8))] <- "LGE_InN" # "Gad2" "Slit3"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(6, 18))] <- "CPU_pre_MSN" # "Ikzf1", "Dach1", "Cntnap3"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(0))] <- "CPU_MSN" # "Rxrg", "Rarb", "Caln1", "Myo3b", "Oprm1", "Brinp1", "Kcnh1"

tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(27))] <- "CLA_ExN" # "Nr4a2", "Sv2b"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(16))] <- "PIR_ExN" # "Tafa1", "Ppp1r14c", "Thsd7a"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(11))] <- "PAL" # "Acvr2a" "Pcdh10"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(1))] <- "PAL" # "Mgat4c", "Tenm1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(3))] <- "PAL" # "Lhx6", "Rgma", "Kcnc2", "Pde1a"

tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(17))] <- "LS_IPC" # "Zic2", "Zic4", "Slit1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(4))] <- "LS" # "Trpc4", "Fgd5", "Zic4", "Zic1", "Peg3"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(9))] <- "MS" # "Ntng1"

tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(5))] <- "DB" # "Nefm", "Samd5", "Sulf2", "Galntl6"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(29))] <- "AVPV" #

tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(12))] <- "Endo" # "Flt1", "Adgrl4", "Kdr"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(7))] <- "Stromal" # "Dcn", "Col6a3"

tmp_object$cell_identity <- factor(tmp_object$cell_identity, levels = c("Cortex_AP", "Cortex_IPC_MigN", "Cortex_IPC", "Cortex_MigN", "Cortex_Layer_6b", "Cortex_CThPN", "Cortex_SCPN", "Cortex_UL_CPN", "Cortex_CR", "MGE_LGE_LS_AP", "MGE_InN", "LGE_IPC", "LGE_InN", "CPU_pre_MSN", "CPU_MSN", "CLA_ExN", "PIR_ExN", "PAL", "LS_IPC", "LS", "MS", "DB", "AVPV", "Endo", "Stromal"))
Seurat::SpatialPlot(tmp_object, group.by = "cell_identity", shape = 22, cols = cluster_colors, crop = F)
#!------------------------------
Idents(tmp_object) <- "cell_identity"
all_markers <- FindAllMarkers(tmp_object, min.pct = 0.1, only.pos = T)
de_markers <- all_markers[which(all_markers$avg_log2FC > 1 & all_markers$p_val_adj < 0.01 & all_markers$pct.1 > 0.3),]
write.csv(de_markers, "brain_20um_analysis/E18_5_ST_cell_identity_DE_markers.csv", quote = F)

key_markers <- c("Satb2", "Top2a", "Mki67", "Eomes", "Pax6", "Veph1", "Ccdc80", "Aldoc", "Neurog2", "Dct", "Pter", "Hopx", "Ltbp1", "Sema3c", "Slc17a6", "Pou3f2", "Nrn1", "Prox1", "Ndst3", "Lypd6", "Adamts2", "Fst", "Calb2", "Cemip", "Nr3c2", "Abca8a", "Tmem132d", "Pappa2", "Nfe2l3", "Lmo7", "Necab1", "Cdh9", "Sv2b", "Klhl1", "Npy", "Rspo3", "Npr3", "Pde1a", "Nr4a3", "Etv5", "Nwd2", "Tafa1", "Tcerg1l", "Adamts20", "Rorb", "L3mbtl4", "Rasgrf2", "Arhgap28", "Mlip", "Mkx", "Ntf3", "Usp43", "Inhba", "Tafa2", "Bcl6", "Pdzrn4", "Mafb", "Nxph2",  "Npas1", "Nxph1", "Lhx6", "Dgkg", "Lhcgr", "Vcam1", "Pdzph1", "Hepacam", "Mlc1", "Dlx2", "Sp8", "Akna", "Smoc1", "St18", "Cdca7", "Sp9", "Dlx6", "Dlx1", "Mob3b", "Six3", "Drd2", "Dchs2", "Ikzf1", "Cntnap3", "Dclk3", "Ano3", "Rxrg", "Oprm1", "Myo3b", "Tgfa", "Ngef", "Gpr88", "Cldn5", "Tmem72", "Folr1", "Clic6", "Gmnc", "Tfap2d", "Pdgfra", "Cspg4", "Pstpip2", "Prkcq", "Shc4", "Dynlrb2", "Daw1", "Crocc2", "Deup1", "Wdr63", "Zic4", "Zic1", "Snhg11", "Fgd5", "Flt1", "Eng", "Epas1", "Igfbp7", "Adgrl4", "C1qb", "Cd180", "Dock2", "Spp1")

DotPlot(tmp_object, features = key_markers, col.min = 0.1, dot.min = 0.1) + RotatedAxis() + labs(y = "cell identity", x = "Features") + scale_color_gradientn(colours = paletteContinuous(set = "whitePurple"))
ggsave("brain_20um_analysis/E18_5_cell_identity_key_markers.pdf")
#!----------------------------
tmp_object$Cluster <- ClusterID[match(as.character(tmp_object$cell_identity), ClusterNames)]
tmp_object$Cluster <- factor(tmp_object$Cluster, levels = intersect(ClusterID, unique(tmp_object$Cluster)))
Seurat::SpatialPlot(tmp_object, group.by = "Cluster", shape = 22, cols = ClusterIDColors, crop = F)

Idents(tmp_object) <- "cell_identity"
seurat_list$E18_5 <- tmp_object
#!-----------------------------------------------------------------------
seurat_object <- merge(seurat_list[[1]], seurat_list[2:length(seurat_list)], add.cell.ids = names(seurat_list))

DefaultAssay(seurat_object) <- "spARC_SCT"
seurat_object <- JoinLayers(seurat_object)

DefaultAssay(seurat_object) <- "RNA"
seurat_object$orig.ident <- factor(seurat_object$orig.ident, levels = names(seurat_list))
Idents(seurat_object) <- "orig.ident"

VlnPlot(seurat_object, features = c('nCount_RNA', 'nFeature_RNA'), pt.size = 0, ncol = 2) + NoLegend()
#!----------------------------------------------------------------------
seurat_object <- SCTransform(seurat_object, vst.flavor = "v2", method = "glmGamPoi", vars.to.regress = c("S.Score", "G2M.Score"))

seurat_object <- RunPCA(seurat_object, npcs = 50, verbose = FALSE)
seurat_object <- RunUMAP(seurat_object, reduction = "pca", dims = 1:50, reduction.name = "umap_SCT", return.model = T)

seurat_object <- IntegrateLayers(object = seurat_object, method = CCAIntegration, orig.reduction = "pca", new.reduction = 'cca', normalization.method = "SCT", verbose = FALSE)
seurat_object <- RunUMAP(seurat_object, reduction = "cca", dims = 1:50, reduction.name = "umap_cca", return.model = T)

# one-liner to run Integration
seurat_object <- IntegrateLayers(object = seurat_object, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = 'harmony', normalization.method = "SCT", verbose = FALSE)
seurat_object <- RunUMAP(seurat_object, reduction = "harmony", dims = 1:50, reduction.name = "umap_harmony", return.model = T)

p0 <- DimPlot(seurat_object, cols = mypal, reduction = "umap_SCT") + NoLegend()
p1 <- DimPlot(seurat_object, cols = mypal, reduction = "umap_cca") + NoLegend()
p2 <- DimPlot(seurat_object, cols = mypal, reduction = "umap_harmony")

p0 + p1 + p2

seurat_object <- FindNeighbors(seurat_object, reduction = "harmony", dims = 1:50, annoy.metric = "cosine", graph.name = c('harmony_nn', 'harmony_snn'))
seurat_object <- FindClusters(seurat_object, resolution = seq(0.1, 2, 0.1), graph.name = 'harmony_snn')
#!---------------------------------------------------------------------
library(SeuratWrappers)
library(Banksy)
library(harmony)

locs <- c()
for(i in names(seurat_object@images)){
	locs <- rbind(locs, Seurat::GetTissueCoordinates(seurat_object, image = i)[,seq_len(2)])
}
all(rownames(locs) == colnames(seurat_object))
colnames(locs) <- c("sdimy", "sdimx")
seurat_object@meta.data <- cbind(seurat_object@meta.data, locs)
seurat_object$orig.ident <- factor(seurat_object$orig.ident, levels = unique(seurat_object$orig.ident))

resolution <- seq(0.1, 2, 0.1)
DefaultAssay(seurat_object) <- "SCT"
seurat_object <- RunBanksy(seurat_object, lambda = 0.2, assay = 'SCT', slot = 'data', dimx = "sdimx", dimy = "sdimy", features = 'variable', group = 'orig.ident', split.scale = TRUE, k_geom = 15)
DefaultAssay(seurat_object) <- "BANKSY"
seurat_object <- RunPCA(seurat_object, assay = "BANKSY", reduction.name = "pca_banksy", features = rownames(seurat_object), npcs = 30)
seurat_object <- RunHarmony(seurat_object, group.by.vars = "orig.ident",  reduction.name = 'pca_banksy', reduction.save = 'banksy_harmony')

seurat_object <- RunUMAP(seurat_object, reduction = "pca_banksy", dims = 1:30, reduction.name = "umap_banksy", return.model = T)
seurat_object <- RunUMAP(seurat_object, reduction = "banksy_harmony", dims = 1:30, reduction.name = "umap_banksy_harmony", return.model = T)

p3 <- DimPlot(seurat_object, cols = mypal, group.by = "orig.ident", reduction = "umap_banksy") + NoLegend()
p4 <- DimPlot(seurat_object, cols = mypal, group.by = "orig.ident", reduction = "umap_banksy_harmony")

#seurat_object <- FindNeighbors(seurat_object, reduction = "pca", dims = 1:50, graph.name = c('banksy_nn', 'banksy_snn'), annoy.metric = "cosine")
#seurat_object <- FindClusters(seurat_object, resolution = resolution, graph.name = "banksy_snn")

seurat_object <- FindNeighbors(seurat_object, dims = 1:30, reduction = 'harmony', annoy.metric = "cosine", graph.name = c('banksy_harmony_nn', 'banksy_harmony_snn'))
seurat_object <- FindClusters(seurat_object, resolution = resolution, graph.name = 'banksy_harmony_snn')

DefaultAssay(seurat_object) <- "RNA"
seurat_object <- JoinLayers(seurat_object)

DefaultAssay(seurat_object) <- "SCT"
library(purrr)
wrap_plots(map(seq(0.6, 2, 0.1), function(x) DimPlot(seurat_object, reduction = "umap_banksy_harmony", group.by = paste0("banksy_harmony_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)
#!-----------------------------------------------------------------
# development merge annotation
wrap_plots(map(seq(0.6, 2, 0.1), function(x) DimPlot(seurat_object, reduction = "umap_harmony", group.by = paste0("harmony_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)

Idents(seurat_object) <- "harmony_snn_res.1.6"
dir.create("SCT_harmony_20um_merge_Res")
selected_cells <- CellsByIdentities(seurat_object)
for(i in 1:length(selected_cells)){
	p0 <- Seurat::SpatialPlot(seurat_object, cells.highlight = selected_cells[i], cols.highlight = c("#FFFF00", alpha("lightgrey", 0)), crop = FALSE, shape = 22, combine = T, ncol = 3, alpha = NULL) & NoLegend()
	print(p0)
	ggsave(paste0("SCT_harmony_20um_merge_Res/", names(selected_cells)[i], ".png"))
}

Idents(seurat_object) <- "banksy_harmony_snn_res.1.3"
dir.create("SCT_banksy_harmony_20um_merge_Res")
selected_cells <- CellsByIdentities(seurat_object)
for(i in 1:length(selected_cells)){
	p0 <- Seurat::SpatialPlot(seurat_object, cells.highlight = selected_cells[i], cols.highlight = c("#FFFF00", alpha("lightgrey", 0)), crop = FALSE, shape = 22, combine = T, ncol = 3, alpha = NULL) & NoLegend()
	print(p0)
	ggsave(paste0("SCT_banksy_harmony_20um_merge_Res/", names(selected_cells)[i], ".png"))
}

Seurat::SpatialPlot(seurat_object, features = "Sftpc", crop = FALSE, shape = 22, alpha = c(0.3, 1)) + scale_fill_gradientn(colours = paletteContinuous(set = "sambaNight"))
Seurat::SpatialPlot(seurat_object, features = i, crop = FALSE, shape = 22, alpha = c(0.3, 1)) + scale_fill_gradientn(colours = paletteContinuous(set = "sambaNight")) + ggtitle("log(UMI)")

seurat_object$seurat_clusters <- seurat_object$harmony_snn_res.1.6
Idents(seurat_object) <- "seurat_clusters"

all_markers <- FindAllMarkers(seurat_object, only.pos = T, min.pct = 0.1)
de_markers <- all_markers[which(all_markers$avg_log2FC > 1 & all_markers$p_val_adj < 0.05),]
Seurat::SpatialPlot(seurat_object, shape = 22, crop = F, ncol = 5, features = de_markers$gene[which(de_markers$cluster == 26)][1:15]) & scale_fill_gradientn(colours = paletteContinuous(set = "sambaNight"))

seurat_object$cell_identity <- as.character(seurat_object$seurat_clusters)
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(11))] <- "Cortex_AP"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(30))] <- "Cortex_IPC_MigN"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(2))] <- "Cortex_IPC"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(10))] <- "Cortex_MigN"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(7))] <- "Cortex_CThPN"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(19, 20))] <- "Cortex_SCPN"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(9))] <- "Cortex_UL_CPN"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(23))] <- "Cortex_CR"

seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(25, 27))] <- "MGE_LGE_LS_AP"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(12))] <- "MGE_LGE_IPC"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(14))] <- "MGE_InN"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(15))] <- "LGE_InN"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(5))] <- "CPU_pre_MSN"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(4))] <- "CPU_MSN"

seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(26))] <- "CLA_ExN"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(13))] <- "PIR_ExN"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(22))] <- "LCS_IMN"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(31))] <- "Epd_ExN"

seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(3))] <- "PLA"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(8))] <- "PLA"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(16))] <- "PLA"

seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(32))] <- "LS_IPC"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(6))] <- "LS"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(21))] <- "LS"

seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(28))] <- "MS"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(17))] <- "DB"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(0))] <- "AVPV"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(24))] <- "POA"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(29))] <- "IIIV"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(33))] <- "Endo"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(1, 34, 35))] <- "Stromal"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(18))] <- "others"

seurat_object$cell_identity <- factor(seurat_object$cell_identity, levels = c("Cortex_AP", "Cortex_IPC_MigN", "Cortex_IPC", "Cortex_MigN", "Cortex_CThPN", "Cortex_SCPN", "Cortex_UL_CPN", "Cortex_CR", "MGE_LGE_LS_AP", "MGE_LGE_IPC", "MGE_InN", "LGE_InN", "CPU_pre_MSN", "CPU_MSN", "CLA_ExN", "PIR_ExN", "LCS_IMN", "Epd_ExN", "PLA", "LS_IPC", "LS", "MS", "DB", "AVPV", "POA", "IIIV", "Endo", "Stromal", "others"))
seurat_object$cell_identity0 <- seurat_object$cell_identity
#!-----------------------------------------------------------------------------------------------------------------------------------------------------------
seurat_object@misc$cluster_colors <- cluster_colors
seurat_object@misc$ClusterIDColors <- ClusterIDColors

tmp <- seurat_object@meta.data[,c("staggered_sdimx", "staggered_sdimy")]
colnames(tmp) <- c("spatial_umap_1", "spatial_umap_2")
tmp$spatial_umap_2 <- 1080 - tmp$spatial_umap_2
seurat_object[["spatial_umap"]] <- CreateDimReducObject(embeddings = as.matrix(tmp), assay = "SCT")

E14_embeding <- seurat_list$E14_5@reductions$SCT_umap@cell.embeddings
rownames(E14_embeding) <- paste0("E14_5_", rownames(E14_embeding))

E16_embeding <- seurat_list$E16_5@reductions$SCT_umap@cell.embeddings
rownames(E16_embeding) <- paste0("E16_5_", rownames(E16_embeding))
E16_embeding[,1] <- E16_embeding[,1] + max(E14_embeding[,1]) + 30

E18_embeding <- seurat_list$E18_5@reductions$SCT_umap@cell.embeddings
rownames(E18_embeding) <- paste0("E18_5_", rownames(E18_embeding))
E18_embeding[,1] <- E18_embeding[,1] + max(E16_embeding[,1]) + 30

tmp_embeding <- rbind(E14_embeding, E16_embeding, E18_embeding)
colnames(tmp_embeding) <- c("spatial_per_umap_1", "spatial_per_umap_2")
seurat_object[["spatial_per_umap"]] <- CreateDimReducObject(embeddings = as.matrix(tmp_embeding), assay = "SCT")

Idents(seurat_object) <- "cell_identity"
DefaultAssay(seurat_object) <- "SCT"

library(CytoTRACE2)
seurat_object <- cytotrace2(seurat_object, is_seurat = TRUE, slot_type = "counts", species = 'mouse')

p1 <- plot_cytotrace2(seurat_object, reduction = "umap_harmony")

(p1$CytoTRACE2_UMAP + p1$CytoTRACE2_Potency_UMAP + p1$CytoTRACE2_Relative_UMAP) & SetAxes()
ggsave("brain_20um_analysis/ST_development_CytoTRACE2.pdf")
Seurat::SpatialPlot(seurat_object, shape = 22, features = "CytoTRACE2_Relative", crop = F)  & scale_fill_gradientn(colours = paletteContinuous(set = "beach"))
ggsave("brain_20um_analysis/ST_development_CytoTRACE2_Spatial.pdf")
qsave(seurat_object, "output/brain_development_ST_final_tutorial_20um.qs")
#!------------------------------------------------------------------------------------------------------------------------------------------------------------
key_markers <- c("Pax6", "Eomes", "Hmga2", "Lrp4", "Dct", "Adamts19", "Eya1", "Nrn1", "Slc17a6", "Sema3c", "Scrt2", "Lhx2", "Bcar1", "Adamts2", "Pou3f2", "Unc5d", "Fgfr1", "Cemip", "Tbr1", "Hs3st4", "Nfe2l3", "Tmem132d", "Pappa2", "Kcnab1", "Abca8a", "Npas2", "Wnt7b", "Klhl1", "Npy", "Kctd12", "Nr4a3", "Slc24a2", "Pde1a", "Fezf2", "Kcnk1", "Cntn6", "Tafa1", "Tcerg1l", "Rorb", "Mical2", "Trmt9b", "Ntf3", "Id2", "Tafa2", "Bcl6", "Inhba", "Reln", "Dgkg", "Elfn1", "Lhx6", "Rai2", "Lhx8", "Nkx2-1", "Pdzph1", "Daam2", "Hes5", "Ascl1", "Dlx6", "Dlx1", "Dlx2", "Sp9", "Slit3", "Six3", "Drd2", "Sox8", "Nebl", "Ikzf1", "Drd1", "Tac1", "Ngef", "Gmnc", "Lmx1a", "Aqp1", "Tfap2d", "Nr2f1", "Nr4a2", "Sv2b", "Slit2", "Zic4", "Trhde", "Galnt13", "Lama4", "Ets1", "Ptprb", "Col3a1", "Col5a1", "Dcn", "Ptprc", "Ly86", "Csf1r")

DotPlot(seurat_object, features = key_markers, cols = c("lightgrey", "red"), col.min = 0.1, dot.min = 0.1) + RotatedAxis() + labs(y = "cell identity", x = "Features") + NoLegend()
#!----------------------------------------------------------------------------------------------------------------------------------------------------------
seurat_object <- subset(seurat_object, subset = cell_identity %in% levels(seurat_object$cell_identity)[1:11])

seurat_object <- shiny_st(seurat = seurat_object, isVisium = F, assay = "RNA", image = "E14_5")
seurat_object <- subset(seurat_object, subset = cell_identity != "removed_cells")
seurat_object$cell_identity <- factor(seurat_object$cell_identity, levels = intersect(levels(seurat_object$cell_identity), unique(seurat_object$cell_identity)))

seurat_object$Region <- "CP"
seurat_object$Region[which(seurat_object$cell_identity %in% c("Cortex_AP", "Cortex_IPC", "Cortex_MigN"))] <- "GZ"
seurat_object$Region <- factor(seurat_object$Region, levels = c("GZ", "CP"))

Seurat::SpatialPlot(seurat_object, shape = 22, crop = F, cols = ClusterIDColors, group.by = "Cluster")
ggsave("brain_20um_analysis/ST_Cortex_banksy_cluster.pdf")
Seurat::SpatialPlot(seurat_object, shape = 22, crop = F, group.by = "Region", cols = cluster_colors)
ggsave("brain_20um_analysis/ST_Cortex_Region.pdf")

qsave(seurat_object, "output/brain_development_ST_final_tutorial_cortex_20um.qs")
#!----------------------------------------------------------------------------------------------------------------------------------------------------------
hic_object <- qread("output/brain_development_Spatial_HiC_final_tutorial_20um.qs")
hic_object <- Split_Layers(hic_object, assay = "scAB500kb_scale", split.by = "orig.ident")

library(reticulate)
use_python("/home/yiyelinfeng/softwares/miniconda3/envs/R4/bin")
SPARC <- import("SPARC", convert = FALSE)

data_matrix <- Seurat::GetAssayData(object = hic_object, slot = "data.E14_5", assay = "scAB500kb_scale")
data_sparc <- SPARC$spARC()$fit_transform(expression_X = t(x = as.matrix(data_matrix)), spatial_X = hic_object@meta.data[colnames(data_matrix),c("sdimx", "sdimy")])
E14_matrix <- t(as.matrix(x = data_sparc))
colnames(E14_matrix) <- colnames(data_matrix)
rownames(E14_matrix) <- rownames(data_matrix)

data_matrix <- Seurat::GetAssayData(object = hic_object, slot = "data.E16_5", assay = "scAB500kb_scale")
data_sparc <- SPARC$spARC()$fit_transform(expression_X = t(x = as.matrix(data_matrix)), spatial_X = hic_object@meta.data[colnames(data_matrix),c("sdimx", "sdimy")])
E16_matrix <- t(as.matrix(x = data_sparc))
colnames(E16_matrix) <- colnames(data_matrix)
rownames(E16_matrix) <- rownames(data_matrix)

data_matrix <- Seurat::GetAssayData(object = hic_object, slot = "data.E18_5", assay = "scAB500kb_scale")
data_sparc <- SPARC$spARC()$fit_transform(expression_X = t(x = as.matrix(data_matrix)), spatial_X = hic_object@meta.data[colnames(data_matrix),c("sdimx", "sdimy")])
E18_matrix <- t(as.matrix(x = data_sparc))
colnames(E18_matrix) <- colnames(data_matrix)
rownames(E18_matrix) <- rownames(data_matrix)

data_matrix <- cbind(E14_matrix, E16_matrix, E18_matrix)
hic_object[["spARC_scAB500kb_scale"]] <- CreateAssay5Object(data = data_matrix[,colnames(hic_object)])
hic_object <- JoinLayers(hic_object, assay = "scAB500kb_scale")

DefaultAssay(seurat_object) <- "spARC_SCT"
DefaultAssay(hic_object) <- "spARC_scAB500kb_scale"
seurat_object$shape = "RNA"
hic_object$shape = "HiC"
#!-----------------------------------------------------------------------------------------------------------------------------------------
DefaultAssay(seurat_object) <- "spARC_SCT"
DefaultAssay(hic_object) <- "spARC_scAB500kb_scale"
seurat_object$shape = "RNA"
hic_object$shape = "HiC"

for(i in c("Sox2", "Pax6", "Notch2", "Rbfox1", "Ank3", "Celf2", "Eomes", "Satb2", "Neurod2", "Dlx2", "Gad2", "Foxp1", "Pde10a")){
	p0 <- Seurat::SpatialPlot(seurat_object, features = i, shape = 22, crop = F) & scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")))
	p1 <- Seurat::SpatialPlot(hic_object, shape = 22, crop = F, features = unique(gene_infor_500k$gene_region[which(gene_infor_500k$gene_name == i & gene_infor_500k$is_promoter == "promoter")])) & scale_fill_gradientn(values = seq(0, 1, 0.01), colours = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"))
	print(p0/p1)
	ggsave(paste0("brain_20um_analysis/", i, "_Spatial.pdf"))

	p0 <- FeaturePlot(seurat_object, reduction = "spatial_umap", features = i, shape.by = "shape", pt.size = 0.8) & scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) & scale_shape_manual(values = 15) & ylim(c(0, 1080)) & xlim(c(0, 5500))
	p1 <- FeaturePlot(hic_object, features = unique(gene_infor_500k$gene_region[which(gene_infor_500k$gene_name == i & gene_infor_500k$is_promoter == "promoter")]), reduction = "spatial_umap", shape.by = "shape", pt.size = 0.8) & scale_colour_gradientn(values = seq(0, 1, 0.01), colours = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), na.value = "white") & scale_shape_manual(values = 15) & ylim(c(0, 1080)) & xlim(c(0, 5500))
	print(p0/p1)
	ggsave(paste0("brain_20um_analysis/", i, "_Spatial0.pdf"))
}
#!-----------------------------------------------------------------------------------------------------------------------------------------
DefaultAssay(seurat_object) <- "spARC_SCT"
DefaultAssay(hic_object) <- "spARC_scAB500kb_scale"
seurat_object$shape = "RNA"
hic_object$shape = "HiC"

for(i in c("Tmem132d", "Tafa1", "Cdh13", "Gnal", "Unc5d", "Sox2", "Pax6", "Tafa2", "Meis2")){
	p0 <- Seurat::SpatialPlot(seurat_object, features = i, shape = 22, crop = F) & scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")))
	p1 <- Seurat::SpatialPlot(hic_object, shape = 22, crop = F, features = unique(gene_infor_500k$gene_region[which(gene_infor_500k$gene_name == i & gene_infor_500k$is_promoter == "promoter")])) & scale_fill_gradientn(values = seq(0, 1, 0.01), colours = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"))
	print(p0/p1)
	ggsave(paste0("brain_20um_analysis/", i, "_Spatial.pdf"))

	p0 <- FeaturePlot(seurat_object, reduction = "spatial_umap", features = i, shape.by = "shape", pt.size = 0.8) & scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) & scale_shape_manual(values = 15) & ylim(c(0, 1080)) & xlim(c(0, 5500))
	p1 <- FeaturePlot(hic_object, features = unique(gene_infor_500k$gene_region[which(gene_infor_500k$gene_name == i & gene_infor_500k$is_promoter == "promoter")]), reduction = "spatial_umap", shape.by = "shape", pt.size = 0.8) & scale_colour_gradientn(values = seq(0, 1, 0.01), colours = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), na.value = "white") & scale_shape_manual(values = 15) & ylim(c(0, 1080)) & xlim(c(0, 5500))
	print(p0/p1)
	ggsave(paste0("brain_20um_analysis/", i, "_Spatial0.pdf"))
}
#!-----------------------------------------------------------------------------------------------------------------------------------------
i = "Tafa1"

p0 <- Seurat::SpatialPlot(seurat_object, features = i, shape = 22, crop = F) & scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")))
p1 <- Seurat::SpatialPlot(hic_object, shape = 22, crop = F, features = unique(gene_infor_500k$gene_region[which(gene_infor_500k$gene_name == i & gene_infor_500k$is_promoter == "promoter")])) & scale_fill_gradientn(values = seq(0, 1, 0.01), colours = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"))
print(p0/p1)
ggsave(paste0("brain_20um_analysis/", i, "_Spatial.pdf"))

p0 <- FeaturePlot(seurat_object, reduction = "spatial_umap", features = i, shape.by = "shape", pt.size = 0.8) & scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) & scale_shape_manual(values = 15) & ylim(c(0, 1080)) & xlim(c(0, 5500))
p1 <- FeaturePlot(hic_object, features = unique(gene_infor_500k$gene_region[which(gene_infor_500k$gene_name == i & gene_infor_500k$is_promoter == "promoter")]), reduction = "spatial_umap", shape.by = "shape", pt.size = 0.8) & scale_colour_gradientn(values = seq(0, 1, 0.01), colours = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), na.value = "white") & scale_shape_manual(values = 15) & ylim(c(0, 1080)) & xlim(c(0, 5500))
print(p0/p1)
ggsave(paste0("brain_20um_analysis/", i, "_Spatial0.pdf"))

DefaultAssay(seurat_object) <- "SCT"
DefaultAssay(hic_object) <- "scAB500kb_scale"

p0 <- VlnPlot(seurat_object, features = i, pt.size = 0, group.by = "orig.ident", split.by = "Region", cols = c("#66C2A5", "#FC8D62"))
p1 <- VlnPlot(hic_object, pt.size = 0, features = unique(gene_infor_500k$gene_region[which(gene_infor_500k$gene_name == i & gene_infor_500k$is_promoter == "promoter")]), group.by = "orig.ident", split.by = "Region", cols = c("#66C2A5", "#FC8D62")) + ylab("scAB values") & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.3, position = position_dodge(0.9), show.legend = F)
print(p0 + p1)
ggsave(paste0("brain_20um_analysis/", i, "_vlnplot.pdf"))

P_values <- c()
for(j in levels(hic_object$orig.ident)){
	GZ <- GetAssayData(hic_object, assay = "scAB500kb_scale", layer = "data")[unique(gene_infor_500k$gene_region[which(gene_infor_500k$gene_name == i & gene_infor_500k$is_promoter == "promoter")]), which(hic_object$orig.ident == j & hic_object$Region == "GZ")]
	CP <- GetAssayData(hic_object, assay = "scAB500kb_scale", layer = "data")[unique(gene_infor_500k$gene_region[which(gene_infor_500k$gene_name == i & gene_infor_500k$is_promoter == "promoter")]), which(hic_object$orig.ident == j & hic_object$Region == "CP")]
	P_valves_t <- t.test(GZ, CP)[3]
	P_valves_w <- wilcox.test(GZ, CP)[3]
	P_values <- rbind(P_values, c(i, j, P_valves_t, P_valves_w))
}
colnames(P_values) <- c("feature", "group", "t_test_P_value", "wilcox_test_P_value")
write.csv(P_values, paste0("brain_20um_analysis/", i, "_HiC_mean_scAB_values_P_values.csv"), quote = F, row.names = F)

VlnPlot(seurat_object, features = i, pt.size = 0, group.by = "orig.ident", split.by = "Cluster", cols = ClusterIDColors)
ggsave(paste0("brain_20um_analysis/", i, "_SCT_Cluster_per_time_vlnplot.pdf"))
VlnPlot(seurat_object, features = i, pt.size = 0, group.by = "Cluster", cols = ClusterIDColors) + NoLegend()
ggsave(paste0("brain_20um_analysis/", i, "_SCT_Cluster_vlnplot.pdf"))
VlnPlot(hic_object, pt.size = 0, features = unique(gene_infor_500k$gene_region[which(gene_infor_500k$gene_name == i & gene_infor_500k$is_promoter == "promoter")]), group.by = "banksy_cluster", cols = hic_object@misc$banksy_cluster_colors) + ylab("scAB values") + NoLegend() & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.3, show.legend = F)
ggsave(paste0("brain_20um_analysis/", i, "_scAB_banksy_cluster_vlnplot.pdf"))















