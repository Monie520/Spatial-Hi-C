source("/home/yiyelinfeng/scripts/Rscripts/scRNASeqSeuratAnalysisFunction.R")
source("/home/yiyelinfeng/scripts/Rscripts/R.analysis.Functions.R")
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
library(purrr)
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

cluster_colors <- c("E14_5" = "#34D916", "E16_5" = "#00D4E6", "E18_5" = "#1E90FF", "Adult15" = "#B312A6", "GZ" = "#00D4E6", "CP" = "#00998F",
"Cortex" = "#59B375", "MGE" = "#F97B72", "LGE" = "#89288F", "ChP" = "#FF26A8", "LCS" = "#C1FF73", "CLA" = "#33FF00", "LS" = "#E68316", "Epd" = "#94FFB5", "PIR" = "#B5EFB5", "AEP" = "#653EB3", "others" = "#FEE52C", "Others" = "#dde2e6",
"AP" = "#03FFF4", "Cortex_DP_AP" = "#03FFF4", "Cortex_MP_AP" = "#036DF4", "IPC" = "#0BD3B1", "Cortex_IPC" = "#0BD3B1", "MigN" = "#62CFE8", "Cortex_DP_MigN" = "#62CFE8", "Cortex_MP_MigN" = "#2F7DD1",
"Cortex_Layer_6b" = "#99FFBC", "Layer_6b" = "#99FFBC", "CThPN" = "#7EC136", "Cortex_CThPN" = "#7EC136", "Cortex_DP_CThPN" = "#00CC14", "Cortex_MP_CThPN" = "#2D739E", "SCPN" = "#34A047", "Cortex_SCPN" = "#34A047", "Cortex_DP_SCPN" = "#00991F", "Cortex_MP_SCPN" = "#5C5CA6", "UL_CPN" = "#01545a", "Cortex_UL_CPN" = "#01545a",
"CTGL_ITL6GL" = "#66C5CC", "PTGL" = "#C9DB74", "ITL5GL" = "#87C55F", "ITL4GL" = "#54990F", "ITL23GL" = "#017351", "ITL1" = "#3DCCB1", "Astrocytes" = "#16F2F2", "CR" = "#F2F318", "Cortex_CR" = "#F2F318", "MGE_derived_InN" = "#D4E727", "Cortex_MGE_derived_InN" = "#D4E727",
"Layer_I" = "#3DCCB1", "Layer_II_III" = "#017351", "Layer_IV" = "#54990F", "Layer_V" = "#87C55F", "Layer_VI" = "#66C5CC", "CC" = "#64C2FC", "CPU" = "#D38B5C",
"VZ0" = "#4DFF99", "VZ" = "#03FFF4", "SVZ" = "#52FF16", "IZ" = "#00CC14", "Layer_VIb" = "#99FFBC", "Layer_VI" = "#66C5CC", "Layer_V" = "#1FCCCC", "Layer_II_III_IV" = "#006B99", "CP" = "#009E73",
"AP0" = "#4DFF99", "AP" = "#03FFF4", "IPC_MigN" = "#52FF16", "CThPN_SCPN" = "#59B375","SCPN_UL_CPN" = "#009E73",
"ODC" = "", "ODC1" = "#9EFF99", "ODC2" = "#64C2FC", "SSTGA" = "#CC79A7", "PVGA" = "#FF26CB", "VIPGA" = "#267DFF", 
"LGE_MGE_AP" = "#9E7BFF", "LGE_AP" = "#9E7BFF", "MGE_IPC" = "#DEA0FD", "MGE_InN" = "#BE50FF", "LGE_IPC" = "#D85FF7", "LGE_InN" = "#AA0DFE", "CPU_pre_MSN" = "#F07F92", "CPU_MSN" = "#D38B5C", "CPU_Endo" = "#FF9D4B",
"ChP_MCC" = "#FF26A8", "LCS_IMN" = "#C1FF73", "CLA_ExN" = "#33FF00", "Epd_ExN" = "#94FFB5", "PIR_ExN" = "#B5EFB5", "LS_AP" = "#000080", "LS_IPC" = "#CC79A7", "AEP_AP" = "#653EB3",
"VLMC" = "8C8C8C", "Endo" = "#994567", "Fibrob" = "#AA9F0D", "Stromal" = "#886C00", "MG" = "#3D3D3D", "unknown" = "#dde2e6")

ClusterNames <- c("Cortex_DP_AP", "Cortex_MP_AP", "Cortex_IPC", "Cortex_DP_MigN", "Cortex_MP_MigN", "Cortex_Layer_6b", "Cortex_CThPN", "Cortex_DP_CThPN", "Cortex_MP_CThPN", "Cortex_SCPN", "Cortex_DP_SCPN", "Cortex_MP_SCPN", "Cortex_UL_CPN", "Cortex_CR", "LGE_MGE_AP", "LGE_AP", "MGE_IPC", "MGE_InN", "Cortex_MGE_derived_InN", "LGE_IPC", "LGE_InN", "CPU_pre_MSN", "CPU_MSN", "CPU_Endo", "ChP_MCC", "LCS_IMN", "Epd_ExN", "PIR_ExN", "CLA_ExN", "LS_AP", "LS_IPC", "AEP_AP", "Endo", "Fibrob", "MG", "Stromal")
ClusterID <- paste0("C", 1:length(ClusterNames))
ClusterIDColors <- cluster_colors[match(ClusterNames, names(cluster_colors))]
names(ClusterIDColors) <- ClusterID
names(ClusterNames) <- ClusterID

domainNames <- c("Cortex_1", "Cortex_2", "Cortex_3", "Cortex_4",  "Cortex_5", "MigN", "MGE", "LGE0", "LGE", "CPU", "ChP", "CLA", "LS", "End")
domainID <- paste0("D", 1:length(domainNames))
domainIDColors <- c("#03FFF4", "#61E2A4", "#3CB44B", "#35586D", "#B3A726", "#2F7DD1", "#6D32E6", "#BE50FF", "#AA0DFE", "#D38B5C", "#FC5151", "#26FF3E", "#450099", "#B38B5C")
names(domainIDColors) <- domainID
#!---------------------------------------------------------------------------------------------------------------------------------------------------
source("/home/yiyelinfeng/scripts/Rscripts/lung_project/IPF/spatial-lung-fibrosis/scripts/custom_colors.R")
seurat_list <- qread("output/seurat_hic_devp_list.qs")
seurat_RNA <- qread("output/seurat_list_SCT.qs")
options(scipen = 999)
gene_infor <- read.table("GRCm38_102_gtf_infor", header = T, sep = "\t")

bin = 250000

gene_infor$start <- floor((gene_infor$start - 1)/bin) * bin
gene_infor$end <- (floor((gene_infor$end - 1)/bin) + 1) * bin

gene_infor_tmp <- c()
for(i in 1:nrow(gene_infor)){
	if((gene_infor$end[i] - gene_infor$start[i]) > bin){
		tmp <- (gene_infor$end[i] - gene_infor$start[i]) / bin
		gene_infor_tmp <- rbind(gene_infor_tmp, c(gene_infor$chr[i], gene_infor$start[i], gene_infor$start[i] + bin, gene_infor$gene_id[i], gene_infor$gene_name[i], gene_infor$gene_type[i]))
		for(j in 2:tmp){
			gene_infor_tmp <- rbind(gene_infor_tmp, c(gene_infor$chr[i], gene_infor$start[i] + (j - 1) * bin, gene_infor$start[i] + j * bin,  gene_infor$gene_id[i], gene_infor$gene_name[i], gene_infor$gene_type[i]))
		}
	}else{
		gene_infor_tmp <- rbind(gene_infor_tmp, gene_infor[i,])
	}
}
gene_infor <- gene_infor_tmp
gene_infor$gene_region <- paste0(gene_infor$chr, ":", as.character(gene_infor$start), "-", as.character(gene_infor$end))
write.csv(gene_infor, "GRCm38_102_gtf_infor_hic_250kb_reference")
gene_infor_250k <- read.csv("GRCm38_102_gtf_infor_hic_250kb_reference")
gene_infor_500k <- read.csv("GRCm38_102_gtf_infor_hic_500kb_reference")
#!---------------------------------------------------------------------------
seurat_object <- seurat_list$E14_5
locs <- Seurat::GetTissueCoordinates(seurat_object)[,seq_len(2)]
colnames(locs) <- c("array_row", "array_col")
seurat_object@meta.data <- cbind(seurat_object@meta.data, locs)

DefaultAssay(seurat_object) <- "scAB500k"
VariableFeatures(seurat_object) <- rownames(seurat_object)

seurat_object <- RunBanksy(seurat_object, lambda = 0.3, assay = 'scAB500k', slot = 'data', dimx = "array_row", dimy = "array_col", features = 'all', k_geom = 24, assay_name = "BANKSY_scAB500k")
npcs = 30
seurat_object <- RunPCA(seurat_object, assay = "BANKSY_scAB500k", reduction.name = "pca_banksy", npcs = npcs, features = rownames(seurat_object))
seurat_object <- RunUMAP(seurat_object, reduction = "pca_banksy", dims = 1:npcs, min.dist = 0.3, reduction.name = "banksy_umap", return.model = T, n.neighbors = 24)
seurat_object <- FindNeighbors(seurat_object, reduction = "pca_banksy", dims = 1:npcs, annoy.metric = "cosine", graph.name = c('banksy_nn', 'banksy_snn'), k.param = 24)
seurat_object <- FindClusters(seurat_object, resolution = seq(0.1, 2, 0.1), graph.name = "banksy_snn")
wrap_plots(map(seq(0.6, 2, 0.1), function(x) DimPlot(seurat_object, reduction = "banksy_umap", group.by = paste0("banksy_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)
#!----------------------------------------------------------------------------
seurat_object <- FindClusters(seurat_object, cluster.name = "banksy_cluster", resolution = 1.7, graph.name = "banksy_snn")
seurat_object$banksy_cluster <- as.character(seurat_object$banksy_cluster)
seurat_object$banksy_cluster[intersect(which(seurat_object$banksy_cluster == 19), which(seurat_object@reductions$banksy_umap@cell.embeddings[,2] < -5))] <- 21
seurat_object$banksy_cluster[which(seurat_object$banksy_snn_res.2 == 21)] <- 22
seurat_object$banksy_cluster <- factor(seurat_object$banksy_cluster, levels = 0:22)

DimPlot(seurat_object, label = T, cols = selected_colors, group.by = "banksy_cluster", reduction = "banksy_umap") + NoLegend()
Seurat::SpatialPlot(seurat_object, cols = selected_colors, group.by = "banksy_cluster", crop = FALSE, shape = 22)
Idents(seurat_object) <- "banksy_cluster"
selected_cells <- CellsByIdentities(seurat_object)
Seurat::SpatialPlot(seurat_object, cells.highlight = selected_cells[setdiff(names(selected_cells), "NA")], cols.highlight = c("#FFFF00", alpha("lightgrey", 0)), crop = FALSE, shape = 22, facet.highlight = T, combine = T, ncol = 6, alpha = NULL) & NoLegend()

seurat_object$seurat_clusters <- seurat_object$banksy_cluster
seurat_object$banksy_cluster <- as.character(seurat_object$banksy_cluster)
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(2))] <- "C1"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(3))] <- "C2"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(16))] <- "C3"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(8))] <- "C4"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(1))] <- "C5"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(6, 20))] <- "C6"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(11, 22))] <- "C7"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(14))] <- "C8"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(0))] <- "C9"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(5))] <- "C10"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(7))] <- "C11"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(9))] <- "C12"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(4, 17))] <- "C13"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(19))] <- "C14"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(13))] <- "C15"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(12))] <- "C16"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(15))] <- "C17"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(10))] <- "C18"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(18, 21))] <- "C19" # Stromal

seurat_object$banksy_cluster <- factor(seurat_object$banksy_cluster, levels = paste0("C", 1:19))

Idents(seurat_object) <- "banksy_cluster"

E14_colors <- c("C1" = "#4DFF99", "C2" = "#1FCCCC", "C3" = "#0075DC", "C4" = "#7ACC00", "C5" = "#1FCC1F", "C6" = "#459967", "C7" = "#38B7FF", "C8" = "#967ACC", "C9" = "#AA0DFE", "C10" = "#782AB6", "C11" = "#7F00FF", "C12" = "#CC7A88", "C13" = "#B33E52", "C14" = "#00FFBE", "C15" = "#94FFB5", "C16" = "#16FF32", "C17" = "#B5EFB5", "C18" = "#3D0F99", "C19" = "#B3823E")
DimPlot(seurat_object, label = T, cols = E14_colors, reduction = "banksy_umap") + NoLegend() + SetAxes()
Seurat::SpatialPlot(seurat_object, cols = E14_colors, crop = FALSE, shape = 22)
seurat_object@misc$E14_colors <- E14_colors

seurat_rna <- seurat_RNA$E14_5
DefaultAssay(seurat_rna) <- "SCT"
all_markers <- FindAllMarkers(seurat_rna, only.pos = T, min.pct = 0.1)
write.csv(all_markers, "E14_5_ST_cell_identity.csv", quote = F)
de_markers <- all_markers[which(all_markers$avg_log2FC > 0.5 & all_markers$p_val_adj < 0.01 & all_markers$pct.1 > 0.3),]
marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:50]
	count = count + 1
}
names(marker_genes) <- unique(de_markers$cluster)

DefaultAssay(seurat_object) <- "scAB250k"
for(i in names(marker_genes)){
	marker_genes[[i]] <- intersect(rownames(seurat_object), unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name %in% marker_genes[[i]])]))
	seurat_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(seurat_object, assay = "scAB250k", slot = "data")[marker_genes[[i]],]))
}

tmp_meta <- as.matrix(seurat_object@meta.data[,levels(seurat_rna)])
rownames(tmp_meta) <- paste0(seurat_object$banksy_cluster)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- paste0(levels(seurat_rna$Cluster), "_", colnames(tmp_meta))
tmp_meta0 <- tmp_meta0[levels(seurat_object$banksy_cluster),]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[1], "_hic_250k_top50_zscore.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[1], "_hic_250k_top50_zscore1.pdf"), width = 10, height = 10)

marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:50]
	count = count + 1
}
names(marker_genes) <- unique(de_markers$cluster)
DefaultAssay(seurat_object) <- "scAB500k"
for(i in names(marker_genes)){
	marker_genes[[i]] <- intersect(rownames(seurat_object), unique(gene_infor_500k$gene_region[which(gene_infor_500k$gene_name %in% marker_genes[[i]])]))
	seurat_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(seurat_object, assay = "scAB500k", slot = "data")[marker_genes[[i]],]))
}

tmp_meta <- as.matrix(seurat_object@meta.data[,levels(seurat_rna)])
rownames(tmp_meta) <- paste0(seurat_object$banksy_cluster)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- paste0(levels(seurat_rna$Cluster), "_", colnames(tmp_meta))
tmp_meta0 <- tmp_meta0[levels(seurat_object$banksy_cluster),]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[1], "_hic_500k_top50_zscore.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[1], "_hic_500k_top50_zscore1.pdf"), width = 10, height = 10)

seurat_rna <- readRDS("output/E14_brain_mouse_scRNA_2021_Nature.rds")
Idents(seurat_rna) <- "cell_identity"
selected_clusters <- paste0("E14_5_", levels(seurat_rna$cell_identity))

all_markers <- FindAllMarkers(seurat_rna, only.pos = T, min.pct = 0.1)
write.csv(all_markers, "E14_5_single_cells_cell_identity_2021_Nature.csv", quote = F)
de_markers <- all_markers[which(all_markers$avg_log2FC > 0.5 & all_markers$p_val_adj < 0.01 & all_markers$pct.1 > 0.3),]
marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:50]
	count = count + 1
}
names(marker_genes) <- paste0("E14_5_", unique(de_markers$cluster))

DefaultAssay(seurat_object) <- "scAB250k"
for(i in names(marker_genes)){
	marker_genes[[i]] <- intersect(rownames(seurat_object), unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name %in% marker_genes[[i]])]))
	seurat_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(seurat_object, assay = "scAB250k", slot = "data")[marker_genes[[i]],]))
}

tmp_meta <- as.matrix(seurat_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(seurat_object$banksy_cluster)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta <- tmp_meta[levels(seurat_object$banksy_cluster)[c(1:11)],]
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[1], "_hic_250k_top50_zscore_2021_Nature_E14_5.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[1], "_hic_250k_top50_zscore1_2021_Nature_E14_5.pdf"), width = 10, height = 10)

marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:50]
	count = count + 1
}
names(marker_genes) <- paste0("E14_5_", unique(de_markers$cluster))

DefaultAssay(seurat_object) <- "scAB500k"
for(i in names(marker_genes)){
	marker_genes[[i]] <- intersect(rownames(seurat_object), unique(gene_infor_500k$gene_region[which(gene_infor_500k$gene_name %in% marker_genes[[i]])]))
	seurat_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(seurat_object, assay = "scAB500k", slot = "data")[marker_genes[[i]],]))
}

tmp_meta <- as.matrix(seurat_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(seurat_object$banksy_cluster)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta <- tmp_meta[levels(seurat_object$banksy_cluster)[c(1:11)],]
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[1], "_hic_500k_top50_zscore_2021_Nature_E14_5.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[1], "_hic_500k_top50_zscore1_2021_Nature_E14_5.pdf"), width = 10, height = 10)

seurat_list$E14_5 <- seurat_object
#!---------------------------------------------------------------------------------------------------------------------------------------------------
seurat_object <- seurat_list$E16_5
locs <- Seurat::GetTissueCoordinates(seurat_object)[,seq_len(2)]
colnames(locs) <- c("array_row", "array_col")
seurat_object@meta.data <- cbind(seurat_object@meta.data, locs)

DefaultAssay(seurat_object) <- "scAB500k"
VariableFeatures(seurat_object) <- rownames(seurat_object)

seurat_object <- RunBanksy(seurat_object, lambda = 0.3, assay = 'scAB500k', slot = 'data', dimx = "array_row", dimy = "array_col", features = 'all', k_geom = 24, assay_name = "BANKSY_scAB500k")
npcs = 30
seurat_object <- RunPCA(seurat_object, assay = "BANKSY_scAB500k", reduction.name = "pca_banksy", npcs = npcs, features = rownames(seurat_object))
seurat_object <- RunUMAP(seurat_object, reduction = "pca_banksy", dims = 1:npcs, min.dist = 0.3, reduction.name = "banksy_umap", return.model = T, n.neighbors = 24)
seurat_object <- FindNeighbors(seurat_object, reduction = "pca_banksy", dims = 1:npcs, annoy.metric = "cosine", graph.name = c('banksy_nn', 'banksy_snn'), k.param = 24)
seurat_object <- FindClusters(seurat_object, resolution = seq(0.1, 2, 0.1), graph.name = "banksy_snn")
wrap_plots(map(seq(0.6, 2, 0.1), function(x) DimPlot(seurat_object, reduction = "banksy_umap", group.by = paste0("banksy_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)
#!----------------------------------------------------------------------------
seurat_object <- FindClusters(seurat_object, cluster.name = "banksy_cluster", resolution = 2, graph.name = "banksy_snn")
seurat_object <- FindSubCluster(seurat_object, cluster = 3, graph.name = "banksy_snn", subcluster.name = "banksy_sub_cluster")
seurat_object$banksy_sub_cluster[which(seurat_object$banksy_sub_cluster == "3_0")] <- 3
seurat_object$banksy_sub_cluster[which(seurat_object$banksy_sub_cluster == "3_1")] <- 21
seurat_object$banksy_cluster <- factor(seurat_object$banksy_sub_cluster, levels = 0:21)
seurat_object$banksy_cluster[which(seurat_object$banksy_snn_res.1.2 == 6)] <- 7
DimPlot(seurat_object, label = T, cols = selected_colors, group.by = "banksy_cluster", reduction = "banksy_umap") + NoLegend()
Seurat::SpatialPlot(seurat_object, cols = selected_colors, group.by = "banksy_cluster", crop = FALSE, shape = 22)
Idents(seurat_object) <- "banksy_cluster"
selected_cells <- CellsByIdentities(seurat_object)
Seurat::SpatialPlot(seurat_object, cells.highlight = selected_cells[setdiff(names(selected_cells), "NA")], cols.highlight = c("#FFFF00", alpha("lightgrey", 0)), crop = FALSE, shape = 22, facet.highlight = T, combine = T, ncol = 6, alpha = NULL) & NoLegend()

seurat_object$seurat_clusters <- seurat_object$banksy_cluster
seurat_object$banksy_cluster <- as.character(seurat_object$banksy_cluster)
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(20))] <- "C1"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(3))] <- "C2"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(21))] <- "C3"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(13, 17))] <- "C4"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(9))] <- "C5"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(8))] <- "C6"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(0))] <- "C7"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(1, 4))] <- "C8"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(11))] <- "C9"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(19))] <- "C10"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(12))] <- "C11"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(5))] <- "C12"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(15))] <- "C13"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(2))] <- "C14"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(16))] <- "C15"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(18))] <- "C16"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(7))] <- "C17"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(6, 10))] <- "C18"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(14))] <- "C19" # Stromal

seurat_object$banksy_cluster <- factor(seurat_object$banksy_cluster, levels = paste0("C", 1:19))

Idents(seurat_object) <- "banksy_cluster"

E16_colors <- c("C1" = "#4DFF99", "C2" = "#1FCCCC", "C3" = "#7ACC00", "C4" = "#1FCC1F", "C5" = "#459967", "C6" = "#8CA252", "C7" = "#00998F", "C8" = "#005300", "C9" = "#0073C2", "C10" = "#BFB1D5", "C11" = "#AA0DFE", "C12" = "#782AB6", "C13" = "#CC7A88", "C14" = "#B33E52", "C15" = "#C1FF73", "C16" = "#325A98", "C17" = "#16FF32", "C18" = "#CC79A7", "C19" = "#B3823E")
DimPlot(seurat_object, label = T, cols = E16_colors, reduction = "banksy_umap") + NoLegend() + SetAxes()
Seurat::SpatialPlot(seurat_object, cols = E16_colors, crop = FALSE, shape = 22)
seurat_object@misc$E16_colors <- E16_colors

seurat_rna <- seurat_RNA$E16_5
DefaultAssay(seurat_rna) <- "SCT"
all_markers <- FindAllMarkers(seurat_rna, only.pos = T, min.pct = 0.1)
write.csv(all_markers, "E16_5_ST_cell_identity.csv", quote = F)
de_markers <- all_markers[which(all_markers$avg_log2FC > 0.5 & all_markers$p_val_adj < 0.01 & all_markers$pct.1 > 0.3),]

marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:30]
	count = count + 1
}
names(marker_genes) <- unique(de_markers$cluster)
DefaultAssay(seurat_object) <- "scAB500k"
for(i in names(marker_genes)){
	marker_genes[[i]] <- intersect(rownames(seurat_object), unique(gene_infor_500k$gene_region[which(gene_infor_500k$gene_name %in% marker_genes[[i]])]))
	seurat_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(seurat_object, assay = "scAB500k", slot = "data")[marker_genes[[i]],]))
}

tmp_meta <- as.matrix(seurat_object@meta.data[,levels(seurat_rna)])
rownames(tmp_meta) <- paste0(seurat_object$banksy_cluster)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- paste0(levels(seurat_rna$Cluster), "_", colnames(tmp_meta))
tmp_meta0 <- tmp_meta0[levels(seurat_object$banksy_cluster),]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[2], "_hic_500k_top30_zscore.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[2], "_hic_500k_top30_zscore1.pdf"), width = 10, height = 10)

marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:30]
	count = count + 1
}
names(marker_genes) <- unique(de_markers$cluster)

DefaultAssay(seurat_object) <- "scAB250k"
for(i in names(marker_genes)){
	marker_genes[[i]] <- intersect(rownames(seurat_object), unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name %in% marker_genes[[i]])]))
	seurat_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(seurat_object, assay = "scAB250k", slot = "data")[marker_genes[[i]],]))
}

tmp_meta <- as.matrix(seurat_object@meta.data[,levels(seurat_rna)])
rownames(tmp_meta) <- paste0(seurat_object$banksy_cluster)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- paste0(levels(seurat_rna$Cluster), "_", colnames(tmp_meta))
tmp_meta0 <- tmp_meta0[levels(seurat_object$banksy_cluster),]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[2], "_hic_250k_top30_zscore.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[2], "_hic_250k_top30_zscore1.pdf"), width = 10, height = 10)

seurat_rna <- readRDS("output/E16_brain_mouse_scRNA_2021_Nature.rds")
selected_clusters <- paste0("E16_5_", levels(seurat_rna$cell_identity))

all_markers <- FindAllMarkers(seurat_rna, only.pos = T, min.pct = 0.1)
write.csv(all_markers, "E16_5_single_cells_cell_identity_2021_Nature.csv", quote = F)
de_markers <- all_markers[which(all_markers$avg_log2FC > 0.5 & all_markers$p_val_adj < 0.01 & all_markers$pct.1 > 0.3),]

marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:30]
	count = count + 1
}
names(marker_genes) <- paste0("E16_5_", unique(de_markers$cluster))

DefaultAssay(seurat_object) <- "scAB500k"
for(i in names(marker_genes)){
	marker_genes[[i]] <- intersect(rownames(seurat_object), unique(gene_infor_500k$gene_region[which(gene_infor_500k$gene_name %in% marker_genes[[i]])]))
	seurat_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(seurat_object, assay = "scAB500k", slot = "data")[marker_genes[[i]],]))
}

tmp_meta <- as.matrix(seurat_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(seurat_object$banksy_cluster)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta <- tmp_meta[levels(seurat_object$banksy_cluster)[c(1:12)],]
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[2], "_hic_500k_top30_zscore_2021_Nature_E16_5.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[2], "_hic_500k_top30_zscore1_2021_Nature_E16_5.pdf"), width = 10, height = 10)

marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:30]
	count = count + 1
}
names(marker_genes) <- paste0("E16_5_", unique(de_markers$cluster))

DefaultAssay(seurat_object) <- "scAB250k"
for(i in names(marker_genes)){
	marker_genes[[i]] <- intersect(rownames(seurat_object), unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name %in% marker_genes[[i]])]))
	seurat_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(seurat_object, assay = "scAB250k", slot = "data")[marker_genes[[i]],]))
}

tmp_meta <- as.matrix(seurat_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(seurat_object$banksy_cluster)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta <- tmp_meta[levels(seurat_object$banksy_cluster)[c(1:12)],]
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[2], "_hic_250k_top30_zscore_2021_Nature_E16_5.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[2], "_hic_250k_top30_zscore1_2021_Nature_E16_5.pdf"), width = 10, height = 10)

seurat_list$E16_5 <- seurat_object
#!---------------------------------------------------------------------------------------------------------------------------------------------------
seurat_object <- seurat_list$E18_5
locs <- Seurat::GetTissueCoordinates(seurat_object)[,seq_len(2)]
colnames(locs) <- c("array_row", "array_col")
seurat_object@meta.data <- cbind(seurat_object@meta.data, locs)

DefaultAssay(seurat_object) <- "scAB500k"
VariableFeatures(seurat_object) <- rownames(seurat_object)

seurat_object <- RunBanksy(seurat_object, lambda = 0.3, assay = 'scAB500k', slot = 'data', dimx = "array_row", dimy = "array_col", features = 'all', k_geom = 15, assay_name = "BANKSY_scAB500k")
npcs = 30
seurat_object <- RunPCA(seurat_object, assay = "BANKSY_scAB500k", reduction.name = "pca_banksy", npcs = npcs, features = rownames(seurat_object))
seurat_object <- RunUMAP(seurat_object, reduction = "pca_banksy", dims = 1:npcs, min.dist = 0.3, reduction.name = "banksy_umap", return.model = T, n.neighbors = 24)
seurat_object <- FindNeighbors(seurat_object, reduction = "pca_banksy", dims = 1:npcs, annoy.metric = "cosine", graph.name = c('banksy_nn', 'banksy_snn'), k.param = 24)
seurat_object <- FindClusters(seurat_object, resolution = seq(0.1, 2, 0.1), graph.name = "banksy_snn")
wrap_plots(map(seq(0.6, 2, 0.1), function(x) DimPlot(seurat_object, reduction = "banksy_umap", group.by = paste0("banksy_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)
#!----------------------------------------------------------------------------
seurat_object <- FindClusters(seurat_object, cluster.name = "banksy_cluster", resolution = 1.5, graph.name = "banksy_snn")
seurat_object$banksy_cluster <- as.character(seurat_object$banksy_cluster)
seurat_object$banksy_cluster[which(seurat_object$banksy_snn_res.1.8 == 11)] <- 12
seurat_object$banksy_cluster <- factor(seurat_object$banksy_cluster, levels = 0:12)

DimPlot(seurat_object, label = T, cols = selected_colors, group.by = "banksy_cluster", reduction = "banksy_umap") + NoLegend()
Seurat::SpatialPlot(seurat_object, cols = selected_colors, group.by = "banksy_cluster", crop = FALSE, shape = 22)
Idents(seurat_object) <- "banksy_cluster"
selected_cells <- CellsByIdentities(seurat_object)
Seurat::SpatialPlot(seurat_object, cells.highlight = selected_cells[setdiff(names(selected_cells), "NA")], cols.highlight = c("#FFFF00", alpha("lightgrey", 0)), crop = FALSE, shape = 22, facet.highlight = T, combine = T, ncol = 5, alpha = NULL) & NoLegend()

seurat_object$seurat_clusters <- seurat_object$banksy_cluster
seurat_object$banksy_cluster <- as.character(seurat_object$banksy_cluster)
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(9))] <- "C1"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(7))] <- "C2"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(13))] <- "C3"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(3))] <- "C4"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(2))] <- "C5"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(10))] <- "C6"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(1))] <- "C7"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(6))] <- "C8"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(8))] <- "C9"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(12, 4))] <- "C10"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(0))] <- "C11"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(5))] <- "C12"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(11))] <- "C13"

seurat_object$banksy_cluster <- factor(seurat_object$banksy_cluster, levels = paste0("C", 1:13))

Idents(seurat_object) <- "banksy_cluster"

E18_colors <- c("C1" = "#4DFF99", "C2" = "#1FCCCC", "C3" = "#7ACC00", "C4" = "#1FCC1F", "C5" = "#5CCC89", "C6" = "#00991F", "C7" = "#00998F", "C8" = "#005300", "C9" = "#AA0DFE", "C10" = "#CC7A88", "C11" = "#B33E52", "C12" = "#CC79A7", "C13" = "#B3823E")
DimPlot(seurat_object, label = T, cols = E18_colors, reduction = "banksy_umap") + NoLegend() + SetAxes()
Seurat::SpatialPlot(seurat_object, cols = E18_colors, crop = FALSE, shape = 22)
seurat_object@misc$E18_colors <- E18_colors

seurat_rna <- seurat_RNA$E18_5

DefaultAssay(seurat_rna) <- "SCT"
all_markers <- FindAllMarkers(seurat_rna, only.pos = T, min.pct = 0.1)
write.csv(all_markers, "E18_5_ST_cell_identity.csv", quote = F)
de_markers <- all_markers[which(all_markers$avg_log2FC > 0.5 & all_markers$p_val_adj < 0.01 & all_markers$pct.1 > 0.3),]

marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:30]
	count = count + 1
}
names(marker_genes) <- unique(de_markers$cluster)
DefaultAssay(seurat_object) <- "scAB500k"
for(i in names(marker_genes)){
	marker_genes[[i]] <- intersect(rownames(seurat_object), unique(gene_infor_500k$gene_region[which(gene_infor_500k$gene_name %in% marker_genes[[i]])]))
	seurat_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(seurat_object, assay = "scAB500k", slot = "data")[marker_genes[[i]],]))
}

tmp_meta <- as.matrix(seurat_object@meta.data[,levels(seurat_rna)])
rownames(tmp_meta) <- paste0(seurat_object$banksy_cluster)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- paste0(levels(seurat_rna$Cluster), "_", colnames(tmp_meta))
tmp_meta0 <- tmp_meta0[levels(seurat_object$banksy_cluster),]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[3], "_hic_500k_top30_zscore.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[3], "_hic_500k_top30_zscore1.pdf"), width = 10, height = 10)

marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:30]
	count = count + 1
}
names(marker_genes) <- unique(de_markers$cluster)

DefaultAssay(seurat_object) <- "scAB250k"
for(i in names(marker_genes)){
	marker_genes[[i]] <- intersect(rownames(seurat_object), unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name %in% marker_genes[[i]])]))
	seurat_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(seurat_object, assay = "scAB250k", slot = "data")[marker_genes[[i]],]))
}

tmp_meta <- as.matrix(seurat_object@meta.data[,levels(seurat_rna)])
rownames(tmp_meta) <- paste0(seurat_object$banksy_cluster)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- paste0(levels(seurat_rna$Cluster), "_", colnames(tmp_meta))
tmp_meta0 <- tmp_meta0[levels(seurat_object$banksy_cluster),]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[3], "_hic_250k_top30_zscore.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[3], "_hic_250k_top30_zscore1.pdf"), width = 10, height = 10)

seurat_rna <- readRDS("output/E18_brain_mouse_scRNA_2021_Nature.rds")
Idents(seurat_rna) <- "cell_identity"
selected_clusters <- paste0("E18_5_", levels(seurat_rna$cell_identity))

all_markers <- FindAllMarkers(seurat_rna, only.pos = T, min.pct = 0.1)
write.csv(all_markers, "E18_5_single_cells_cell_identity_2021_Nature.csv", quote = F)
de_markers <- all_markers[which(all_markers$avg_log2FC > 0.5 & all_markers$p_val_adj < 0.01 & all_markers$pct.1 > 0.3),]

marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:30]
	count = count + 1
}
names(marker_genes) <- paste0("E18_5_", unique(de_markers$cluster))

DefaultAssay(seurat_object) <- "scAB500k"
for(i in names(marker_genes)){
	marker_genes[[i]] <- intersect(rownames(seurat_object), unique(gene_infor_500k$gene_region[which(gene_infor_500k$gene_name %in% marker_genes[[i]])]))
	seurat_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(seurat_object, assay = "scAB500k", slot = "data")[marker_genes[[i]],]))
}

tmp_meta <- as.matrix(seurat_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(seurat_object$banksy_cluster)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta <- tmp_meta[levels(seurat_object$banksy_cluster)[c(1:9)],]
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[3], "_hic_500k_top30_zscore_2021_Nature_E18_5.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[3], "_hic_500k_top30_zscore1_2021_Nature_E18_5.pdf"), width = 10, height = 10)

marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:30]
	count = count + 1
}
names(marker_genes) <- paste0("E18_5_", unique(de_markers$cluster))

DefaultAssay(seurat_object) <- "scAB250k"
for(i in names(marker_genes)){
	marker_genes[[i]] <- intersect(rownames(seurat_object), unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name %in% marker_genes[[i]])]))
	seurat_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(seurat_object, assay = "scAB250k", slot = "data")[marker_genes[[i]],]))
}

tmp_meta <- as.matrix(seurat_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(seurat_object$banksy_cluster)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta <- tmp_meta[levels(seurat_object$banksy_cluster)[c(1:9)],]
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[3], "_hic_250k_top30_zscore_2021_Nature_E18_5.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_", names(seurat_list)[3], "_hic_250k_top30_zscore1_2021_Nature_E18_5.pdf"), width = 10, height = 10)

seurat_list$E18_5 <- seurat_object
#!---------------------------------------------------------------------------------------------------------------------------------------------------
seurat_list <- qread("output/seurat_hic_devp_list.qs")
seurat_object <- merge(seurat_list[[1]], seurat_list[2:3], add.cell.ids = names(seurat_list))
seurat_object <- JoinLayers(seurat_object, assay = "scAB500k")
seurat_object <- JoinLayers(seurat_object, assay = "scAB250k")

locs <- c()
for(i in names(seurat_object@images)){
	locs <- rbind(locs, Seurat::GetTissueCoordinates(seurat_object, image = i)[,seq_len(2)])
}
all(rownames(locs) == colnames(seurat_object))
colnames(locs) <- c("sdimy", "sdimx")
seurat_object@meta.data <- cbind(seurat_object@meta.data, locs)
seurat_object$orig.ident <- factor(seurat_object$orig.ident, levels = unique(seurat_object$orig.ident))

library(SeuratWrappers)
library(Banksy)
library(harmony)

dir.create("banksy_cluster_hic_mouse_brain_Spatial")
setwd("banksy_cluster_hic_mouse_brain_Spatial")

for(lambda in seq(0.1, 0.8, 0.05)){
	for(k in c(8, 15, 24)){
		#k = 24
		#lambda = 0.2
		resolution <- 1.5
		DefaultAssay(seurat_object) <- "scAB500k"
		tmp_object <- RunBanksy(seurat_object, lambda = lambda, assay = 'scAB500k', slot = 'data', dimx = "sdimx", dimy = "sdimy", features = 'all', group = 'orig.ident', split.scale = TRUE, k_geom = k)
		DefaultAssay(tmp_object) <- "BANKSY"
		tmp_object <- RunPCA(tmp_object, assay = "BANKSY", reduction.name = "pca", features = rownames(tmp_object), npcs = 50)
		tmp_object <- RunHarmony(tmp_object, group.by.vars = "orig.ident",  reduction.name = 'pca', reduction.save = 'harmony')

		tmp_object <- RunUMAP(tmp_object, reduction = "harmony", dims = 1:50, min.dist = 0.3, reduction.name = "banksy_umap", return.model = T, n.neighbors = 24)

		tmp_object <- FindNeighbors(tmp_object, dims = 1:50, reduction = 'harmony', annoy.metric = "cosine", graph.name = c('banksy_nn', 'banksy_snn'), k.param = 24)
		tmp_object <- FindClusters(tmp_object, resolution = 2, cluster.name = "banksy_harmony_cluster", graph.name = 'banksy_snn')
		p1 <- DimPlot(tmp_object, reduction = "banksy_umap", group.by = "banksy_harmony_cluster", cols = selected_colors, label = T) + NoLegend()
		p2 <- DimPlot(tmp_object, reduction = "banksy_umap", group.by = "orig.ident", label = F, cols = cluster_colors)
		print(p1 + p2)
		ggsave(paste0("banksy_lambda_", lambda, "_k_", k, "_banksy_umap_scAB500k.png"), width = 14.5, height = 7)

		p0 <- Seurat::SpatialPlot(tmp_object, crop = FALSE, shape = 22, group.by = "banksy_harmony_cluster", cols = selected_colors) & NoLegend()
		print(p0)
		ggsave(paste0("banksy_cluster_lambda_", lambda, "_k_", k, "_scAB500k.png"), width = 21, height = 7)
	}
}

DefaultAssay(seurat_object) <- "scAB500k"
VariableFeatures(seurat_object) <- rownames(seurat_object)
seurat_object <- RunBanksy(seurat_object, lambda = 0.3, assay = 'scAB500k', slot = 'data', dimx = "sdimx", dimy = "sdimy", features = 'all', group = 'orig.ident', split.scale = TRUE, k_geom = 24, assay_name = "BANKSY")
DefaultAssay(seurat_object) <- "BANKSY"
seurat_object <- RunPCA(seurat_object, assay = "BANKSY", reduction.name = "pca", features = rownames(seurat_object), npcs = 50)
seurat_object <- RunHarmony(seurat_object, group.by.vars = "orig.ident",  reduction.name = 'pca', reduction.save = 'harmony')

seurat_object <- RunUMAP(seurat_object, reduction = "harmony", dims = 1:50, min.dist = 0.3, reduction.name = "banksy_umap", return.model = T, n.neighbors = 24)
seurat_object <- FindNeighbors(seurat_object, dims = 1:50, reduction = 'harmony', annoy.metric = "cosine", graph.name = c('banksy_nn', 'banksy_snn'), k.param = 24)
seurat_object <- FindClusters(seurat_object, resolution = seq(0.1, 2, 0.1), graph.name = "banksy_snn")

wrap_plots(map(seq(0.6, 2, 0.1), function(x) DimPlot(seurat_object, reduction = "banksy_umap", group.by = paste0("banksy_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)
seurat_object <- FindClusters(seurat_object, cluster.name = "banksy_cluster", resolution = 1.5, graph.name = 'banksy_snn')

DimPlot(seurat_object, label = T, cols = selected_colors, group.by = "banksy_cluster", reduction = "banksy_umap") + NoLegend()
Seurat::SpatialPlot(seurat_object, crop = FALSE, shape = 22, cols = selected_colors) & NoLegend()

Idents(seurat_object) <- "banksy_cluster"
selected_cells <- CellsByIdentities(seurat_object)
Seurat::SpatialPlot(seurat_object, cells.highlight = selected_cells[setdiff(names(selected_cells), "NA")], cols.highlight = c("#FFFF00", alpha("lightgrey", 0)), crop = FALSE, shape = 22, facet.highlight = T, combine = T, ncol = 5, alpha = NULL) & NoLegend()

DefaultAssay(seurat_object) <- "scAB500k"
library(purrr)
wrap_plots(map(seq(0.1, 1, 0.1), function(x) DimPlot(seurat_object, reduction = "banksy_umap", group.by = paste0("banksy_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)

seurat_object$seurat_clusters <- seurat_object$banksy_cluster

seurat_object$banksy_cluster <- NA
seurat_object$banksy_cluster[match(paste0("E14_5_", colnames(seurat_list$E14_5)), colnames(seurat_object))] <- paste0("E14_5_", as.character(seurat_list$E14_5$banksy_cluster))
seurat_object$banksy_cluster[match(paste0("E16_5_", colnames(seurat_list$E16_5)), colnames(seurat_object))] <- paste0("E16_5_", as.character(seurat_list$E16_5$banksy_cluster))
seurat_object$banksy_cluster[match(paste0("E18_5_", colnames(seurat_list$E18_5)), colnames(seurat_object))] <- paste0("E18_5_", as.character(seurat_list$E18_5$banksy_cluster))
names(E14_colors) <- paste0("E14_5_", names(E14_colors))
names(E16_colors) <- paste0("E16_5_", names(E16_colors))
names(E18_colors) <- paste0("E18_5_", names(E18_colors))
seurat_object@misc$banksy_cluster_colors <- c(E14_colors, E16_colors, E18_colors)

E14_embed <- seurat_list$E14_5@reductions$banksy_umap@cell.embeddings
rownames(E14_embed) <- paste0("E14_5_", rownames(E14_embed))
E16_embed <- seurat_list$E16_5@reductions$banksy_umap@cell.embeddings
rownames(E16_embed) <- paste0("E16_5_", rownames(E16_embed))
E18_embed <- seurat_list$E18_5@reductions$banksy_umap@cell.embeddings
rownames(E18_embed) <- paste0("E18_5_", rownames(E18_embed))

E16_embed[,1] <- E16_embed[,1] + max(E14_embed[,1]) + 20
E18_embed[,1] <- E18_embed[,1] + max(E16_embed[,1]) + 20
tmp_embed <- rbind(E14_embed, E16_embed, E18_embed)
colnames(tmp_embed) <- c("banksyperumap_1", "banksyperumap_2")
all(rownames(tmp_embed) == colnames(seurat_object))
seurat_object[["banksy_per_umap"]] <- CreateDimReducObject(embeddings = tmp_embed, assay = "scAB500k")

seurat_object$staggered_sdimy <- 1080 - seurat_object$staggered_sdimy
tmp <- seurat_object@meta.data[,c("staggered_sdimx", "staggered_sdimy")]
colnames(tmp) <- c("spatial_umap_1", "spatial_umap_2")
seurat_object[["spatial_umap"]] <- CreateDimReducObject(embeddings = as.matrix(tmp), assay = "scAB250k")

Idents(seurat_object) <- "banksy_cluster"
DimPlot(seurat_object, reduction = "banksy_per_umap", cols = seurat_object@misc$banksy_cluster_colors)
qsave(seurat_object, "output/brain_development_Spatial_HiC_final_tutorial.qs")
#!----------------------------------------------------------------------------------------------------------------
seurat2scanpy(seurat_list$E14_5, assay = "scAB250k", manual_color = seurat_list$E14_5@misc$E14_colors, major_umap = "banksy_umap", savefile = paste0("Trajectory_analysis/anndata/brain_development_hic_E14_5_final.h5ad"))
seurat2scanpy(seurat_list$E16_5, assay = "scAB250k", manual_color = seurat_list$E16_5@misc$E16_colors, major_umap = "banksy_umap", savefile = paste0("Trajectory_analysis/anndata/brain_development_hic_E16_5_final.h5ad"))
seurat2scanpy(seurat_list$E18_5, assay = "scAB250k", manual_color = seurat_list$E18_5@misc$E18_colors, major_umap = "banksy_umap", savefile = paste0("Trajectory_analysis/anndata/brain_development_hic_E18_5_final.h5ad"))
seurat2scanpy(seurat_object, assay = "scAB250k", manual_color = seurat_object@misc$banksy_cluster_colors, major_umap = "banksy_umap", savefile = paste0("Trajectory_analysis/anndata/brain_development_hic_final.h5ad"))

seurat_object$Region <- as.character(seurat_object$banksy_cluster)
source("~/scripts/seurat2scanpy/shiny_st.R")
options(browser = "/usr/bin/firefox")
seurat_object <- shiny_st(seurat = seurat_object, isVisium = F, assay = "scAB250k", image = "E14_5")
seurat_object$Region <- factor(seurat_object$Region, levels = c("Cortex", "MGE", "LGE", "CPU", "LCS", "Epd", "PIR", "CLA", "LS", "Others"))
Idents(seurat_object) <- "Region"
VlnPlot(seurat_object, features = "Short_Long_Ratio", idents = c("Cortex", "LGE", "CPU", "LS", "Others"), pt.size = 0, split.by = "orig.ident", cols = c("#34D916", "#00D4E6", "#1E90FF")) & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.3, position = position_dodge(0.9), show.legend = F)

DefaultAssay(seurat_object) <- "dist"
key_markers <- c("20kb", "103.75kb", "493.51kb", "1.97Mb", "10.24Mb", "48.71Mb")

plot_doheatmap(dataset = seurat_object, slot = "data", markers = rownames(seurat_object), sort_var = "Region", anno_var = c("Region", "banksy_cluster", "orig.ident"), hm_limit = seq(2, 6, length = 256), hm_colors = paletteContinuous(set = "solarExtra"), anno_colors = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))], seurat_object@misc$banksy_cluster_colors, c("#34D916", "#00D4E6", "#1E90FF")), column_split = "Region",
				row_font_size = 10, label_markers = key_markers, cor_row_label_line = "blue", lwd_row_label_line = 0.3)	
#!----------------------------------------------------------------------------------------------------------------
temp_object <- seurat_object
seurat_object <- subset(seurat_object, subset = banksy_cluster %in% c(paste0("E14_5_C", 1:7), paste0("E16_5_C", 1:9), paste0("E18_5_C", 1:8)))
seurat_object <- shiny_st(seurat = seurat_object, isVisium = F, assay = "scAB250k", image = "E14_5")
seurat_object <- subset(seurat_object, subset = banksy_cluster != "removed_cells")
seurat_object$banksy_cluster <- factor(seurat_object$banksy_cluster, levels = c(paste0("E14_5_C", 1:7), paste0("E16_5_C", 1:9), paste0("E18_5_C", 1:8)))
seurat_object@meta.data <- seurat_object@meta.data[,c(1:7, 29)]
Idents(seurat_object) <- "banksy_cluster"

seurat_object$Region <- as.character(seurat_object$banksy_cluster)
seurat_object <- shiny_st(seurat = seurat_object, isVisium = F, assay = "scAB250k", image = "E14_5")
seurat_object$Region <- factor(seurat_object$Region, levels = c("VZ0", "VZ", "SVZ", "IZ", "CP", "Layer_VI", "Layer_V", "Layer_II_III_IV"))

seurat_object$banksy_enrich <- as.character(seurat_object$banksy_cluster)
seurat_object$banksy_enrich[which(seurat_object$banksy_enrich %in% c("E14_5_C1"))] <- "AP0"
seurat_object$banksy_enrich[which(seurat_object$banksy_enrich %in% c("E14_5_C2", "E16_5_C1", "E16_5_C2", "E18_5_C1", "E18_5_C2"))] <- "AP"
seurat_object$banksy_enrich[which(seurat_object$banksy_enrich %in% c("E14_5_C4", "E14_5_C5", "E16_5_C3", "E16_5_C4", "E16_5_C5", "E16_5_C6", "E18_5_C3", "E18_5_C4"))] <- "IPC_MigN"
seurat_object$banksy_enrich[which(seurat_object$banksy_enrich %in% c("E14_5_C6", "E16_5_C7", "E16_5_C8"))] <- "CThPN_SCPN"
seurat_object$banksy_enrich[which(seurat_object$banksy_enrich %in% c( "E18_5_C5", "E18_5_C6"))] <- "CThPN"
seurat_object$banksy_enrich[which(seurat_object$banksy_enrich %in% c( "E18_5_C7", "E18_5_C8"))] <- "SCPN_UL_CPN"
seurat_object$banksy_enrich[which(seurat_object$banksy_enrich %in% c( "C2"))] <- "Layer_II_III"
seurat_object$banksy_enrich[which(seurat_object$banksy_enrich %in% c( "C3"))] <- "Layer_IV"
seurat_object$banksy_enrich[which(seurat_object$banksy_enrich %in% c( "C4"))] <- "Layer_V"
seurat_object$banksy_enrich[which(seurat_object$banksy_enrich %in% c( "C5"))] <- "Layer_VI"
seurat_object$banksy_enrich <- factor(seurat_object$banksy_enrich, levels = c("AP0", "AP", "IPC_MigN", "CThPN_SCPN", "CThPN", "SCPN_UL_CPN", "Layer_VI", "Layer_V", "Layer_IV", "Layer_II_III"))
qsave(seurat_object, "output/brain_Spatial_HiC_final_tutorial_cortex_new.qs")
seurat2scanpy(seurat_object, assay = "scAB250k", manual_color = list(banksy_cluster_colors = seurat_object@misc$banksy_cluster_colors, cluster_colors = cluster_colors), major_umap = "banksy_umap", savefile = paste0("Trajectory_analysis/anndata/brain_development_hic_final_cortex_new.h5ad"))
#!------------------------------------------------------------------------------------------------------------------
seurat_object <- temp_object
seurat_object <- subset(seurat_object, subset = banksy_cluster %in% c(paste0("E14_5_C", c(9, 11, 12, 13)), paste0("E16_5_C", 11:14), paste0("E18_5_C", 9:11)))
seurat_object <- shiny_st(seurat = seurat_object, isVisium = F, assay = "scAB250k", image = "E14_5")
seurat_object <- subset(seurat_object, subset = banksy_cluster != "removed_cells")
seurat_object$banksy_cluster <- factor(seurat_object$banksy_cluster, levels = c(paste0("E14_5_C", c(9, 11, 12, 13)), paste0("E16_5_C", 11:14), paste0("E18_5_C", 9:11)))
seurat_object@meta.data <- seurat_object@meta.data[,c(1:7, 29)]
Idents(seurat_object) <- "banksy_cluster"

qsave(seurat_object, "output/brain_development_Spatial_HiC_final_tutorial_MSN.qs")
seurat2scanpy(seurat_object, assay = "scAB250k", manual_color = seurat_object@misc$banksy_cluster_colors, major_umap = "banksy_umap", savefile = paste0("Trajectory_analysis/anndata/brain_development_hic_final_MSN.h5ad"))
#!---------------------------------------------------------------------------------------------------------------------------------------------------
seurat_object <- qread("output/brain_Adult16_HiC_final_tutorial.qs")
locs <- Seurat::GetTissueCoordinates(seurat_object)[,seq_len(2)]
colnames(locs) <- c("array_row", "array_col")
seurat_object@meta.data <- cbind(seurat_object@meta.data, locs)

DefaultAssay(seurat_object) <- "scAB250k"
VariableFeatures(seurat_object) <- rownames(seurat_object)

for(lambda in seq(0.1, 0.4, 0.05)){
	for(k in c(8, 15, 24)){
		#k = 24
		#lambda = 0.2
		resolution <- 1.5
		DefaultAssay(seurat_object) <- "scAB250k"
		tmp_object <- RunBanksy(seurat_object, lambda = lambda, assay = 'scAB250k', slot = 'data', dimx = "array_row", dimy = "array_col", features = 'all', group = 'orig.ident', split.scale = TRUE, k_geom = k)
		DefaultAssay(tmp_object) <- "BANKSY"
		tmp_object <- RunPCA(tmp_object, assay = "BANKSY", reduction.name = "pca_banksy", npcs = npcs, features = rownames(tmp_object))
		tmp_object <- RunUMAP(tmp_object, reduction = "pca_banksy", dims = 1:npcs, min.dist = 0.1, reduction.name = "banksy_umap", return.model = T, n.neighbors = 24)
		tmp_object <- FindNeighbors(tmp_object, reduction = "pca_banksy", dims = 1:npcs, annoy.metric = "cosine", graph.name = c('banksy_nn', 'banksy_snn'), k.param = 24)
		tmp_object <- FindClusters(tmp_object, resolution = 1.5, cluster.name = "banksy_cluster", graph.name = 'banksy_snn')
		p0 <- DimPlot(tmp_object, reduction = "banksy_umap", group.by = "banksy_cluster", cols = selected_colors, label = T) + NoLegend()
		p1 <- Seurat::SpatialPlot(tmp_object, crop = FALSE, shape = 22, group.by = "banksy_cluster", cols = selected_colors) & NoLegend()
		print(p0 + p1)
		ggsave(paste0("Adult16_banksy_cluster_lambda_", lambda, "_k_", k, "_scAB250k.png"), width = 14, height = 7)
		Idents(tmp_object) <- "banksy_cluster"
		selected_cells <- CellsByIdentities(tmp_object)
		p0 <- Seurat::SpatialPlot(tmp_object, cells.highlight = selected_cells[setdiff(names(selected_cells), "NA")], cols.highlight = c("#FFFF00", alpha("lightgrey", 0)), crop = FALSE, shape = 22, facet.highlight = T, combine = T, ncol = 7, alpha = NULL) & NoLegend()
		print(p0)
		ggsave(paste0("Adult16_banksy_per_cluster_lambda_", lambda, "_k_", k, "_scAB250k.png"), width = 28, height = 20)
		
	}
}

npcs = 30
seurat_object <- RunBanksy(seurat_object, lambda = 0.4, assay = 'scAB500k', slot = 'data', dimx = "array_row", dimy = "array_col", features = 'all', k_geom = 24, assay_name = "BANKSY_scAB500k")
seurat_object <- RunPCA(seurat_object, assay = "BANKSY_scAB500k", reduction.name = "pca_banksy", npcs = 30, features = rownames(seurat_object))
seurat_object <- RunUMAP(seurat_object, reduction = "pca_banksy", dims = 1:npcs, min.dist = 0.1, reduction.name = "banksy_umap", return.model = T, n.neighbors = 24)
seurat_object <- FindNeighbors(seurat_object, reduction = "pca_banksy", dims = 1:npcs, annoy.metric = "cosine", graph.name = c('banksy_nn', 'banksy_snn'), k.param = 24)
seurat_object <- FindClusters(seurat_object, resolution = seq(0.1, 2, 0.1), graph.name = "banksy_snn")
wrap_plots(map(seq(0.6, 2, 0.1), function(x) DimPlot(seurat_object, reduction = "banksy_umap", group.by = paste0("banksy_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)
#!----------------------------------------------------------------------------
seurat_object <- FindClusters(seurat_object, cluster.name = "banksy_cluster", resolution = 1.1, graph.name = "banksy_snn")
DimPlot(seurat_object, label = T, cols = selected_colors, group.by = "banksy_cluster", reduction = "banksy_umap") + NoLegend()
Seurat::SpatialPlot(seurat_object, cols = selected_colors, group.by = "banksy_cluster", crop = FALSE, shape = 22)
Idents(seurat_object) <- "banksy_cluster"
selected_cells <- CellsByIdentities(seurat_object)
Seurat::SpatialPlot(seurat_object, cells.highlight = selected_cells[setdiff(names(selected_cells), "NA")], cols.highlight = c("#FFFF00", alpha("lightgrey", 0)), crop = FALSE, shape = 22, facet.highlight = T, combine = T, ncol = 5, alpha = NULL) & NoLegend()

tmp_object <- subset(seurat_object, subset = banksy_cluster %in% c(0, 1, 3, 4, 6, 11))
DefaultAssay(tmp_object) <- "scAB500k"
tmp_object <- RunBanksy(tmp_object, lambda = 0.4, assay = 'scAB500k', slot = 'data', dimx = "array_row", dimy = "array_col", features = 'all', k_geom = 24, assay_name = "BANKSY_scAB500k")
tmp_object <- RunPCA(tmp_object, assay = "BANKSY_scAB500k", reduction.name = "pca_banksy", npcs = 30, features = rownames(tmp_object))
tmp_object <- RunUMAP(tmp_object, reduction = "pca_banksy", dims = 1:npcs, min.dist = 0.3, reduction.name = "banksy_umap", return.model = T, n.neighbors = 24)
tmp_object <- FindNeighbors(tmp_object, reduction = "pca_banksy", dims = 1:npcs, annoy.metric = "cosine", graph.name = c('banksy_nn', 'banksy_snn'), k.param = 24)
tmp_object <- FindClusters(tmp_object, cluster.name = "banksy_cluster", resolution = 0.1, graph.name = "banksy_snn")

DimPlot(tmp_object, label = T, cols = selected_colors, group.by = "banksy_cluster", reduction = "banksy_umap") + NoLegend()
Seurat::SpatialPlot(tmp_object, cols = selected_colors, group.by = "banksy_cluster", crop = FALSE, shape = 22)

seurat_object$banksy_cluster <- as.character(seurat_object$banksy_cluster)
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster == 2)] <- 4
seurat_object$banksy_cluster[match(colnames(tmp_object), colnames(seurat_object))] <- as.character(tmp_object$banksy_cluster)
seurat_object$banksy_cluster <- factor(seurat_object$banksy_cluster, levels = c(0:5, 7:10, 12:14))

DimPlot(seurat_object, label = T, cols = selected_colors, group.by = "banksy_cluster", reduction = "banksy_umap") + NoLegend()
Seurat::SpatialPlot(seurat_object, cols = selected_colors, group.by = "banksy_cluster", crop = FALSE, shape = 22)

Idents(seurat_object) <- "banksy_cluster"
selected_cells <- CellsByIdentities(seurat_object)
Seurat::SpatialPlot(seurat_object, cells.highlight = selected_cells[setdiff(names(selected_cells), "NA")], cols.highlight = c("#FFFF00", alpha("lightgrey", 0)), crop = FALSE, shape = 22, facet.highlight = T, combine = T, ncol = 5, alpha = NULL) & NoLegend()

seurat_object$seurat_clusters <- seurat_object$banksy_cluster
seurat_object$banksy_cluster <- as.character(seurat_object$banksy_cluster)
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(3))] <- "C1"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(2))] <- "C2"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(4, 10))] <- "C3"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(1))] <- "C4"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(0))] <- "C5"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(9, 14))] <- "C6"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(5))] <- "C7"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(12))] <- "C8"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(13))] <- "C9"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(7))] <- "C10"
seurat_object$banksy_cluster[which(seurat_object$banksy_cluster %in% c(8))] <- "C11"
seurat_object$banksy_cluster <- factor(seurat_object$banksy_cluster, levels = paste0("C", 1:11))

Idents(seurat_object) <- "banksy_cluster"

Adult16_colors <- c("C1" = "#3DCCB1", "C2" = "#017351", "C3" = "#54990F", "C4" = "#A3CC7A", "C5" = "#66C5CC", "C6" = "#64C2FC", "C7" = "#D38B5C", "C8" = "#D4FF26", "C9" = "#9EFF99", "C10" = "#886C00", "C11" = "#dde2e6")
DimPlot(seurat_object, label = T, cols = Adult16_colors, reduction = "banksy_umap") + NoLegend() + SetAxes()
Seurat::SpatialPlot(seurat_object, cols = Adult16_colors, crop = FALSE, shape = 22)
seurat_object@misc$Adult16_colors <- Adult16_colors

seurat_rna <- qread("output/brain_Adult15_ST_final_tutorial.qs")

DefaultAssay(seurat_rna) <- "SCT"
Idents(seurat_rna) <- "banksy_domains"
all_markers <- FindAllMarkers(seurat_rna, only.pos = T, min.pct = 0.1)
write.csv(all_markers, "Adult15_ST_banksy_domains.csv", quote = F)
de_markers <- all_markers[which(all_markers$avg_log2FC > 0.5 & all_markers$p_val_adj < 0.01 & all_markers$pct.1 > 0.3),]

marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:30]
	count = count + 1
}
names(marker_genes) <- unique(de_markers$cluster)

DefaultAssay(seurat_object) <- "scAB500k"
for(i in names(marker_genes)){
	marker_genes[[i]] <- intersect(rownames(seurat_object), unique(gene_infor_500k$gene_region[which(gene_infor_500k$gene_name %in% marker_genes[[i]])]))
	seurat_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(seurat_object, assay = "scAB500k", slot = "data")[marker_genes[[i]],]))
}

tmp_meta <- as.matrix(seurat_object@meta.data[,levels(seurat_rna)])
rownames(tmp_meta) <- paste0(seurat_object$banksy_cluster)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- paste0(levels(seurat_rna$banksy_Cluster), "_", colnames(tmp_meta))
tmp_meta0 <- tmp_meta0[levels(seurat_object$banksy_cluster),]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_Adult16_hic_500k_top30_zscore.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_Adult16_hic_500k_top30_zscore1.pdf"), width = 10, height = 10)

marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:30]
	count = count + 1
}
names(marker_genes) <- unique(de_markers$cluster)

DefaultAssay(seurat_object) <- "scAB250k"
for(i in names(marker_genes)){
	marker_genes[[i]] <- intersect(rownames(seurat_object), unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name %in% marker_genes[[i]])]))
	seurat_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(seurat_object, assay = "scAB250k", slot = "data")[marker_genes[[i]],]))
}

tmp_meta <- as.matrix(seurat_object@meta.data[,levels(seurat_rna)])
rownames(tmp_meta) <- paste0(seurat_object$banksy_cluster)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- paste0(levels(seurat_rna$banksy_Cluster), "_", colnames(tmp_meta))
tmp_meta0 <- tmp_meta0[levels(seurat_object$banksy_cluster),]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_Adult16_hic_250k_top30_zscore.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_Adult16_hic_250k_top30_zscore1.pdf"), width = 10, height = 10)

seurat_rna <- readRDS("output/FC_mouse_snRNA_2023_NSMB.rds")
Idents(seurat_rna) <- "cell_identity"
selected_clusters <- paste0("2023_NSMB_", levels(seurat_rna$cell_identity))

all_markers <- FindAllMarkers(seurat_rna, only.pos = T, min.pct = 0.1)
write.csv(all_markers, "single_cells_cell_identity_2023_NSMB.csv", quote = F)
de_markers <- all_markers[which(all_markers$avg_log2FC > 0.5 & all_markers$p_val_adj < 0.01 & all_markers$pct.1 > 0.3),]

marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:50]
	count = count + 1
}
names(marker_genes) <- paste0("2023_NSMB_", unique(de_markers$cluster))

DefaultAssay(seurat_object) <- "scAB500k"
for(i in names(marker_genes)){
	marker_genes[[i]] <- intersect(rownames(seurat_object), unique(gene_infor_500k$gene_region[which(gene_infor_500k$gene_name %in% marker_genes[[i]])]))
	seurat_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(seurat_object, assay = "scAB500k", slot = "data")[marker_genes[[i]],]))
}

tmp_meta <- as.matrix(seurat_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(seurat_object$banksy_cluster)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)
tmp_meta0 <- tmp_meta0[levels(seurat_object$banksy_cluster),]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_Adult16_hic_500k_top50_zscore_2023_NSMB.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_Adult16_hic_500k_top50_zscore1_2023_NSMB.pdf"), width = 10, height = 10)

marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:50]
	count = count + 1
}
names(marker_genes) <- paste0("2023_NSMB_", unique(de_markers$cluster))

DefaultAssay(seurat_object) <- "scAB250k"
for(i in names(marker_genes)){
	marker_genes[[i]] <- intersect(rownames(seurat_object), unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name %in% marker_genes[[i]])]))
	seurat_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(seurat_object, assay = "scAB250k", slot = "data")[marker_genes[[i]],]))
}

tmp_meta <- as.matrix(seurat_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(seurat_object$banksy_cluster)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)
tmp_meta0 <- tmp_meta0[levels(seurat_object$banksy_cluster),]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_Adult16_hic_250k_top50_zscore_2023_NSMB.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_Adult16_hic_250k_top50_zscore1_2023_NSMB.pdf"), width = 10, height = 10)

tmp_meta <- as.matrix(seurat_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(seurat_object$banksy_cluster)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta <- tmp_meta[levels(seurat_object$banksy_cluster)[2:5],selected_clusters[c(2:8)]]
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_Adult16_hic_250k_top50_zscore_2023_NSMB_cortex.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_2_density_Adult16_hic_250k_top50_zscore1_2023_NSMB_cortex.pdf"), width = 10, height = 10)
#!----------------------------------------------------------------------------------------------------
source("~/scripts/seurat2scanpy/shiny_st.R")
options(browser = "/usr/bin/firefox")
seurat_object$Region <- seurat_object$banksy_cluster
seurat_object <- shiny_st(seurat = seurat_object, isVisium = F, assay = "scAB250k", image = "Adult16")
seurat_object$Region <- factor(seurat_object$Region, levels = c("Layer_I", "Layer_II_III", "Layer_IV", "Layer_V", "Layer_VI", "CC", "CPU"))
seurat_object$banksy_type <- "excitory_neuron"
seurat_object$banksy_type[which(seurat_object$banksy_cluster %in% c("C6", "C9"))] <- "ODC"
seurat_object$banksy_type[which(seurat_object$banksy_cluster %in% c("C8"))] <- "Astrocytes"
seurat_object$banksy_type[which(seurat_object$banksy_cluster %in% c("C7"))] <- "MSN"
seurat_object$banksy_type[which(seurat_object$banksy_cluster %in% c("C10", "C11"))] <- "others"
seurat_object$banksy_type <- factor(seurat_object$banksy_type, levels = rev(c("excitory_neuron", "ODC", "Astrocytes", "MSN", "others")))

ggplot(seurat_object@meta.data[which(seurat_object$banksy_type %in% c("excitory_neuron", "ODC", "Astrocytes", "others") & seurat_object$Region %in% c("Layer_I", "Layer_II_III", "Layer_IV", "Layer_V", "Layer_VI")),], aes(Region, fill = banksy_type)) + 
	geom_bar(position = "fill") + ylab("percent of cells") + xlab("") + 
    theme(plot.background = element_rect(fill = NA), legend.key = element_rect(size = 2), legend.background = element_rect(color = NA), axis.text.y = element_text(size = 10, face='bold'), 
    axis.text.x = element_text(size = 10, face='bold', angle = 45, hjust = 1, vjust = 1)) +  
    scale_fill_manual(values = rev(c("#009E73", "#0072B2", "#16F2F2", "#886C00")))
ggsave("Development_QC_cell_identity/Adult16_Region_banksy_cluster_barplot.pdf")

tmp_data <- table(seurat_object$Region, seurat_object$banksy_cluster)
tmp_data <- tmp_data/rowSums(tmp_data)
meta_data <- data.frame(celltype = rep(colnames(tmp_data), each = nrow(tmp_data)), values = as.vector(tmp_data), group = rep(rownames(tmp_data), ncol(tmp_data)))
meta_data$celltype <- factor(meta_data$celltype, levels = colnames(tmp_data))
meta_data$group <- factor(meta_data$group, levels = rownames(tmp_data))
meta_data <- meta_data[which(meta_data$celltype %in% c("C1", "C2", "C3", "C4", "C5", "C6", "C8", "C9", "C10") & meta_data$group %in% c("Layer_I", "Layer_II_III", "Layer_IV", "Layer_V", "Layer_VI")),]

ggplot(meta_data, aes(x = group, y = values, fill = celltype)) + geom_bar(stat = "identity", position = "dodge") + ylab("percent of cells") + xlab("") + 
	coord_cartesian(ylim = c(0, 0.13)) +
	geom_text(aes(label = paste0(round(meta_data$values * 100, 1), "%")), position = position_dodge(0.9), vjust = -0.25) +
    theme(plot.background = element_rect(fill = NA), legend.position = "right", legend.key = element_rect(size = 2), legend.background = element_rect(color = NA), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12)) +
    scale_fill_manual(values = seurat_object@misc$Adult16_colors)

    ggplot(meta_data, aes(x = group, y = values, fill = celltype)) + geom_bar(stat = "identity", position = "dodge") + ylab("percent of cells") + xlab("") + 
	coord_cartesian(ylim = c(0.18, 0.8)) +
	geom_text(aes(label = paste0(round(meta_data$values * 100, 1), "%")), position = position_dodge(0.9), vjust = -0.25) +
    theme(plot.background = element_rect(fill = NA), legend.position = "right", legend.key = element_rect(size = 2), legend.background = element_rect(color = NA), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12)) +
    scale_fill_manual(values = seurat_object@misc$Adult16_colors)
ggsave("Development_QC_cell_identity/Adult16_cortex_Region_banksy_cluster_barplot.pdf")
#!------------------------------------------------------------------------------------------------
library(semla)
tmp_object <- qread("output/brain_development_ST_final_tutorial_cortex.qs")
seurat_list <- qread("output/seurat_list_SCT.qs")
seurat_list$E14_5 <- RenameCells(seurat_list$E14_5, add.cell.id = "E14_5")
seurat_list$E16_5 <- RenameCells(seurat_list$E16_5, add.cell.id = "E16_5")
seurat_list$E18_5 <- RenameCells(seurat_list$E18_5, add.cell.id = "E18_5")

seurat_object <- UpdateSeuratForSemla(seurat_list$E14_5)
seurat_object <- LoadImages(seurat_object, image_height = 1080)
temp <- seurat_object@meta.data[,c("col", "row")]
temp <- temp[seurat_object@tools$Staffli@meta_data$barcode,]
seurat_object@tools$Staffli@meta_data$x = temp$row
seurat_object@tools$Staffli@meta_data$y = temp$col
seurat_object@tools$Staffli@meta_data <- seurat_object@tools$Staffli@meta_data[match(colnames(seurat_object), seurat_object@tools$Staffli@meta_data$barcode),]

seurat_object$Region <- 'others'
seurat_object$Region[which(colnames(seurat_object) %in% colnames(tmp_object))] <- "cortex"
seurat_object@meta.data[which(seurat_object$Region == "others"),53:81] <- 0

MapMultipleFeatures(seurat_object, image_use = "raw", shape = "tile", pt_size = 6, max_cutoff = 0.95, min_cutoff = 0.05, override_plot_dims = TRUE, colors = c("#03FFF4", "#0BD3B1", "#62CFE8", "#D4FF26", "#7EC136","#34A047"), features = c("E14_5_Apical_progenitors", "E14_5_Intermediate_progenitors", "E14_5_Migrating_neurons", "E14_5_Immature_neurons", "E14_5_CThPN", "E14_5_SCPN"))
MapMultipleFeatures(seurat_object, image_use = "raw", shape = "tile", pt_size = 6, max_cutoff = 0.95, min_cutoff = 0.05, override_plot_dims = TRUE, colors = c("#03FFF4", "#0BD3B1", "#62CFE8", "#00CC14", "#7EC136","#34A047", "#01545a"), features = c("E16_5_Apical_progenitors", "E16_5_Intermediate_progenitors", "E16_5_Migrating_neurons", "E16_5_CThPN", "E16_5_SCPN", "E16_5_DL_CPN", "E16_5_UL_CPN"))

MapMultipleFeatures(seurat_object, image_use = "raw", shape = "tile", pt_size = 6, max_cutoff = 0.95, min_cutoff = 0.05, override_plot_dims = TRUE, colors = c("#0BD3B1", "#62CFE8", "#99FFBC", "#B6CC5C", "#00CC14", "#7EC136","#34A047", "#267dFF", "#01545a"), features = c("E18_5_Intermediate_progenitors", "E18_5_Migrating_neurons", "E18_5_Layer_6b", "NP", "E18_5_CThPN", "E18_5_SCPN", "E18_5_DL_CPN", "E18_5_Layer_4", "E18_5_UL_CPN"))
#!----------------------------
MapMultipleFeatures(seurat_object, image_use = "raw", shape = "tile", pt_size = 6, max_cutoff = 0.95, min_cutoff = 0.05, override_plot_dims = TRUE,  colors = c("#017351", "#332288", "#44AA99", "#C9DB74", "#DDCC77", "#CC6677","#AA4499", "#16F2F2"), features = c("ITL1", "ITL23GL", "ITL4GL", "ITL5GL", "PTGL", "CTGL_ITL6GL", "CPU_MSN", "Astrocytes"))
VlnPlot(seurat_object, cols = c("#017351", "#332288", "#44AA99", "#C9DB74", "#DDCC77", "#CC6677","#AA4499", "#16F2F2"), features = c("ITL23GL", "ITL45GL", "ITL5GL", "PTGL", "CTGL", "ITL6GL", "D12MSN", "ASC"), pt.size = 0, stack = T, flip = T) + ylab("cell2location score")
#!-------------------------------------------------------------------------------------------------
FeaturePlot(seurat_object, reduction = "spatial_umap", features = "Short_Long_Ratio", shape.by = "shape", pt.size = 1) & scale_shape_manual(values = 15) & scale_color_gradientn(colours = paletteContinuous(set = "sambaNight"))
#!-------------------------------------------------------------------------------------------------
E14_contacts <- read.csv("../E145-brain-hic-4_output_matrix_heatmap.csv", row.names = 1)
length(intersect(seurat_object$pixel[which(seurat_object$orig.ident == "E14_5")], colnames(E14_contacts)))
E14_contacts <- E14_contacts[rev(rownames(E14_contacts)),seurat_object$pixel[which(seurat_object$orig.ident == "E14_5")]]

E16_contacts <- read.csv("../E165-brain-hic-17_output_matrix_heatmap.csv", row.names = 1)
length(intersect(seurat_object$pixel[which(seurat_object$orig.ident == "E16_5")], colnames(E16_contacts)))
E16_contacts <- E16_contacts[rev(rownames(E16_contacts)),seurat_object$pixel[which(seurat_object$orig.ident == "E16_5")]]

E18_contacts <- read.csv("../E185-brain-hic-15_output_matrix_heatmap.csv", row.names = 1)
length(intersect(seurat_object$pixel[which(seurat_object$orig.ident == "E18_5")], colnames(E18_contacts)))
E18_contacts <- E18_contacts[rev(rownames(E18_contacts)),seurat_object$pixel[which(seurat_object$orig.ident == "E18_5")]]
#!----------------------------------------------------------------------
Adult_contacts <- read.csv("../Adult-HiC-16_output_matrix_heatmap.csv", row.names = 1)
length(intersect(seurat_object$pixel, colnames(Adult_contacts)))
Adult_contacts <- Adult_contacts[rev(rownames(Adult_contacts)),seurat_object$pixel]
#!----------------------------------------------------------------------
tmp <- as.numeric(rownames(E14_contacts))
tmp <- round(tmp / 1000, 2)
tmp1 <- paste0(tmp[64:133], "kb")
tmp <- round(tmp / 1000, 2)
tmp2 <- paste0(tmp[1:63], "Mb")
tmp <- c(tmp2, tmp1)
megacontacts <- cbind(E14_contacts, E16_contacts, E18_contacts)
rownames(megacontacts) <- tmp
all(colnames(megacontacts) == seurat_object$pixel)
colnames(megacontacts) <- colnames(seurat_object)

seurat_object[["dist"]] <- CreateAssay5Object(counts = megacontacts)
DefaultAssay(seurat_object) <- "dist"


plot_doheatmap(dataset = seurat_object, markers = rownames(seurat_object), slot = "", sort_var = "banksy_cluster", anno_var = c("banksy_cluster", "orig.ident"), anno_colors = list(seurat_object@misc$banksy_cluster_colors, c("#34D916", "#00D4E6", "#1E90FF")), column_split = "banksy_cluster", hm_limit = seq(-2, 2, length = 100))

rownames(tmp1$E14_scAB250kb_scale) <- paste0("E14_5_", rownames(tmp1$E14_scAB250kb_scale))
rownames(tmp1$E14_scAB500kb_scale) <- paste0("E14_5_", rownames(tmp1$E14_scAB500kb_scale))
rownames(tmp1$E16_scAB500kb_scale) <- paste0("E16_5_", rownames(tmp1$E16_scAB500kb_scale))
rownames(tmp1$E16_scAB250kb_scale) <- paste0("E16_5_", rownames(tmp1$E16_scAB250kb_scale))
rownames(tmp1$E18_scAB250kb_scale) <- paste0("E18_5_", rownames(tmp1$E18_scAB250kb_scale))
rownames(tmp1$E18_scAB500kb_scale) <- paste0("E18_5_", rownames(tmp1$E18_scAB500kb_scale))
#!----------------------------------------------------------------------------------------------------------------------------------------------------------
# correlation clustering analysis
seurat_object <- qread("output/brain_development_Spatial_HiC_final_tutorial.qs")

seurat_object$Region_Cluster <- NA
seurat_object$Region_Cluster[which(seurat_object$banksy_cluster %in% paste0("E14_5_C", c(1:5)))] <- "GZ"
seurat_object$Region_Cluster[which(seurat_object$banksy_cluster %in% paste0("E14_5_C", c(6, 7)))] <- "CP"
seurat_object$Region_Cluster[which(seurat_object$banksy_cluster %in% paste0("E14_5_C", c(10)))] <- "MGE"
seurat_object$Region_Cluster[which(seurat_object$banksy_cluster %in% paste0("E14_5_C", c(8, 9, 11, 14)))] <- "LGE"
seurat_object$Region_Cluster[which(seurat_object$banksy_cluster %in% paste0("E14_5_C", c(12, 13)))] <- "CPU"
seurat_object$Region_Cluster[which(seurat_object$banksy_cluster %in% paste0("E14_5_C", c(15:17)))] <- "PIR"
seurat_object$Region_Cluster[which(seurat_object$banksy_cluster %in% paste0("E14_5_C", c(18)))] <- "LS"
seurat_object$Region_Cluster[which(seurat_object$banksy_cluster %in% paste0("E14_5_C", c(19)))] <- "Others"

seurat_object$Region_Cluster[which(seurat_object$banksy_cluster %in% paste0("E16_5_C", c(1:6)))] <- "GZ"
seurat_object$Region_Cluster[which(seurat_object$banksy_cluster %in% paste0("E16_5_C", c(7:10)))] <- "CP"
seurat_object$Region_Cluster[which(seurat_object$banksy_cluster %in% paste0("E16_5_C", c(11,12)))] <- "LGE"
seurat_object$Region_Cluster[which(seurat_object$banksy_cluster %in% paste0("E16_5_C", c(13, 14)))] <- "CPU"
seurat_object$Region_Cluster[which(seurat_object$banksy_cluster %in% paste0("E16_5_C", c(15:17)))] <- "PIR"
seurat_object$Region_Cluster[which(seurat_object$banksy_cluster %in% paste0("E16_5_C", c(18)))] <- "LS"
seurat_object$Region_Cluster[which(seurat_object$banksy_cluster %in% paste0("E16_5_C", c(19)))] <- "Others"

seurat_object$Region_Cluster[which(seurat_object$banksy_cluster %in% paste0("E18_5_C", c(1:4)))] <- "GZ"
seurat_object$Region_Cluster[which(seurat_object$banksy_cluster %in% paste0("E18_5_C", c(5:8)))] <- "CP"
seurat_object$Region_Cluster[which(seurat_object$banksy_cluster %in% paste0("E18_5_C", c(9)))] <- "LGE"
seurat_object$Region_Cluster[which(seurat_object$banksy_cluster %in% paste0("E18_5_C", c(10, 11)))] <- "CPU"
seurat_object$Region_Cluster[which(seurat_object$banksy_cluster %in% paste0("E18_5_C", c(12)))] <- "LS"
seurat_object$Region_Cluster[which(seurat_object$banksy_cluster %in% paste0("E18_5_C", c(13)))] <- "Others"

seurat_object$Region_Cluster <- factor(seurat_object$Region_Cluster, levels = c("GZ", "CP", "MGE", "LGE", "CPU", "PIR", "LS", "Others"))

library(SingleCellExperiment)
library(MetaNeighbor)
library(gplots)
library(RColorBrewer)
cols = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0 , 1,length = 101)

study_ID = "orig.ident"
cell_type = "banksy_cluster"
DefaultAssay(seurat_object) <- "scAB500kb_scale"
seurat_sce <- as.SingleCellExperiment(seurat_object, assay = "scAB500kb_scale")

study_ID = seurat_sce[[study_ID]]
cell_type = seurat_sce[[cell_type]]
global_hvgs = rownames(seurat_object)
aurocs = MetaNeighborUS(var_genes = global_hvgs, dat = seurat_sce, study_id = study_ID, cell_type = cell_type, fast_version = T)

tmp <- as.data.frame(strsplit(colnames(aurocs), split = '\\|'))

pdf(file = paste0("MetaNeighborUS_Spatial_HiC_banksy_cluster_heatmap.pdf"), width = 8, height = 8)
	rownames(aurocs) <- substr(rownames(aurocs), start = nchar(unique(as.character(study_ID))[1]) + 2, stop = nchar(rownames(aurocs)))
	colnames(aurocs) <- substr(colnames(aurocs), start = nchar(unique(as.character(study_ID))[1]) + 2, stop = nchar(colnames(aurocs)))
	MetaNeighbor::plotHeatmap(aurocs, cex = 0.6)
dev.off()

sample <- tmp[1,]
names(sample) <- NULL
major_types <- tmp[2,]
names(major_types) <- NULL

sample <- factor(sample, levels(seurat_object$orig.ident))
major_types <- factor(major_types, levels(seurat_object$banksy_cluster))

Region <- unique(seurat_object@meta.data[,c("Region_Cluster", "banksy_cluster")])
Region <- Region$Region_Cluster[match(major_types, Region$banksy_cluster)]

S_L_ratio <- c()
for(i in rownames(aurocs)){
	S_L_ratio <- c(S_L_ratio, mean(seurat_object$Short_Long_Ratio[which(seurat_object$banksy_cluster == i)]))
}
S_L_Ratio <- as.character(S_L_ratio)
S_L_Ratio[which(S_L_ratio <= 1.4)] <- "<=1.4"
S_L_Ratio[which(S_L_ratio > 1.4 & S_L_ratio <= 1.6)] <- "1.4~1.6"
S_L_Ratio[which(S_L_ratio > 1.6 & S_L_ratio <= 1.8)] <- "1.6~1.8"
S_L_Ratio[which(S_L_ratio > 1.8 & S_L_ratio <= 2)] <- "1.8~2"
S_L_Ratio[which(S_L_ratio > 2 & S_L_ratio <= 2.2)] <- "2~2.2"
S_L_Ratio[which(S_L_ratio > 2.2 & S_L_ratio <= 2.4)] <- "2.2~2.4"
S_L_Ratio[which(S_L_ratio >= 2.4 )] <- ">=2.4"
S_L_Ratio <- factor(S_L_Ratio, levels = c("<=1.4", "1.4~1.6", "1.6~1.8", "1.8~2", "2~2.2", "2.2~2.4", ">=2.4"))

S_L_Ratio_color <- viridis::turbo(length(levels(S_L_Ratio)))
S_L_Ratio_color <- paletteContinuous(set = "sambaNight", n = length(levels(S_L_Ratio)))
names(S_L_Ratio_color) <- levels(S_L_Ratio)

ann_row_1 <- rowAnnotation(a = major_types, border = TRUE, col = list(a = seurat_object@misc$banksy_cluster_colors[match(levels(major_types), names(seurat_object@misc$banksy_cluster_colors))]), annotation_label = "banksy_cluster", annotation_legend_param = list(title = "banksy_cluster"))
ann_row_2 <- rowAnnotation(a = sample, border = TRUE, col = list(a = cluster_colors[match(levels(sample), names(cluster_colors))]), annotation_label = "sample", annotation_legend_param = list(title = "sample"))
ann_row_3 <- rowAnnotation(a = Region, border = TRUE, col = list(a = cluster_colors[match(levels(Region), names(cluster_colors))]), annotation_label = "Region", annotation_legend_param = list(title = "Region"))
ann_row_4 <- rowAnnotation(a = S_L_Ratio, border = TRUE, col = list(a = S_L_Ratio_color), annotation_label = "S_L_Ratio", annotation_legend_param = list(title = "S_L_Ratio"))

ann_row <- c(ann_row_2, ann_row_1, ann_row_4, ann_row_3)
#ann_row@gap <- rep(unit(1, "mm"), length(ann_row))

ann_col_1 <- HeatmapAnnotation(a = major_types, border = TRUE, col = list(a = seurat_object@misc$banksy_cluster_colors[match(levels(major_types), names(seurat_object@misc$banksy_cluster_colors))]), annotation_label = "banksy_cluster", show_legend = F, annotation_name_side = "left")
ann_col_2 <- HeatmapAnnotation(a = sample, border = TRUE, col = list(a = cluster_colors[match(levels(sample), names(cluster_colors))]), annotation_label = "sample", show_legend = F, annotation_name_side = "left")
ann_col_3 <- HeatmapAnnotation(a = Region, border = TRUE, col = list(a = cluster_colors[match(levels(Region), names(cluster_colors))]), annotation_label = "Region", show_legend = F, annotation_name_side = list(title = "left"))
ann_col_4 <- HeatmapAnnotation(a = S_L_Ratio, border = TRUE, col = list(a = S_L_Ratio_color), annotation_label = "S_L_Ratio", show_legend = F, annotation_name_side = list(title = "left"))

ann_col <- c(ann_col_2, ann_col_1, ann_col_4, ann_col_3)
#ann_col@gap <- rep(unit(1, "mm"), length(ann_col))

ordering <- stats::as.dendrogram(orderCellTypes(aurocs), hang = 1)

ident_labels <- substr(rownames(aurocs), start = nchar(unique(as.character(study_ID))[1]) + 2, stop = nchar(rownames(aurocs)))

pdf(paste0("MetaNeighborUS_Spatial_HiC_banksy_cluster_heatmap1.pdf"), width = 14, height = 13)

Heatmap(aurocs, col = cols, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8), row_dend_width = unit(70, 'mm'), column_dend_height = unit(70, 'mm'), cluster_rows = ordering, cluster_columns = ordering, 
		left_annotation = ann_row, top_annotation = ann_col, row_labels = ident_labels, column_labels = ident_labels,
		# show_row_names = F, show_column_names = F,
		heatmap_legend_param = list(title_position = "topcenter", legend_width = unit(5, "cm"), legend_height = unit(3, "cm"), title = "AUROC"))

dev.off()

study_ID = "orig.ident"
cell_type = "Region"
DefaultAssay(seurat_object) <- "scAB500kb_scale"
seurat_sce <- as.SingleCellExperiment(seurat_object, assay = "scAB500kb_scale")

study_ID = seurat_sce[[study_ID]]
cell_type = seurat_sce[[cell_type]]
global_hvgs = rownames(seurat_object)
aurocs = MetaNeighborUS(var_genes = global_hvgs, dat = seurat_sce, study_id = study_ID, cell_type = cell_type, fast_version = T)

pdf(file = paste0("MetaNeighborUS_Spatial_HiC_Region_heatmap.pdf"), width = 8, height = 8)
	MetaNeighbor::plotHeatmap(aurocs, cex = 0.6)
dev.off()
#!----------------------------------------------------------------------------------------------------------------------------------------------------------
seurat_object <- qread("output/brain_development_Spatial_HiC_final_tutorial_cortex_new.qs")
seurat_object[["scAB250kb"]] <- CreateAssayObject(GetAssayData(seurat_object, assay = "scAB250kb_scale"))
DefaultAssay(seurat_object) <- "scAB250kb"
seurat_object <- FindVariableFeatures(seurat_object, assay = "scAB250kb", selection.method = "mvp")
gsea_file <- file("GO_BP_neuron_associated_terms", open="rt")
gsea_gene_set <- list()
count <- 1
while(TRUE){ 
        line <- readLines(gsea_file, n = 1)
        if(length(line) == 0) break
        gene_set <- strsplit(line, "\t")[[1]]
        gsea_gene_set[[count]] <- gene_set[3:length(gene_set)]
        names(gsea_gene_set)[count] <- gene_set[1]
        count <- count + 1
}
close(gsea_file)

DefaultAssay(seurat_object) <- "scAB250kb"
for(i in names(gsea_gene_set)){
	gsea_gene_set[[i]] <- intersect(rownames(seurat_object), unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name %in% gsea_gene_set[[i]])]))
	seurat_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(seurat_object, assay = "scAB250kb", slot = "data")[gsea_gene_set[[i]],]))
}

pdf("GO_BP_neuron_associated_terms_cortex_mean_scAB.pdf", width = 16)
for(i in names(gsea_gene_set)){
	p0 <- VlnPlot(seurat_object, group.by = "orig.ident", split.by = "banksy_cluster", pt.size = 0, cols = seurat_object@misc$banksy_cluster_colors[match(levels(seurat_object$banksy_cluster), names(seurat_object@misc$banksy_cluster_colors))], features = i) & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.1, position = position_dodge(0.9), show.legend = F) & NoLegend()
	print(p0)
}
dev.off()

pdf("GO_BP_neuron_associated_terms_cortex_mean_scAB_banksy_enrich.pdf", width = 10)
for(i in names(gsea_gene_set)){
	p0 <- VlnPlot(seurat_object, group.by = "banksy_enrich", split.by = "orig.ident", pt.size = 0, cols = c("#34D916", "#00D4E6", "#1E90FF"), features = i) & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.1, position = position_dodge(0.9), show.legend = F) & NoLegend()
	print(p0)
}
dev.off()

#AP analysis-----------------------------------------------------------------------------------------------------
Idents(seurat_object) <- "banksy_cluster"

E14_5_C1_C2 <- FindMarkers(seurat_object, ident.1 = "E14_5_C1", ident.2 = "E14_5_C2", min.pct = 0.1, test.use = 'LR', fc.name = "avg_diff")
de_markers <- E14_5_C1_C2[which(abs(E14_5_C1_C2$avg_diff) > 0.25 & E14_5_C1_C2$p_val_adj < 0.05),]
de_markers$gene <- rownames(de_markers)
de_markers$clusters <- "C1"
de_markers$clusters[which(de_markers$avg_diff < 0)] <- "C2"
de_markers$clusters <- factor(de_markers$clusters, levels = c("C1", "C2"))
write.csv(de_markers, "E14_5_C1_vs_E14_5_C2_scAB.csv", quote = F, row.names = F)

C1_features <- unique(seurat_object@misc$gene_infor_250k$gene_name[which(seurat_object@misc$gene_infor_250k$gene_region %in% de_markers$gene[which(de_markers$clusters == "C1")] & seurat_object@misc$gene_infor_250k$gene_type == "protein_coding")])
C2_features <- unique(seurat_object@misc$gene_infor_250k$gene_name[which(seurat_object@misc$gene_infor_250k$gene_region %in% de_markers$gene[which(de_markers$clusters == "C2")] & seurat_object@misc$gene_infor_250k$gene_type == "protein_coding")])

seurat_rna <- qread("output/brain_development_ST_final_tutorial_cortex.qs")
de_conserved_markers <- read.csv("Cortex_de_conserved_markers_cell_type.csv")
de_cell_type_markers <- read.csv("Cortex_cell_type_DE_markers.csv")
C1_features <- intersect(C1_features, c(de_conserved_markers$gene, de_cell_type_markers$gene))
C2_features <- intersect(C2_features, c(de_conserved_markers$gene, de_cell_type_markers$gene))

meta_data <- data.frame(genes = c(C1_features, C2_features), clusters = "C1")
meta_data$clusters[match(C2_features, meta_data$genes)] <- "C2"
ht <- SCP::GroupHeatmap(seurat_rna, features = meta_data$genes, feature_split = meta_data$clusters,  slot = "counts", exp_method = "zscore", group.by = "cell_type", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_rna$cell_type), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = rev(RColorBrewer::brewer.pal(11, "RdYlBu")), features_label = meta_data$gene, nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot

de_markers <- de_markers[which(de_markers$gene %in% seurat_object@misc$gene_infor_250k$gene_region[which(seurat_object@misc$gene_infor_250k$gene_name %in% c(C1_features, C2_features))]),]

tmp_object <- subset(seurat_object, subset = orig.ident == "E14_5")
tmp_object$banksy_cluster <- factor(tmp_object$banksy_cluster, levels = intersect(levels(seurat_object$banksy_cluster), unique(tmp_object$banksy_cluster)))
ht <- SCP::GroupHeatmap(tmp_object, features = de_markers$gene, feature_split = de_markers$clusters, slot = "data", exp_method = "zscore", group.by = "banksy_cluster", group_palcolor = list(seurat_object@misc$banksy_cluster_colors[match(levels(seurat_object$banksy_cluster), names(seurat_object@misc$banksy_cluster_colors))]), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot

meta_data <- seurat_object@misc$gene_infor_250k[which(seurat_object@misc$gene_infor_250k$gene_region %in% de_markers$gene & seurat_object@misc$gene_infor_250k$gene_name %in% c(C1_features, C2_features) & seurat_object@misc$gene_infor_250k$gene_type == "protein_coding"),]
meta_data <- unique(meta_data[,-c(1, which(colnames(meta_data) == "gene_id"))])
meta_data$region_gene <- paste0(meta_data$gene_region, "-", meta_data$gene_name)
meta_data$cluster <- de_markers$clusters[match(meta_data$gene_region, de_markers$gene)]
meta_data$idx <- paste0(meta_data$gene_name, "-", 1:nrow(meta_data))

scAB_matrix <- GetAssayData(tmp_object, assay = "scAB250kb", slot = "counts")
scAB_matrix <- scAB_matrix[meta_data$gene_region,]
rownames(scAB_matrix) <- meta_data$idx
tmp_object[["overlap"]] <- CreateAssayObject(scAB_matrix)

rna_matrix <- GetAssayData(seurat_rna, assay = "SCT", slot = "counts")
rna_matrix <- rna_matrix[meta_data$gene_name,]
rownames(rna_matrix) <- meta_data$idx
seurat_rna[["overlap"]] <- CreateAssayObject(rna_matrix)

ht1 <- SCP::GroupHeatmap(tmp_object, features = meta_data$idx, feature_split = meta_data$cluster, assay = "overlap", slot = "data", exp_method = "zscore", group.by = "banksy_cluster", group_palcolor = list(tmp_object@misc$banksy_cluster_colors[match(levels(tmp_object$banksy_cluster), names(tmp_object@misc$banksy_cluster_colors))]), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE), features_label = meta_data$idx)
ht1$plot

ht2 <- SCP::GroupHeatmap(seurat_rna, features = ht1$feature_metadata$features[order(ht1$feature_metadata$index)], feature_split = ht1$feature_metadata$feature_split[order(ht1$feature_metadata$index)], assay = "overlap", slot = "counts", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = FALSE, cluster_columns = FALSE, heatmap_palcolor = rev(RColorBrewer::brewer.pal(11, "RdYlBu")), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht2$plot
#!---------------------------------------------------------------------------------------------------------------
Idents(seurat_object) <- "banksy_enrich_new"
# identity conserved markers cross time
conserved_markers <- list()
for(i in levels(seurat_object$banksy_enrich_new)){
	conserved_markers[[i]] <- FindConservedMarkers(seurat_object, ident.1 = i, grouping.var = "orig.ident", min.pct = 0.1, test.use = 'LR', fc.name = "avg_diff", verbose = FALSE, assay = "scAB250kb")
}

all_conserved_markers <- c()
for(i in levels(seurat_object$banksy_enrich_new)){
	conserved_markers[[i]]$gene <- rownames(conserved_markers[[i]])
	conserved_markers[[i]]$cluster <- i
	all_conserved_markers <- rbind(all_conserved_markers, conserved_markers[[i]])
}
for(i in 1:nrow(all_conserved_markers)){
	all_conserved_markers$max_pval_adj[i] <- max(all_conserved_markers[i,c("E14_5_p_val_adj", "E16_5_p_val_adj", "E18_5_p_val_adj")])
}
all_conserved_markers <- all_conserved_markers[which(all_conserved_markers$max_pval_adj < 0.05),]
write.csv(all_conserved_markers, "Cortex_de_conserved_markers_HiC_banksy_enrich_new.csv", quote = F, row.names = F)

ht <- SCP::GroupHeatmap(seurat_object, features = unique(all_conserved_markers$gene), slot = "data", exp_method = "zscore", group.by = "banksy_enrich_new", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$banksy_enrich_new), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, feature_split_palcolor = list(cluster_colors[match(levels(seurat_object$banksy_enrich_new), names(cluster_colors))]), heatmap_palette = "viridis", nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot
pdf("Cortex_banksy_enrich_new_conserved_features.pdf", width = 12, height = 12)
ht$plot
dev.off()
#!---------------------------------------------------------------------------------------------------
ht <- SCP::GroupHeatmap(seurat_object, features = VariableFeatures(seurat_object), slot = "data", exp_method = "zscore", group.by = "orig.ident", split.by = "banksy_cluster", group_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cell_split_palcolor = list(seurat_object@misc$banksy_cluster_colors[match(levels(seurat_object$banksy_cluster), names(seurat_object@misc$banksy_cluster_colors))]), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot

ht <- SCP::GroupHeatmap(seurat_object, features = VariableFeatures(seurat_object), slot = "data", exp_method = "zscore", group.by = "orig.ident", split.by = "banksy_enrich", group_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cell_split_palcolor = list(cluster_colors[match(levels(seurat_object$banksy_enrich), names(cluster_colors))]), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot

ht <- SCP::GroupHeatmap(seurat_object, features = VariableFeatures(seurat_object), slot = "data", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot

seurat_object$time_Region <- paste0(seurat_object$orig.ident, "_", seurat_object$Region)
seurat_object$time_Region <- factor(seurat_object$time_Region, levels = c("E14_5_GZ", "E14_5_CP", "E16_5_GZ", "E16_5_CP", "E18_5_GZ", "E18_5_CP"))
Idents(seurat_object) <- "time_Region"

E14_5_GZ_CP <- FindMarkers(seurat_object, ident.1 = "E14_5_CP", ident.2 = "E14_5_GZ", min.pct = 0.1, test.use = 'LR', fc.name = "avg_diff")
E16_5_GZ_CP <- FindMarkers(seurat_object, ident.1 = "E16_5_CP", ident.2 = "E16_5_GZ", min.pct = 0.1, test.use = 'LR', fc.name = "avg_diff")
E18_5_GZ_CP <- FindMarkers(seurat_object, ident.1 = "E18_5_CP", ident.2 = "E18_5_GZ", min.pct = 0.1, test.use = 'LR', fc.name = "avg_diff")

E14_5_E16_5_CP <- FindMarkers(seurat_object, ident.1 = "E16_5_CP", ident.2 = "E14_5_CP", min.pct = 0.1, test.use = 'LR', fc.name = "avg_diff")
E14_5_E18_5_CP <- FindMarkers(seurat_object, ident.1 = "E18_5_CP", ident.2 = "E14_5_CP", min.pct = 0.1, test.use = 'LR', fc.name = "avg_diff")
E16_5_E18_5_CP <- FindMarkers(seurat_object, ident.1 = "E18_5_CP", ident.2 = "E16_5_CP", min.pct = 0.1, test.use = 'LR', fc.name = "avg_diff")

E14_5_E16_5_GZ <- FindMarkers(seurat_object, ident.1 = "E16_5_GZ", ident.2 = "E14_5_GZ", min.pct = 0.1, test.use = 'LR', fc.name = "avg_diff")
E14_5_E18_5_GZ <- FindMarkers(seurat_object, ident.1 = "E18_5_GZ", ident.2 = "E14_5_GZ", min.pct = 0.1, test.use = 'LR', fc.name = "avg_diff")
E16_5_E18_5_GZ <- FindMarkers(seurat_object, ident.1 = "E18_5_GZ", ident.2 = "E16_5_GZ", min.pct = 0.1, test.use = 'LR', fc.name = "avg_diff")

E14_5_GZ_CP <- E14_5_GZ_CP %>% filter(abs(avg_diff) >= 0.5, p_val_adj < 0.05)
E16_5_GZ_CP <- E16_5_GZ_CP %>% filter(abs(avg_diff) >= 0.5, p_val_adj < 0.05)
E18_5_GZ_CP <- E18_5_GZ_CP %>% filter(abs(avg_diff) >= 0.5, p_val_adj < 0.05)

E14_5_E16_5_CP <- E14_5_E16_5_CP %>% filter(abs(avg_diff) >= 0.5, p_val_adj < 0.05)
E14_5_E18_5_CP <- E14_5_E18_5_CP %>% filter(abs(avg_diff) >= 0.5, p_val_adj < 0.05)
E16_5_E18_5_CP <- E16_5_E18_5_CP %>% filter(abs(avg_diff) >= 0.5, p_val_adj < 0.05)

E14_5_E16_5_GZ <- E14_5_E16_5_GZ %>% filter(abs(avg_diff) >= 0.5, p_val_adj < 0.05)
E14_5_E18_5_GZ <- E14_5_E18_5_GZ %>% filter(abs(avg_diff) >= 0.5, p_val_adj < 0.05)
E16_5_E18_5_GZ <- E16_5_E18_5_GZ %>% filter(abs(avg_diff) >= 0.5, p_val_adj < 0.05)

features <- unique(c(rownames(E14_5_GZ_CP), rownames(E16_5_GZ_CP), rownames(E18_5_GZ_CP), rownames(E14_5_E16_5_CP), rownames(E14_5_E18_5_CP), rownames(E16_5_E18_5_CP), rownames(E14_5_E16_5_GZ), rownames(E14_5_E18_5_GZ), rownames(E16_5_E18_5_GZ)))

ht <- SCP::GroupHeatmap(seurat_object, features = features, n_split = 8, assay = "scAB250kb", slot = "data", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot
pdf("GZ_CP_HiC_Clusters.pdf", width = 6, height = 8)
ht$plot
dev.off()
write.csv(ht$feature_metadata, "GZ_CP_HiC_Clusters.csv", quote = F, row.names = F)

rna_features <- unique(seurat_object@misc$gene_infor_250k$gene_name[which(seurat_object@misc$gene_infor_250k$gene_region %in% features & seurat_object@misc$gene_infor_250k$gene_type == "protein_coding")])

de_conserved_markers <- read.csv("Cortex_de_conserved_markers_cell_type.csv")
de_cell_type_markers <- read.csv("Cortex_cell_type_DE_markers.csv")
rna_features <- intersect(rna_features, c(de_conserved_markers$gene, de_cell_type_markers$gene))
ht <- SCP::GroupHeatmap(seurat_rna, features = rna_features, slot = "counts", exp_method = "zscore", group.by = "cell_type", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_rna$cell_type), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = rev(RColorBrewer::brewer.pal(11, "RdYlBu")), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot

features <- intersect(features, seurat_object@misc$gene_infor_250k$gene_region[which(seurat_object@misc$gene_infor_250k$gene_name %in% rna_features)])

ht <- SCP::GroupHeatmap(seurat_object, features = features, n_split = 7, assay = "scAB250kb", slot = "data", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot

meta_data <- seurat_object@misc$gene_infor_250k[which(seurat_object@misc$gene_infor_250k$gene_region %in% features & seurat_object@misc$gene_infor_250k$gene_name %in% rna_features & seurat_object@misc$gene_infor_250k$gene_type == "protein_coding"),]
meta_data <- unique(meta_data[,-c(1, which(colnames(meta_data) == "gene_id"))])
meta_data$region_gene <- paste0(meta_data$gene_region, "-", meta_data$gene_name)
meta_data$idx <- paste0(meta_data$gene_name, "-", 1:nrow(meta_data))

key_markers <- c("Calm1","Camk2a", "Slc1a2", "Lrrtm4", "Rasgrf2", "Slc9a9", "Cux2", "Unc5d", "Astn2", "Rorb", "Myo16", "Brinp3", "Satb2", "Slc24a3", "Fat3", "Zfpm2", "Pde1a", "Nos1ap", "Cdh18", "Neto1", "Ncald", "Garnl3", "Plp1", "Mag", "Prr5l", "Dock10", "Mbp", "Mobp", "Fht1", "Tlk1", "Pde10a", "Gng7", "Adcy5", "Penk", "Drd2", "Lsamp", "Mertk", "Plpp3", "Prex2", "Grip1", "Cdh13", "Grin3a", "Reln", "Sst", "Erbb4", "Gad1", "Pvalb", "Syt2", "Adarb2", "Adra1a", "Vip", "Slco1a4", "Atp10a", "Flt1")

meta_data$label <- NA
meta_data$label[match(key_markers, meta_data$gene_name)] <- "T"

scAB_matrix <- GetAssayData(seurat_object, assay = "scAB250kb", slot = "counts")
scAB_matrix <- scAB_matrix[meta_data$gene_region,]
rownames(scAB_matrix) <- meta_data$idx
seurat_object[["overlap"]] <- CreateAssayObject(scAB_matrix)

rna_matrix <- GetAssayData(seurat_rna, assay = "SCT", slot = "counts")
rna_matrix <- rna_matrix[meta_data$gene_name,]
rownames(rna_matrix) <- meta_data$idx
seurat_rna[["overlap"]] <- CreateAssayObject(rna_matrix)

ht1 <- SCP::GroupHeatmap(seurat_object, features = meta_data$idx, n_split = 7, assay = "overlap", slot = "data", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE), features_label = meta_data$idx[which(meta_data$label == "T")])
ht1$plot

ht2 <- SCP::GroupHeatmap(seurat_rna, features = ht1$feature_metadata$features[order(ht1$feature_metadata$index)], feature_split = ht1$feature_metadata$feature_split[order(ht1$feature_metadata$index)], assay = "overlap", slot = "counts", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = T, cluster_columns = FALSE, heatmap_palcolor = rev(RColorBrewer::brewer.pal(11, "RdYlBu")), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE), features_label = meta_data$idx[which(meta_data$label == "T")])
ht2$plot

ht1$plot + ht2$plot
pdf("GZ_CP_HiC_Clusters_0.pdf", width = 10, height = 10)
ht1$plot
dev.off()
pdf("GZ_CP_HiC_Clusters_1.pdf", width = 10, height = 10)
ht2$plot
dev.off()

features <- intersect(features, seurat_object@misc$gene_infor_250k$gene_region[which(seurat_object@misc$gene_infor_250k$gene_name %in% rna_features)])

RNA_features <- read.csv("GZ_CP_DE_Clusters.csv")
HiC_features <- unique(seurat_object@misc$gene_infor_250k$gene_name[which(seurat_object@misc$gene_infor_250k$gene_region %in% features & seurat_object@misc$gene_infor_250k$gene_type == "protein_coding")])
overlap_features <- intersect(RNA_features$features, HiC_features)
idx_features <- intersect(seurat_object@misc$gene_infor_250k$gene_region[which(seurat_object@misc$gene_infor_250k$gene_name %in% overlap_features)], features)
meta_data <- seurat_object@misc$gene_infor_250k[which(seurat_object@misc$gene_infor_250k$gene_region %in% idx_features & seurat_object@misc$gene_infor_250k$gene_name %in% overlap_features & seurat_object@misc$gene_infor_250k$gene_type == "protein_coding"),]
meta_data <- unique(meta_data[,-c(1, which(colnames(meta_data) == "gene_id"))])
meta_data$region_gene <- paste0(meta_data$gene_region, "-", meta_data$gene_name)
meta_data$idx <- paste0(meta_data$gene_name, "-", 1:nrow(meta_data))
meta_data$label <- NA
meta_data$label[match(key_markers, meta_data$gene_name)] <- "T"

scAB_matrix <- GetAssayData(seurat_object, assay = "scAB250kb", slot = "counts")
scAB_matrix <- scAB_matrix[meta_data$gene_region,]
rownames(scAB_matrix) <- meta_data$idx
seurat_object[["overlap"]] <- CreateAssayObject(scAB_matrix)

rna_matrix <- GetAssayData(seurat_rna, assay = "SCT", slot = "counts")
rna_matrix <- rna_matrix[meta_data$gene_name,]
rownames(rna_matrix) <- meta_data$idx
seurat_rna[["overlap"]] <- CreateAssayObject(rna_matrix)

ht1 <- SCP::GroupHeatmap(seurat_object, features = meta_data$idx, n_split = 6, assay = "overlap", slot = "data", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE), features_label = meta_data$idx[which(meta_data$label == "T")])
ht1$plot

ht2 <- SCP::GroupHeatmap(seurat_rna, features = ht1$feature_metadata$features[order(ht1$feature_metadata$index)], feature_split = ht1$feature_metadata$feature_split[order(ht1$feature_metadata$index)], assay = "overlap", slot = "counts", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = T, cluster_columns = FALSE, heatmap_palcolor = rev(RColorBrewer::brewer.pal(11, "RdYlBu")), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE), features_label = meta_data$idx[which(meta_data$label == "T")])
ht2$plot

pdf("GZ_CP_overlap_Clusters_0.pdf", width = 10, height = 10)
ht1$plot
dev.off()

pdf("GZ_CP_overlap_Clusters_1.pdf", width = 10, height = 10)
ht2$plot
dev.off()

meta_data$cluster <- ht1$feature_metadata$feature_split[match(meta_data$idx, ht1$feature_metadata$features)]
write.csv(meta_data, "GZ_CP_overlap_Clusters_0.csv", quote = F, row.names = F)
dir.create("Cortex_Region_overlap_features")
for(i in unique(meta_data$cluster)){
        gene_list <- unique(as.character(meta_data[which(meta_data$cluster == i),'gene_name']))
        enrichRes <- enrichAnalysis(genelist = gene_list, geneType = "SYMBOL", species = "mmu", database = c("go", "kegg", "MSigDb"), GO.model = "BP", MSigDb.signature = "H", sampleName = paste0("enrich_", i), minGSSize = 5, outpath = "Cortex_Region_overlap_features")
}

DefaultAssay(seurat_rna) <- "spARC_SCT"
DefaultAssay(seurat_object) <- "spARC_scAB250kb_scale"
dir.create("Cortex_Region_overlap_features/features")
for(i in unique(meta_data$gene_name)){
	try({
		p0 <- Seurat::SpatialPlot(seurat_rna, features = i, shape = 22, crop = F) & scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")))
		print(p0)
		ggsave(paste0("Cortex_Region_overlap_features/features/", i, "_rna.png"), width = 24, height = 24)
		p1 <- Seurat::SpatialPlot(seurat_object, shape = 22, crop = F, features = seurat_object@misc$gene_infor_250k$gene_region[which(seurat_object@misc$gene_infor_250k$gene_name == i)]) & scale_fill_gradientn(values = seq(0, 1, 0.01), colours = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"))
		print(p1)
		ggsave(paste0("Cortex_Region_overlap_features/features/", i, "_hic.png"), width = 24, height = 24)
	})
}
#! moscot analysis------------------------------------------------------------------------------------------------------------------------------------------
transition_E14_E16 <- as.matrix(read.csv("Trajectory_analysis/moscot/tutorials/figures_HiC_all_cortex/E14_5_to_E16_5_transition_matrix_moscot_TemporalProblem_banksy_enrich.csv", row.names = 1))
transition_E14_E16[which(transition_E14_E16 <= 0.1)] <- 0
transition_E16_E18 <- as.matrix(read.csv("Trajectory_analysis/moscot/tutorials/figures_HiC_all_cortex/E16_5_to_E18_5_transition_matrix_moscot_TemporalProblem_banksy_enrich.csv", row.names = 1))
transition_E16_E18[which(transition_E16_E18 <= 0.1)] <- 0
transition_E18_P56 <- as.matrix(read.csv("Trajectory_analysis/moscot/tutorials/figures_HiC_all_cortex/E18_5_to_P56_transition_matrix_moscot_TemporalProblem_banksy_enrich.csv", row.names = 1))
transition_E18_P56[which(transition_E18_P56 <= 0.1)] <- 0

meta_data <- table(seurat_object$banksy_cluster, seurat_object$orig.ident)
meta_data <- apply(meta_data, 2, function(x){x/sum(x)})
meta_data[,1] <- round(meta_data[,1] * 10000)
colSums(meta_data)
meta_data0 <- data.frame("E14_5" = rep(rownames(meta_data), meta_data[, 1]))
#! deal transition_E14_E16------------------------------- 
tmp <- meta_data[,1][which(meta_data[,1] != 0)]
all(names(tmp) == rownames(transition_E14_E16))

transition_E14_E16 <- transition_E14_E16 * tmp
transition_E14_E16 <- floor(transition_E14_E16)
meta_data0$E16_5 <- NA
for(i in rownames(transition_E14_E16)){
	meta_data0$E16_5[which(meta_data0$E14_5 == i)][1:sum(transition_E14_E16[i,])] <- rep(colnames(transition_E14_E16), transition_E14_E16[i,])
}
#!-------------------------------------------------------
meta_data[,2] <- round(meta_data[,2] * 10000)
tmp <- meta_data[,2][which(meta_data[,2] != 0)]
all(names(tmp) == rownames(transition_E16_E18))

for(i in names(tmp)){
	tmp_data <- tmp[i] - length(which(meta_data0$E16_5 == i))
	if(tmp_data > 0){
		tmp_data <- data.frame("E14_5" = NA, "E16_5" = rep(i, tmp_data))
		meta_data0 <- rbind(meta_data0, tmp_data)
	}
}

transition_E16_E18 <- transition_E16_E18 * tmp
transition_E16_E18 <- floor(transition_E16_E18)
meta_data0$E18_5 <- NA
for(i in rownames(transition_E16_E18)){
	meta_data0$E18_5[which(meta_data0$E16_5 == i)][1:sum(transition_E16_E18[i,])] <- rep(colnames(transition_E16_E18), transition_E16_E18[i,])
}
#!-----------------------------------------------------
meta_data[,3] <- round(meta_data[,3] * 10000)
tmp <- meta_data[,3][which(meta_data[,3] != 0)]

for(i in names(tmp)){
	tmp_data <- tmp[i] - length(which(meta_data0$E18_5 == i))
	if(tmp_data > 0){
		tmp_data <- data.frame("E14_5" = NA, "E16_5" = NA, "E18_5" = rep(i, tmp_data))
		meta_data0 <- rbind(meta_data0, tmp_data)
	}
}

transition_E18_P56 <- transition_E18_P56 * tmp
transition_E18_P56 <- floor(transition_E18_P56)
meta_data0$P56 <- NA
for(i in rownames(transition_E18_P56)){
	meta_data0$P56[which(meta_data0$E18_5 == i)][1:sum(transition_E18_P56[i,])] <- rep(colnames(transition_E18_P56), transition_E18_P56[i,])
}
#!-----------------------------------------------------
meta_data[,4] <- round(meta_data[,4] * 10000)
tmp <- meta_data[,4][which(meta_data[,4] != 0)]

for(i in names(tmp)){
	tmp_data <- tmp[i] - length(which(meta_data0$P56 == i))
	if(tmp_data > 0){
		tmp_data <- data.frame("E14_5" = NA, "E16_5" = NA, "E18_5" = NA, "P56" = rep(i, tmp_data))
		meta_data0 <- rbind(meta_data0, tmp_data)
	}
}

meta_data <- meta_data0

library(ggsankey)

meta_data <- meta_data %>% make_long(E14_5, E16_5, E18_5, P56)
meta_data$node <- factor(meta_data$node, levels = c(levels(seurat_object$banksy_cluster)[1:24], "C5", "C4", "C3", "C2"))
meta_data$next_node <- factor(meta_data$next_node, levels = c(levels(seurat_object$banksy_cluster)[1:24], "C5", "C4", "C3", "C2"))

ggplot(meta_data, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) + geom_sankey(color = "black", lwd = 0.1, flow.alpha = 0.5, show.legend = FALSE, na.rm = FALSE) + scale_x_discrete(expand = c(0, 0.8)) + 
		scale_fill_manual(values = seurat_object@misc$banksy_cluster_colors, drop = FALSE) + 
		geom_sankey_text(size = 3, color = "black") + theme_void() + theme(axis.text.x = element_text())
ggsave("Trajectory_analysis/moscot/tutorials/figures_HiC_all_cortex/cortex_transition_matrix_moscot_TemporalProblem_sankeyplot.pdf")

meta_data <- read.csv("Trajectory_analysis/moscot/tutorials/figures_HiC_all_cortex/brain_development_HiC_final_all_cortex_moscot_cellrank.csv", row.names = 1)
meta_data <- meta_data[colnames(seurat_object),]
seurat_object@meta.data <- cbind(seurat_object@meta.data, meta_data[,14:27])
Seurat::SpatialPlot(seurat_object, shape = 22, crop = F, features = "IPC_MigN_mapping_14_5_to_16_5", images = c("E14_5", "E16_5")) & scale_fill_gradientn(colours = viridis::viridis(15, direction = -1))
#! E14_5 IPC_MigN push--------------------------------------------------------------------------
Seurat::SpatialPlot(seurat_object, shape = 22, crop = F, features = "IPC_MigN_mapping_14_5_to_16_5", images = c("E16_5")) & scale_fill_gradientn(colours = viridis::viridis(15, direction = -1))

seurat_object$E14_IPC_MigN_push <- NA
seurat_object$E14_IPC_MigN_push[which(seurat_object$banksy_enrich == "IPC_MigN" & seurat_object$orig.ident == "E14_5")] <- "source"
seurat_object$E14_IPC_MigN_push[which(seurat_object$IPC_MigN_mapping_14_5_to_16_5 > 0.0001 & seurat_object$orig.ident == "E16_5"  & seurat_object$banksy_enrich == "IPC_MigN")] <- "IPC_MigN"
seurat_object$E14_IPC_MigN_push[which(seurat_object$IPC_MigN_mapping_14_5_to_16_5 > 0.0001 & seurat_object$orig.ident == "E16_5"  & seurat_object$banksy_enrich == "CThPN_SCPN")] <- "CThPN_SCPN"
seurat_object$E14_IPC_MigN_push <- factor(seurat_object$E14_IPC_MigN_push, levels = c("IPC_MigN", "CThPN_SCPN", "source"))
Seurat::SpatialPlot(seurat_object, shape = 22, crop = F, group.by = "E14_IPC_MigN_push", images = c("E16_5"), cols = cluster_colors)

E14_IPC_MigN_push <- list()
for(i in c("IPC_MigN", "CThPN_SCPN")){
	seurat_object$tmp <- as.character(seurat_object$E14_IPC_MigN_push)
	seurat_object$tmp[which(seurat_object$tmp == "source")] <- i
	Idents(seurat_object) <- "tmp"
	E14_IPC_MigN_push[[i]] <- FindConservedMarkers(seurat_object, ident.1 = i, grouping.var = "orig.ident", min.pct = 0.1, only.pos = T, verbose = FALSE, assay = "scAB250kb_scale", test.use = 'LR', fc.name = "avg_diff")
}
E14_IPC_MigN_push_conserved_markers <- c()
for(i in c("IPC_MigN", "CThPN_SCPN")){
	E14_IPC_MigN_push[[i]]$gene <- rownames(E14_IPC_MigN_push[[i]])
	E14_IPC_MigN_push[[i]]$cluster <- i
	E14_IPC_MigN_push_conserved_markers <- rbind(E14_IPC_MigN_push_conserved_markers, E14_IPC_MigN_push[[i]])
}
E14_IPC_MigN_push_conserved_markers$max_pval_adj <- NA
for(i in 1:nrow(E14_IPC_MigN_push_conserved_markers)){
	E14_IPC_MigN_push_conserved_markers$max_pval_adj[i] <- max(E14_IPC_MigN_push_conserved_markers[i,c("E14_5_p_val_adj", "E16_5_p_val_adj")])
}

E14_IPC_MigN_push_conserved_markers <- E14_IPC_MigN_push_conserved_markers[which(E14_IPC_MigN_push_conserved_markers$max_pval_adj < 0.05),]
write.csv(E14_IPC_MigN_push_conserved_markers, "E14_IPC_MigN_push_conserved_markers_HiC.csv", quote = F, row.names = F)

ht <- SCP::GroupHeatmap(seurat_object, features = E14_IPC_MigN_push_conserved_markers$gene[which(E14_IPC_MigN_push_conserved_markers$cluster == "CThPN_SCPN")], slot = "data", exp_method = "zscore", group.by = "banksy_enrich", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$banksy_enrich), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palette = "viridis", nlabel = 0, show_row_names = T, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot

HiC_features <- unique(seurat_object@misc$gene_infor_250k$gene_name[which(seurat_object@misc$gene_infor_250k$gene_region %in% E14_IPC_MigN_push_conserved_markers$gene[which(E14_IPC_MigN_push_conserved_markers$cluster == "CThPN_SCPN")] & seurat_object@misc$gene_infor_250k$gene_type == "protein_coding")])
#!-------------------------------------------------
TF_file_list <- list.files("Trajectory_analysis/moscot/tutorials/figures_HiC_all_cortex/")
TF_file_list <- TF_file_list[grep("_push.csv", TF_file_list)]
TF_file_list <- gsub(TF_file_list, pattern = ".csv", replacement = "")
TF_file_list <- TF_file_list[-grep("E18_5_P56", TF_file_list)]
core_list <- gsub(TF_file_list, pattern = "_push", replacement = "")
core_list <- gsub(core_list, pattern = "_drivers_tf", replacement = "")
hic_core_regions <- list()
hic_core_genes <- list()
for(i in 1:length(TF_file_list)){
	tmp_data <- read.csv(paste0("Trajectory_analysis/moscot/tutorials/figures_HiC_all_cortex/", TF_file_list[i], ".csv"))
	colnames(tmp_data) <- c("genes", "corr", "pval", "qval", "ci_low", "ci_high")
	tmp_data$FDR <- log10(tmp_data$qval)
	hic_core_regions[[core_list[i]]] <- tmp_data$genes[which(tmp_data$corr >= 0.2)]
	hic_core_genes[[core_list[i]]] <- unique(seurat_object@misc$gene_infor_250k$gene_name[which(seurat_object@misc$gene_infor_250k$gene_region %in% hic_core_regions[[core_list[i]]] & seurat_object@misc$gene_infor_250k$gene_type == "protein_coding")])
	if(length(grep("[0-9]Rik$", hic_core_genes[[core_list[i]]])) > 0){
		hic_core_genes[[core_list[i]]] <- hic_core_genes[[core_list[i]]][-grep("[0-9]Rik$", hic_core_genes[[core_list[i]]])]
	}
	if(length(grep("Gm[0-9]", hic_core_genes[[core_list[i]]])) > 0){
		hic_core_genes[[core_list[i]]] <- hic_core_genes[[core_list[i]]][-grep("Gm[0-9]", hic_core_genes[[core_list[i]]])]
	}
	#p0 <- ggplot(tmp_data[which(tmp_data$corr > 0 & tmp_data$qval < 0.05),], aes(FDR, corr)) + 
	#	geom_point(aes(size = corr, color = qval)) + scale_color_continuous(low = "#e06663", high = "#327eba", name = 'p.adjust', guide = guide_colorbar(reverse = TRUE)) +
	#	scale_size_continuous(range = c(1, 3)) + theme_classic(base_size = 10) + xlab("log10(qval)") + ggtitle(i) +
	#	theme(panel.grid = element_blank(), axis.text = element_text(size = 10, color = "black", face = "bold"), axis.title = element_text(size = 12, color = "black", face = "bold")) + 
	#	geom_text_repel(data = filter(tmp_data, genes %in% core_gene), max.overlaps = getOption("ggrepel.max.overlaps", default = 15), aes(label = genes), size = 4, color = 'black')
	#print(p0)
	#ggsave(paste0("Trajectory_analysis/moscot/tutorials/figures_all_cortex/", i, ".pdf"))
}
TF_file_list <- list.files("Trajectory_analysis/moscot/tutorials/figures_HiC_all_cortex/")
TF_file_list <- TF_file_list[grep("_push_spearman.csv", TF_file_list)]
TF_file_list <- gsub(TF_file_list, pattern = ".csv", replacement = "")
TF_file_list <- TF_file_list[-grep("E18_5_P56", TF_file_list)]
core_list <- gsub(TF_file_list, pattern = "_push_spearman", replacement = "")
core_list <- gsub(core_list, pattern = "_drivers_tf", replacement = "")
for(i in 1:length(TF_file_list)){
	tmp_data <- read.csv(paste0("Trajectory_analysis/moscot/tutorials/figures_HiC_all_cortex/", TF_file_list[i], ".csv"))
	colnames(tmp_data) <- c("genes", "corr", "pval", "qval", "ci_low", "ci_high")
	tmp_data$FDR <- log10(tmp_data$qval)
	hic_core_regions[[core_list[i]]] <- unique(c(hic_core_regions[[core_list[i]]], tmp_data$genes[which(tmp_data$corr >= 0.2)]))
	hic_core_genes[[core_list[i]]] <- unique(c(hic_core_genes[[core_list[i]]], seurat_object@misc$gene_infor_250k$gene_name[which(seurat_object@misc$gene_infor_250k$gene_region %in% hic_core_regions[[core_list[i]]] & seurat_object@misc$gene_infor_250k$gene_type == "protein_coding")]))
	if(length(grep("[0-9]Rik$", hic_core_genes[[core_list[i]]])) > 0){
		hic_core_genes[[core_list[i]]] <- hic_core_genes[[core_list[i]]][-grep("[0-9]Rik$", hic_core_genes[[core_list[i]]])]
	}
	if(length(grep("Gm[0-9]", hic_core_genes[[core_list[i]]])) > 0){
		hic_core_genes[[core_list[i]]] <- hic_core_genes[[core_list[i]]][-grep("Gm[0-9]", hic_core_genes[[core_list[i]]])]
	}
	#p0 <- ggplot(tmp_data[which(tmp_data$corr > 0 & tmp_data$qval < 0.05),], aes(FDR, corr)) + 
	#	geom_point(aes(size = corr, color = qval)) + scale_color_continuous(low = "#e06663", high = "#327eba", name = 'p.adjust', guide = guide_colorbar(reverse = TRUE)) +
	#	scale_size_continuous(range = c(1, 3)) + theme_classic(base_size = 10) + xlab("log10(qval)") + ggtitle(i) +
	#	theme(panel.grid = element_blank(), axis.text = element_text(size = 10, color = "black", face = "bold"), axis.title = element_text(size = 12, color = "black", face = "bold")) + 
	#	geom_text_repel(data = filter(tmp_data, genes %in% core_gene), max.overlaps = getOption("ggrepel.max.overlaps", default = 15), aes(label = genes), size = 4, color = 'black')
	#print(p0)
	#ggsave(paste0("Trajectory_analysis/moscot/tutorials/figures_all_cortex/", i, ".pdf"))
}

#!---------------------------------------------------------------------------
library(SeuratWrappers)
library(Banksy)
library(harmony)

VariableFeatures(tmp_object) <- rownames(tmp_object)
tmp_object <- ScaleData(tmp_object, features = rownames(tmp_object))
tmp_object <- RunPCA(tmp_object, assay = "scAB500k", reduction.name = "pca", features = rownames(tmp_object), npcs = 30)
tmp_object <- RunHarmony(tmp_object, group.by.vars = "sample", reduction.name = 'pca', reduction.save = 'harmony')
tmp_object <- RunUMAP(tmp_object, reduction = "harmony", dims = 1:30, min.dist = 0.2, reduction.name = "harmony_umap", return.model = T, n.neighbors = 24)

tmp_object <- FindNeighbors(tmp_object, reduction = "harmony", dims = 1:30, annoy.metric = "cosine", k.param = 24)
tmp_object <- FindClusters(tmp_object, resolution = seq(0.1, 1.5, 0.1))
wrap_plots(map(seq(0.1, 1.5, 0.1), function(x) DimPlot(tmp_object, reduction = "harmony_umap", group.by = paste0("scAB500k_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)


seurat_object <- qread("output/brain_Adult16_HiC_final_tutorial.qs")

anchors <- FindTransferAnchors(reference = seurat_object, query = tmp_object, dims = 1:30, reference.reduction = "pca_banksy", npcs = NULL, query.assay = "scAB500k", reduction = "cca")

tmp_object <- MapQuery(anchorset = anchors, reference = seurat_object, query = tmp_object, refdata = list(celltype = "banksy_cluster"), reference.reduction = "pca_banksy", new.reduction.name = "ref.pca", reduction.model = "banksy_umap", projectumap.args = list(reduction.name = "ref.umap"), integrateembeddings.args = list(k.weight = 24), transferdata.args = list(k.weight = 24))

p0 <- DimPlot(tmp_object, reduction = "ref.umap", group.by = "predicted.celltype")
p1 <- DimPlot(seurat_object, label = T, cols = seurat_object@misc$Adult16_colors)
p0 + p1
#!----------------------------------------------------------------------------
library(ggsankey)

meta_data <- hic_object@meta.data[which(hic_object$orig.ident == "E18_5"),] %>% make_long(Region, banksy_cluster)
meta_data$node <- factor(meta_data$node, levels = rev(c(levels(hic_object$Region), levels(hic_object$banksy_cluster))))
meta_data$next_node <- factor(meta_data$next_node, levels = rev(c(levels(hic_object$Region), levels(hic_object$banksy_cluster))))

ggplot(meta_data, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) + geom_sankey(color = "black", lwd = 0.2, flow.alpha = 0.5, show.legend = FALSE, na.rm = FALSE) + scale_x_discrete(expand = c(0, 0.8)) + 
		scale_fill_manual(values = c(hic_object@misc$banksy_cluster_colors, cluster_colors), drop = FALSE) + geom_sankey_text(size = 3.2, color = "black") + theme_void() + theme(axis.text.x = element_text())
ggsave("Development_QC_cell_identity/E18_development_banksy_cluster_Region_sankey_plot.pdf")




















