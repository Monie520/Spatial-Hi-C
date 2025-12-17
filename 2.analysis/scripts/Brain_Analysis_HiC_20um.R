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

cluster_colors <- c("E14_5" = "#34D916", "E16_5" = "#00D4E6", "E18_5" = "#1E90FF", "Adult15" = "#B312A6", "GZ" = "#66C2A5", "CP" = "#FC8D62",
"Cortex" = "#59B375", "MGE" = "#9E7BFF", "LGE" = "#89288F", "ChP" = "#FF26A8", "LCS" = "#C1FF73", "CLA" = "#33FF00", "LS" = "#CC79A7", "Epd" = "#94FFB5", "PIR" = "#B3823E", "AEP" = "#653EB3", "others" = "#FEE52C", "Others" = "#dde2e6",
"AP" = "#03FFF4", "Cortex_DP_AP" = "#03FFF4", "Cortex_MP_AP" = "#036DF4", "IPC" = "#0BD3B1", "Cortex_IPC" = "#0BD3B1", "MigN" = "#62CFE8", "Cortex_DP_MigN" = "#62CFE8", "Cortex_MP_MigN" = "#2F7DD1",
"Cortex_Layer_6b" = "#99FFBC", "Layer_6b" = "#99FFBC", "CThPN" = "#7EC136", "Cortex_CThPN" = "#7EC136", "Cortex_DP_CThPN" = "#00CC14", "Cortex_MP_CThPN" = "#2D739E", "SCPN" = "#34A047", "Cortex_SCPN" = "#34A047", "Cortex_DP_SCPN" = "#00991F", "Cortex_MP_SCPN" = "#5C5CA6", "UL_CPN" = "#01545a", "Cortex_UL_CPN" = "#01545a",
"CTGL_ITL6GL" = "#66C5CC", "PTGL" = "#C9DB74", "ITL5GL" = "#87C55F", "ITL4GL" = "#54990F", "ITL23GL" = "#017351", "ITL1" = "#3DCCB1", "Astrocytes" = "#16F2F2", "CR" = "#F2F318", "Cortex_CR" = "#F2F318", "MGE_derived_InN" = "#D4E727", "Cortex_MGE_derived_InN" = "#D4E727",
"Layer_I" = "#3DCCB1", "Layer_II_III" = "#017351", "Layer_IV" = "#54990F", "Layer_V" = "#87C55F", "Layer_VI" = "#66C5CC", "CC" = "#64C2FC", "CPU" = "#B33E52",
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
gene_infor_500k <- read.csv("GRCm38_102_gtf_infor_hic_500kb_reference")
#!---------------------------------------------------------------------------
seurat_list <- qread("output/seurat_hic_devp_list_20um.qs")
seurat_RNA <- qread("output/seurat_list_20um.qs")

hic_object <- seurat_list$E14_5
locs <- Seurat::GetTissueCoordinates(hic_object)[,seq_len(2)]
colnames(locs) <- c("array_row", "array_col")
hic_object@meta.data <- cbind(hic_object@meta.data, locs)

DefaultAssay(hic_object) <- "scAB500k"
VariableFeatures(hic_object) <- rownames(hic_object)

hic_object <- RunBanksy(hic_object, lambda = 0.2, assay = 'scAB500k', slot = 'data', dimx = "array_row", dimy = "array_col", features = 'all', k_geom = 24, assay_name = "BANKSY_scAB500k")
npcs = 30
hic_object <- RunPCA(hic_object, assay = "BANKSY_scAB500k", reduction.name = "pca_banksy", npcs = npcs, features = rownames(hic_object))
hic_object <- RunUMAP(hic_object, reduction = "pca_banksy", dims = 1:npcs, min.dist = 0.3, reduction.name = "banksy_umap", return.model = T, n.neighbors = 24)
hic_object <- FindNeighbors(hic_object, reduction = "pca_banksy", dims = 1:npcs, annoy.metric = "cosine", graph.name = c('banksy_nn', 'banksy_snn'), k.param = 24)
hic_object <- FindClusters(hic_object, resolution = seq(0.1, 2, 0.1), graph.name = "banksy_snn")
wrap_plots(map(seq(0.6, 2, 0.1), function(x) Seurat::DimPlot(hic_object, group.by = paste0("banksy_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)
wrap_plots(map(seq(0.6, 2, 0.1), function(x) Seurat::SpatialPlot(hic_object, shape = 22, crop = F, group.by = paste0("banksy_snn_res.", x), cols = selected_colors, label = F) + NoLegend()), ncol = 5)
#!----------------------------------------------------------------------------
hic_object <- FindClusters(hic_object, cluster.name = "banksy_cluster", resolution = 0.9, graph.name = "banksy_snn")
Idents(hic_object) <- "banksy_cluster"
selected_cells <- CellsByIdentities(hic_object)
Seurat::SpatialPlot(hic_object, cells.highlight = selected_cells[setdiff(names(selected_cells), "NA")], cols.highlight = c("#FFFF00", alpha("lightgrey", 0)), crop = FALSE, shape = 22, facet.highlight = T, combine = T, ncol = 6, alpha = NULL) & NoLegend()

hic_object$seurat_clusters <- hic_object$banksy_cluster
hic_object$banksy_cluster <- as.character(hic_object$banksy_cluster)
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(10))] <- "C1"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(5))] <- "C2"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(2,9))] <- "C3"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(4, 17))] <- "C4"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(0))] <- "C5"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(15))] <- "C6"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(6))] <- "C7"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(1))] <- "C8"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(13))] <- "C9"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(11, 14))] <- "C10"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(7))] <- "C11"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(8))] <- "C12"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(3))] <- "C13"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(12))] <- "C14"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(16))] <- "C15" # Stromal

hic_object$banksy_cluster <- factor(hic_object$banksy_cluster, levels = paste0("C", 1:15))

Idents(hic_object) <- "banksy_cluster"

E14_colors <- c("C1" = "#4DFF99", "C2" = "#1FCCCC", "C3" = "#7ACC00", "C4" = "#459967", "C5" = "#9E7BFF", "C6" = "#BE50FF", "C7" = "#D85FF7", "C8" = "#D38B5C", "C9" = "#33FF00", "C10" = "#96C885", "C11" = "#CC79A7", "C12" = "#07F4B3", "C13" = "#FF9B4D", "C14" = "#B14380", "C15" = "#886C00")
DimPlot(hic_object, label = T, cols = E14_colors, reduction = "banksy_umap") + NoLegend() + SetAxes()
Seurat::SpatialPlot(hic_object, cols = E14_colors, crop = FALSE, shape = 22)
hic_object@misc$E14_colors <- E14_colors

seurat_rna <- seurat_RNA$E14_5
DefaultAssay(seurat_rna) <- "SCT"
all_markers <- FindAllMarkers(seurat_rna, only.pos = T, min.pct = 0.1)

de_markers <- all_markers[which(all_markers$avg_log2FC > 1 & all_markers$p_val_adj < 0.01 & all_markers$pct.1 > 0.3),]

marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:100]
	count = count + 1
}
names(marker_genes) <- unique(de_markers$cluster)

DefaultAssay(hic_object) <- "scAB500k"
for(i in names(marker_genes)){
	marker_genes[[i]] <- intersect(rownames(hic_object), unique(gene_infor_500k$gene_region[which(gene_infor_500k$gene_name %in% marker_genes[[i]])]))
	hic_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(hic_object, assay = "scAB500k", slot = "data")[marker_genes[[i]],]))
}

tmp_meta <- as.matrix(hic_object@meta.data[,levels(seurat_rna)])
rownames(tmp_meta) <- paste0(hic_object$banksy_cluster)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- paste0(levels(seurat_rna$Cluster), "_", colnames(tmp_meta))
tmp_meta0 <- tmp_meta0[levels(hic_object$banksy_cluster),]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("brain_20um_analysis/cell_identity_2_density_", names(seurat_list)[1], "_hic_500k_top100_zscore.pdf"), width = 10, height = 10)

seurat_list$E14_5 <- hic_object
#!---------------------------------------------------------------------------------------------------------------------------------------------------
hic_object <- seurat_list$E16_5
locs <- Seurat::GetTissueCoordinates(hic_object)[,seq_len(2)]
colnames(locs) <- c("array_row", "array_col")
hic_object@meta.data <- cbind(hic_object@meta.data, locs)

DefaultAssay(hic_object) <- "scAB500k"
VariableFeatures(hic_object) <- rownames(hic_object)

hic_object <- RunBanksy(hic_object, lambda = 0.4, assay = 'scAB500k', slot = 'data', dimx = "array_row", dimy = "array_col", features = 'all', k_geom = 24, assay_name = "BANKSY_scAB500k")
npcs = 30
hic_object <- RunPCA(hic_object, assay = "BANKSY_scAB500k", reduction.name = "pca_banksy", npcs = npcs, features = rownames(hic_object))
hic_object <- RunUMAP(hic_object, reduction = "pca_banksy", dims = 1:npcs, min.dist = 0.3, reduction.name = "banksy_umap", return.model = T, n.neighbors = 24)
hic_object <- FindNeighbors(hic_object, reduction = "pca_banksy", dims = 1:npcs, annoy.metric = "cosine", graph.name = c('banksy_nn', 'banksy_snn'), k.param = 24)
hic_object <- FindClusters(hic_object, resolution = seq(0.1, 2, 0.1), graph.name = "banksy_snn")
wrap_plots(map(seq(0.6, 2, 0.1), function(x) DimPlot(hic_object, reduction = "banksy_umap", group.by = paste0("banksy_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)
wrap_plots(map(seq(0.6, 2, 0.1), function(x) Seurat::SpatialPlot(hic_object, shape = 22, crop = F, group.by = paste0("banksy_snn_res.", x), cols = selected_colors, label = F) + NoLegend()), ncol = 5)
#!----------------------------------------------------------------------------
hic_object <- FindClusters(hic_object, cluster.name = "banksy_cluster", resolution = 1.5, graph.name = "banksy_snn")
hic_object <- FindSubCluster(hic_object, cluster = 3, graph.name = "banksy_snn", subcluster.name = "banksy_sub_cluster")
Idents(hic_object) <- "banksy_cluster"
selected_cells <- CellsByIdentities(hic_object)
Seurat::SpatialPlot(hic_object, cells.highlight = selected_cells[setdiff(names(selected_cells), "NA")], cols.highlight = c("#FFFF00", alpha("lightgrey", 0)), crop = FALSE, shape = 22, facet.highlight = T, combine = T, ncol = 6, alpha = NULL) & NoLegend()

hic_object$seurat_clusters <- hic_object$banksy_cluster
hic_object$banksy_cluster <- as.character(hic_object$banksy_cluster)
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(8))] <- "C1"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(5, 13))] <- "C2"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(2, 7, 21))] <- "C3"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(16))] <- "C4"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(22))] <- "C5"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(3))] <- "C6"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(10))] <- "C7"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(17))] <- "C8"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(12))] <- "C9"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(4))] <- "C10"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(20))] <- "C11"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(15))] <- "C12"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(6, 23))] <- "C13"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(1, 18))] <- "C14"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(0))] <- "C15"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(11))] <- "C16"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(14))] <- "C17"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(19))] <- "C18"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(9))] <- "C19"

hic_object$banksy_cluster <- factor(hic_object$banksy_cluster, levels = paste0("C", 1:19))

Idents(hic_object) <- "banksy_cluster"

E16_colors <- c("C1" = "#03FFF4", "C2" = "#0BD3B1", "C3" = "#009E73", "C4" = "#007756", "C5" = "#F2F318", "C6" = "#9E7BFF", "C7" = "#8D73DF", "C8" = "#BE50FF", "C9" = "#D85FF7", "C10" = "#D38B5C", "C11" = "#C1FF73", "C12" = "#33FF00", "C13" = "#B5EFB5", "C14" = "#96C885", "C15" = "#CC79A7", "C16" = "#2271A9", "C17" = "#FFD068", "C18" = "#B14380", "C19" = "#886C00")
DimPlot(hic_object, label = T, cols = E16_colors, reduction = "banksy_umap") + NoLegend() + SetAxes()
Seurat::SpatialPlot(hic_object, cols = E16_colors, crop = FALSE, shape = 22)
hic_object@misc$E16_colors <- E16_colors

seurat_rna <- seurat_RNA$E16_5
DefaultAssay(seurat_rna) <- "SCT"
all_markers <- FindAllMarkers(seurat_rna, only.pos = T, min.pct = 0.1)

de_markers <- all_markers[which(all_markers$avg_log2FC > 1 & all_markers$p_val_adj < 0.01 & all_markers$pct.1 > 0.3),]

marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:50]
	count = count + 1
}
names(marker_genes) <- unique(de_markers$cluster)
DefaultAssay(hic_object) <- "scAB500k"
for(i in names(marker_genes)){
	marker_genes[[i]] <- intersect(rownames(hic_object), unique(gene_infor_500k$gene_region[which(gene_infor_500k$gene_name %in% marker_genes[[i]])]))
	hic_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(hic_object, assay = "scAB500k", slot = "data")[marker_genes[[i]],]))
}

tmp_meta <- as.matrix(hic_object@meta.data[,levels(seurat_rna)])
rownames(tmp_meta) <- paste0(hic_object$banksy_cluster)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- paste0(levels(seurat_rna$Cluster), "_", colnames(tmp_meta))
tmp_meta0 <- tmp_meta0[levels(hic_object$banksy_cluster),]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("brain_20um_analysis/cell_identity_2_density_", names(seurat_list)[2], "_hic_500k_top50_zscore.pdf"), width = 10, height = 10)

seurat_list$E16_5 <- hic_object
#!---------------------------------------------------------------------------------------------------------------------------------------------------
hic_object <- seurat_list$E18_5
locs <- Seurat::GetTissueCoordinates(hic_object)[,seq_len(2)]
colnames(locs) <- c("array_row", "array_col")
hic_object@meta.data <- cbind(hic_object@meta.data, locs)

DefaultAssay(hic_object) <- "scAB500k"
VariableFeatures(hic_object) <- rownames(hic_object)

hic_object <- RunBanksy(hic_object, lambda = 0.3, assay = 'scAB500k', slot = 'data', dimx = "array_row", dimy = "array_col", features = 'all', k_geom = 15, assay_name = "BANKSY_scAB500k")
npcs = 30
hic_object <- RunPCA(hic_object, assay = "BANKSY_scAB500k", reduction.name = "pca_banksy", npcs = npcs, features = rownames(hic_object))
hic_object <- RunUMAP(hic_object, reduction = "pca_banksy", dims = 1:npcs, min.dist = 0.3, reduction.name = "banksy_umap", return.model = T, n.neighbors = 24)
hic_object <- FindNeighbors(hic_object, reduction = "pca_banksy", dims = 1:npcs, annoy.metric = "cosine", graph.name = c('banksy_nn', 'banksy_snn'), k.param = 24)
hic_object <- FindClusters(hic_object, resolution = seq(0.1, 2, 0.1), graph.name = "banksy_snn")
wrap_plots(map(seq(0.6, 2, 0.1), function(x) DimPlot(hic_object, reduction = "banksy_umap", group.by = paste0("banksy_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)
wrap_plots(map(seq(0.6, 2, 0.1), function(x) Seurat::SpatialPlot(hic_object, shape = 22, crop = F, group.by = paste0("banksy_snn_res.", x), cols = selected_colors, label = F) + NoLegend()), ncol = 5)
#!----------------------------------------------------------------------------
hic_object <- FindClusters(hic_object, cluster.name = "banksy_cluster", resolution = 1.8, graph.name = "banksy_snn")
hic_object <- FindSubCluster(hic_object, graph.name = "banksy_snn", cluster = 4, subcluster.name = "banksy_sub_cluster", resolution = 0.4)
hic_object$banksy_cluster <- as.character(hic_object$banksy_cluster)
hic_object$banksy_cluster[which(hic_object$banksy_sub_cluster == "4_0")] <- 4
hic_object$banksy_cluster[which(hic_object$banksy_sub_cluster == "4_1")] <- 18
hic_object$banksy_cluster[which(hic_object$banksy_sub_cluster == "4_2")] <- 19
hic_object$banksy_cluster <- factor(hic_object$banksy_cluster, levels= 0:19)
Idents(hic_object) <- "banksy_cluster"
selected_cells <- CellsByIdentities(hic_object)
Seurat::SpatialPlot(hic_object, cells.highlight = selected_cells[setdiff(names(selected_cells), "NA")], cols.highlight = c("#FFFF00", alpha("lightgrey", 0)), crop = FALSE, shape = 22, facet.highlight = T, combine = T, ncol = 7, alpha = NULL) & NoLegend()

hic_object$seurat_clusters <- hic_object$banksy_cluster
hic_object$banksy_cluster <- as.character(hic_object$seurat_clusters)
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(10))] <- "C1"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(15))] <- "C2"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(5))] <- "C3"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(6, 7))] <- "C4"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(3))] <- "C5"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(2))] <- "C6"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(17))] <- "C7"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(4))] <- "C8"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(14))] <- "C9"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(12))] <- "C10"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(11))] <- "C11"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(9, 18))] <- "C12"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(13))] <- "C13"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(0, 19))] <- "C14"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(1))] <- "C15"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(8))] <- "C16"
hic_object$banksy_cluster[which(hic_object$banksy_cluster %in% c(16))] <- "C17"
hic_object$banksy_cluster <- factor(hic_object$banksy_cluster, levels = paste0("C", 1:17))

Idents(hic_object) <- "banksy_cluster"

E18_colors <- c("C1" = "#03FFF4", "C2" = "#0BD3B1", "C3" = "#62CFE8", "C4" = "#1FCC1F", "C5" = "#009E73", "C6" = "#007756", "C7" = "#F2F318", "C8" = "#9E7BFF", "C9" = "#D85FF7", "C10" = "#AA0DFE", "C11" = "#F07F92", "C12" = "#D38B5C", "C13" = "#B5EFB5", "C14" = "#96C885", "C15" = "#CC79A7", "C16" = "#2271A9", "C17" = "#886C00")

DimPlot(hic_object, label = T, cols = E18_colors, reduction = "banksy_umap") + NoLegend() + SetAxes()
Seurat::SpatialPlot(hic_object, cols = E18_colors, crop = FALSE, shape = 22)
hic_object@misc$E18_colors <- E18_colors

seurat_rna <- seurat_RNA$E18_5

DefaultAssay(seurat_rna) <- "SCT"
all_markers <- FindAllMarkers(seurat_rna, only.pos = T, min.pct = 0.1)

de_markers <- all_markers[which(all_markers$avg_log2FC > 1 & all_markers$p_val_adj < 0.01 & all_markers$pct.1 > 0.3),]

marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:100]
	count = count + 1
}
names(marker_genes) <- unique(de_markers$cluster)
DefaultAssay(hic_object) <- "scAB500k"
for(i in names(marker_genes)){
	marker_genes[[i]] <- intersect(rownames(hic_object), unique(gene_infor_500k$gene_region[which(gene_infor_500k$gene_name %in% marker_genes[[i]])]))
	hic_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(hic_object, assay = "scAB500k", slot = "data")[marker_genes[[i]],]))
}

tmp_meta <- as.matrix(hic_object@meta.data[,levels(seurat_rna)])
rownames(tmp_meta) <- paste0(hic_object$banksy_cluster)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- paste0(levels(seurat_rna$Cluster), "_", colnames(tmp_meta))
tmp_meta0 <- tmp_meta0[levels(hic_object$banksy_cluster),]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("brain_20um_analysis/cell_identity_2_density_", names(seurat_list)[3], "_hic_500k_top100_zscore.pdf"), width = 10, height = 10)

seurat_list$E18_5 <- hic_object

qsave(seurat_list, "output/seurat_hic_devp_list_20um.qs")
#!---------------------------------------------------------------------------------------------------------------------------------------------------
seurat_list <- qread("output/seurat_hic_devp_list_20um.qs")
hic_object <- merge(seurat_list[[1]], seurat_list[2:3], add.cell.ids = names(seurat_list))
hic_object$banksy_cluster <- paste0(hic_object$orig.ident, "_", hic_object$banksy_cluster)
hic_object$banksy_cluster <- factor(hic_object$banksy_cluster, levels = c(paste0("E14_5_", levels(seurat_list$E14_5$banksy_cluster)), paste0("E16_5_", levels(seurat_list$E16_5$banksy_cluster)), paste0("E18_5_", levels(seurat_list$E18_5$banksy_cluster))))

hic_object <- JoinLayers(hic_object, assay = "scAB500k")
hic_object <- JoinLayers(hic_object, assay = "scAB500kb_scale")

locs <- c()
for(i in names(hic_object@images)){
	locs <- rbind(locs, Seurat::GetTissueCoordinates(hic_object, image = i)[,seq_len(2)])
}
all(rownames(locs) == colnames(hic_object))
colnames(locs) <- c("sdimy", "sdimx")
hic_object@meta.data <- cbind(hic_object@meta.data, locs)
hic_object$orig.ident <- factor(hic_object$orig.ident, levels = unique(hic_object$orig.ident))

names(E14_colors) <- paste0("E14_5_", names(E14_colors))
names(E16_colors) <- paste0("E16_5_", names(E16_colors))
names(E18_colors) <- paste0("E18_5_", names(E18_colors))
hic_object@misc$banksy_cluster_colors <- c(E14_colors, E16_colors, E18_colors)

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
all(rownames(tmp_embed) == colnames(hic_object))
hic_object[["banksy_per_umap"]] <- CreateDimReducObject(embeddings = tmp_embed, assay = "scAB500k")

hic_object$staggered_sdimy <- 1080 - hic_object$staggered_sdimy
tmp <- hic_object@meta.data[,c("staggered_sdimx", "staggered_sdimy")]
colnames(tmp) <- c("spatial_umap_1", "spatial_umap_2")
hic_object[["spatial_umap"]] <- CreateDimReducObject(embeddings = as.matrix(tmp), assay = "scAB250k")

Idents(hic_object) <- "banksy_cluster"
DimPlot(hic_object, reduction = "banksy_per_umap", cols = hic_object@misc$banksy_cluster_colors)
Seurat::SpatialPlot(hic_object, shape = 22, crop = F, cols = hic_object@misc$banksy_cluster_colors) & NoLegend()
#!----------------------------------------------------------------------------------------------------------------
E14_contacts <- read.csv("contact_heatmap_20um/E14-brain-hic-15_output_matrix_heatmap.csv", row.names = 1)
length(intersect(hic_object$pixel[which(hic_object$orig.ident == "E14_5")], colnames(E14_contacts)))
E14_contacts <- E14_contacts[rev(rownames(E14_contacts)),hic_object$pixel[which(hic_object$orig.ident == "E14_5")]]

E16_contacts <- read.csv("contact_heatmap_20um/E165-brain-hic-20_output_matrix_heatmap.csv", row.names = 1)
length(intersect(hic_object$pixel[which(hic_object$orig.ident == "E16_5")], colnames(E16_contacts)))
E16_contacts <- E16_contacts[rev(rownames(E16_contacts)),hic_object$pixel[which(hic_object$orig.ident == "E16_5")]]

E18_contacts <- read.csv("contact_heatmap_20um/E185-brain-hic-4_output_matrix_heatmap.csv", row.names = 1)
length(intersect(hic_object$pixel[which(hic_object$orig.ident == "E18_5")], colnames(E18_contacts)))
E18_contacts <- E18_contacts[rev(rownames(E18_contacts)),hic_object$pixel[which(hic_object$orig.ident == "E18_5")]]
#!----------------------------------------------------------------------
tmp <- as.numeric(rownames(E14_contacts))
tmp <- round(tmp / 1000, 2)
tmp1 <- paste0(tmp[64:133], "kb")
tmp <- round(tmp / 1000, 2)
tmp2 <- paste0(tmp[1:63], "Mb")
tmp <- c(tmp2, tmp1)
megacontacts <- cbind(E14_contacts, E16_contacts, E18_contacts)
rownames(megacontacts) <- tmp
all(colnames(megacontacts) == hic_object$pixel)
colnames(megacontacts) <- colnames(hic_object)

hic_object[["dist"]] <- CreateAssay5Object(counts = megacontacts)
DefaultAssay(hic_object) <- "scAB500kb_scale"

qsave(hic_object, "output/brain_development_Spatial_HiC_final_tutorial_20um.qs")
#!----------------------------------------------------------------------------------------------------------------
hic_object <- subset(hic_object, subset = banksy_cluster %in% c(paste0("E14_5_", c("C1", "C2", "C3", "C4")), paste0("E16_5_", c("C1", "C2", "C3", "C4")), paste0("E18_5_", c("C1", "C2", "C3", "C4", "C5", "C6"))))

hic_object$Region <- "CP"
hic_object$Region[which(hic_object$banksy_cluster %in% c(paste0("E14_5_", c("C1", "C2", "C3")), paste0("E16_5_", c("C1", "C2")), paste0("E18_5_", c("C1", "C2", "C3"))))] <- "GZ"
hic_object$Region <- factor(hic_object$Region, levels = c("GZ", "CP"))
source("~/scripts/seurat2scanpy/shiny_st.R")
options(browser = "/usr/bin/firefox")
hic_object <- shiny_st(seurat = hic_object, isVisium = F, assay = "scAB500k", image = "E14_5")
hic_object <- subset(hic_object, subset = Region != "removed_cells")
hic_object$Region <- factor(hic_object$Region, levels = c("GZ", "CP"))
hic_object$banksy_cluster <- factor(hic_object$banksy_cluster, levels = intersect(levels(hic_object$banksy_cluster), unique(hic_object$banksy_cluster)))

Seurat::SpatialPlot(hic_object, shape = 22, crop = F, cols = hic_object@misc$banksy_cluster_colors)
ggsave("brain_20um_analysis/HiC_Cortex_Cluster.pdf")
Seurat::SpatialPlot(hic_object, shape = 22, crop = F, group.by = "Region", cols = cluster_colors)
ggsave("brain_20um_analysis/HiC_Cortex_Region.pdf")

Idents(hic_object) <- "banksy_cluster"
qsave(hic_object, "output/brain_development_Spatial_HiC_final_tutorial_cortex_20um.qs")

















