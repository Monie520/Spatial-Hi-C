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

cluster_colors <- c("E14_5" = "#EF5EB4", "E16_5" = "#8A2BD2", "E18_5" = "#0047AB", "Adult15" = "#B312A6", "GZ" = "#66C2A5", "CP" = "#FC8D62",
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
#!------------------------------------------------------------------------------------------------------------------------------------------------------
# perl -lane 'if(/(chr.+?)\t.+?\t(.+?)\t(.+?)\t(.+?).\t(.)\t.\tgene_id\s.(ENSMUSG.+?).\;.+?exon_number\s.(1).\;/){print "$1\t$2\t$3\t$4\t$5\t$6\t$7"}' /date/lvjunjie/Reference/mouse/GRCm38_mm10/ensembl/Mus_musculus.GRCm38.102.chr.gtf > promoter
options(scipen = 999)
gene_infor <- read.table("GRCm38_102_gtf_infor", header = T, sep = "\t")
gene_infor$strand <- "+"

promoter <- read.csv("promoter", header = F)
gene_infor$strand[match(promoter$V6, gene_infor$gene_id)] <- promoter$V5
gene_infor$is_promoter <- NA

gene_infor[which(gene_infor$strand == "+"), 2] <- gene_infor[which(gene_infor$strand == "+"), 2] - 2000
gene_infor[which(gene_infor$strand == "-"), 3] <- gene_infor[which(gene_infor$strand == "-"), 3] + 2000
gene_infor[which(gene_infor$start <= 0), 3] <- 1

bin = 100000

gene_infor$start <- floor((gene_infor$start - 1)/bin) * bin
gene_infor$end <- (floor((gene_infor$end - 1)/bin) + 1) * bin

gene_infor_tmp <- c()
for(i in 1:nrow(gene_infor)){
	if((gene_infor$end[i] - gene_infor$start[i]) > bin){
		tmp <- (gene_infor$end[i] - gene_infor$start[i]) / bin
		if(gene_infor$strand[i] == "+"){
			gene_infor_tmp <- rbind(gene_infor_tmp, c(gene_infor$chr[i], gene_infor$start[i], gene_infor$start[i] + bin, gene_infor$gene_id[i], gene_infor$gene_name[i], gene_infor$gene_type[i], gene_infor$strand[i], "promoter"))
		}else{
			gene_infor_tmp <- rbind(gene_infor_tmp, c(gene_infor$chr[i], gene_infor$start[i], gene_infor$start[i] + bin, gene_infor$gene_id[i], gene_infor$gene_name[i], gene_infor$gene_type[i], gene_infor$strand[i], NA))
		}
		for(j in 2:(tmp-1)){
			gene_infor_tmp <- rbind(gene_infor_tmp, c(gene_infor$chr[i], gene_infor$start[i] + (j - 1) * bin, gene_infor$start[i] + j * bin,  gene_infor$gene_id[i], gene_infor$gene_name[i], gene_infor$gene_type[i], gene_infor$strand[i], NA))
		}
		if(gene_infor$strand[i] == "-"){
			gene_infor_tmp <- rbind(gene_infor_tmp, c(gene_infor$chr[i], gene_infor$start[i] + (tmp - 1) * bin, gene_infor$start[i] + tmp * bin,  gene_infor$gene_id[i], gene_infor$gene_name[i], gene_infor$gene_type[i], gene_infor$strand[i], "promoter"))
		}else{
			gene_infor_tmp <- rbind(gene_infor_tmp, c(gene_infor$chr[i], gene_infor$start[i] + (tmp - 1) * bin, gene_infor$start[i] + tmp * bin,  gene_infor$gene_id[i], gene_infor$gene_name[i], gene_infor$gene_type[i], gene_infor$strand[i], NA))
		}
	}else{
		gene_infor$is_promoter[i] <- "promoter"
		gene_infor_tmp <- rbind(gene_infor_tmp, gene_infor[i,])
	}
}
gene_infor <- gene_infor_tmp
gene_infor$gene_region <- paste0(gene_infor$chr, ":", as.character(gene_infor$start), "-", as.character(gene_infor$end))
write.csv(gene_infor, "GRCm38_102_gtf_infor_hic_100kb_reference", quote = F, row.names = F)
gene_infor_100k <- read.csv("GRCm38_102_gtf_infor_hic_100kb_reference")
gene_infor_250k <- read.csv("GRCm38_102_gtf_infor_hic_250kb_reference")
gene_infor_500k <- read.csv("GRCm38_102_gtf_infor_hic_500kb_reference")
#!--------------------------------------------------------------------------------------------------------------------------------------------------
DefaultAssay(seurat_object) <- "SCT"
seurat_object$time_Region <- paste0(seurat_object$orig.ident, "_", seurat_object$Region)
seurat_object$time_Region <- factor(seurat_object$time_Region, levels = c("E14_5_GZ", "E14_5_CP", "E16_5_GZ", "E16_5_CP", "E18_5_GZ", "E18_5_CP"))
Idents(seurat_object) <- "time_Region"

#pseudo_object <- AggregateExpression(seurat_object, assays = "SCT", return.seurat = T, group.by = c("orig.ident", "Region"))
#pseudo_object$group <- c("E14_5", "E14_5", "E16_5", "E16_5", "E18_5", "E18_5")
#pseudo_object[["pseudobulk"]] <- CreateAssayObject(GetAssayData(pseudo_object, layer = "data"))
#pseudo_object@assays$pseudobulk@scale.data <- GetAssayData(pseudo_object, assay = "SCT", layer = "scale.data")
#pseudo_object$Region <- factor(pseudo_object$Region, levels = c("GZ", "CP"))
#Idents(pseudo_object) <- "Region"

#GZ_CP <- FindMarkers(object = pseudo_object, ident.1 = "CP", ident.2 = "GZ", test.use = "DESeq2")
#GZ_CP <- GZ_CP[which(GZ_CP$p_val_adj < 0.05),]
#GZ_CP <- GZ_CP[-grep("^Gm[0-9]|Rik$|Rik[1-9]$", rownames(GZ_CP)),]


E14_5_GZ_CP <- FindMarkers(seurat_object, ident.1 = "E14_5_CP", ident.2 = "E14_5_GZ", min.pct = 0.1)
E16_5_GZ_CP <- FindMarkers(seurat_object, ident.1 = "E16_5_CP", ident.2 = "E16_5_GZ", min.pct = 0.1)
E18_5_GZ_CP <- FindMarkers(seurat_object, ident.1 = "E18_5_CP", ident.2 = "E18_5_GZ", min.pct = 0.1)

E14_5_GZ_CP <- E14_5_GZ_CP %>% filter(abs(avg_log2FC) >= 0.7, p_val_adj < 0.05)
E16_5_GZ_CP <- E16_5_GZ_CP %>% filter(abs(avg_log2FC) >= 0.7, p_val_adj < 0.05)
E18_5_GZ_CP <- E18_5_GZ_CP %>% filter(abs(avg_log2FC) >= 0.7, p_val_adj < 0.05)

E14_5_GZ_CP <- E14_5_GZ_CP[-which(E14_5_GZ_CP$pct.1 < 0.25 & E14_5_GZ_CP$pct.2 < 0.25),]
E16_5_GZ_CP <- E16_5_GZ_CP[-which(E16_5_GZ_CP$pct.1 < 0.25 & E16_5_GZ_CP$pct.2 < 0.25),]
E18_5_GZ_CP <- E18_5_GZ_CP[-which(E18_5_GZ_CP$pct.1 < 0.25 & E18_5_GZ_CP$pct.2 < 0.25),]

E14_5_GZ_CP <- E14_5_GZ_CP[-grep("^Gm[0-9]|Rik$|Rik[1-9]$", rownames(E14_5_GZ_CP)),]
E16_5_GZ_CP <- E16_5_GZ_CP[-grep("^Gm[0-9]|Rik$|Rik[1-9]$", rownames(E16_5_GZ_CP)),]
E18_5_GZ_CP <- E18_5_GZ_CP[-grep("^Gm[0-9]|Rik$|Rik[1-9]$", rownames(E18_5_GZ_CP)),]

features <- unique(c(rownames(E14_5_GZ_CP), rownames(E16_5_GZ_CP), rownames(E18_5_GZ_CP)))

ht <- SCP::GroupHeatmap(seurat_object, features = features, slot = "counts", assay = "SCT", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palette = "viridis", nlabel = 0, show_row_names = F, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot

GZ_CP <- ht$feature_metadata
GZ_CP <- GZ_CP[order(GZ_CP$index),]
zscore_matrix <- ht$matrix_list$Region
colnames(zscore_matrix) <- c("E14_GZ", "E16_GZ", "E18_GZ", "E14_CP", "E16_CP", "E18_CP")
all(rownames(zscore_matrix) == GZ_CP$gene)
zscore_matrix <- zscore_matrix[GZ_CP$gene,]

GZ_CP$final_cluster <- NA
bin = 0.5
for(i in 1:nrow(GZ_CP)){
	tmp <- max(zscore_matrix[i,])
	tmp_names <- colnames(zscore_matrix)[which(zscore_matrix[i,] == tmp)]
	tmp_sec <- max(zscore_matrix[i,setdiff(colnames(zscore_matrix), tmp_names)])
	tmp_sec_names <- colnames(zscore_matrix)[which(zscore_matrix[i,] == tmp_sec)]
	if(tmp > (tmp_sec + bin) | tmp > (tmp_sec * 2)){
		GZ_CP$final_cluster[i] <- tmp_names
	}else if(all(c(tmp_names, tmp_sec_names) %in% c("E14_GZ", "E16_GZ"))){
		if(tmp > (tmp_sec  + zscore_matrix[i,"E18_GZ"] + bin)){
			GZ_CP$final_cluster[i] <- tmp_names
		}else if(tmp_sec > (zscore_matrix[i,"E18_GZ"] + bin)){
			GZ_CP$final_cluster[i] <- "E14_E16_GZ"
		}else{
			GZ_CP$final_cluster[i] <- "E14_E16_E18_GZ"
		}
	}else if(all(c(tmp_names, tmp_sec_names) %in% c("E16_GZ", "E18_GZ"))){
		if(tmp > (tmp_sec  + zscore_matrix[i,"E14_GZ"] + bin)){
			GZ_CP$final_cluster[i] <- tmp_names
		}else if(tmp_sec > (zscore_matrix[i,"E14_GZ"] + bin)){
			GZ_CP$final_cluster[i] <- "E16_E18_GZ"
		}else{
			GZ_CP$final_cluster[i] <- "E14_E16_E18_GZ"
		}
	}else if(all(c(tmp_names, tmp_sec_names) %in% c("E14_GZ", "E18_GZ"))){
			if(tmp > (tmp_sec  + zscore_matrix[i,"E16_GZ"] + bin)){
				GZ_CP$final_cluster[i] <- tmp_names
			}else{
				GZ_CP$final_cluster[i] <- "E14_E16_E18_GZ"
			}
	}else if(all(c(tmp_names, tmp_sec_names) %in% c("E14_CP", "E16_CP"))){
		if(tmp > (tmp_sec  + zscore_matrix[i,"E18_CP"] + bin)){
			GZ_CP$final_cluster[i] <- tmp_names
		}else if(tmp_sec > (zscore_matrix[i,"E18_CP"] + bin)){
			GZ_CP$final_cluster[i] <- "E14_E16_CP"
		}else{
			GZ_CP$final_cluster[i] <- "E14_E16_E18_CP"
		}
	}else if(all(c(tmp_names, tmp_sec_names) %in% c("E16_CP", "E18_CP"))){
		if(tmp > (tmp_sec  + zscore_matrix[i,"E14_CP"] + bin)){
			GZ_CP$final_cluster[i] <- tmp_names
		}else if(tmp_sec > (zscore_matrix[i,"E14_CP"] + bin)){
			GZ_CP$final_cluster[i] <- "E16_E18_CP"
		}else{
			GZ_CP$final_cluster[i] <- "E14_E16_E18_CP"
		}
	}else if(all(c(tmp_names, tmp_sec_names) %in% c("E14_CP", "E18_CP"))){
			if(tmp > (tmp_sec  + zscore_matrix[i,"E16_CP"] + bin)){
				GZ_CP$final_cluster[i] <- tmp_names
			}else{
				GZ_CP$final_cluster[i] <- "E14_E16_E18_CP"
			}
	}else{
		GZ_CP$final_cluster[i] <- "E14_E16_E18"
	}
}

GZ_CP$final_cluster <- factor(GZ_CP$final_cluster, levels = c("E14_GZ", "E14_E16_GZ", "E16_GZ", "E16_E18_GZ", "E18_GZ", "E14_E16_E18_GZ", "E14_E16_E18", "E14_CP", "E14_E16_CP", "E16_CP", "E16_E18_CP", "E18_CP", "E14_E16_E18_CP"))

ht <- SCP::GroupHeatmap(seurat_object, features = GZ_CP$gene, feature_split = GZ_CP$final_cluster, slot = "counts", assay = "SCT", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palette = "viridis", nlabel = 0, show_row_names = F, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot

GZ_CP$Clusters <- as.character(GZ_CP$final_cluster)
GZ_CP$Clusters[which(GZ_CP$final_cluster == "E14_GZ")] <- "C1"
GZ_CP$Clusters[which(GZ_CP$final_cluster == "E14_E16_GZ")] <- "C2"
GZ_CP$Clusters[which(GZ_CP$final_cluster == "E16_GZ")] <- "C3"
GZ_CP$Clusters[which(GZ_CP$final_cluster == "E16_E18_GZ")] <- "C4"
GZ_CP$Clusters[which(GZ_CP$final_cluster == "E18_GZ")] <- "C5"
GZ_CP$Clusters[which(GZ_CP$final_cluster == "E14_E16_E18_GZ")] <- "C6"
GZ_CP$Clusters[which(GZ_CP$final_cluster == "E14_E16_E18")] <- "C7"
GZ_CP$Clusters[which(GZ_CP$final_cluster == "E14_CP")] <- "C8"
GZ_CP$Clusters[which(GZ_CP$final_cluster == "E14_E16_CP")] <- "C9"
GZ_CP$Clusters[which(GZ_CP$final_cluster == "E16_CP")] <- "C10"
GZ_CP$Clusters[which(GZ_CP$final_cluster == "E16_E18_CP")] <- "C11"
GZ_CP$Clusters[which(GZ_CP$final_cluster == "E18_CP")] <- "C12"
GZ_CP$Clusters[which(GZ_CP$final_cluster == "E14_E16_E18_CP")] <- "C13"

GZ_CP$Clusters <- factor(GZ_CP$Clusters, levels = paste0("C", 1:13))

key_markers <- c("Bcl11b", "Celf4", "Cux1", "Cux2", "Egfr", "Emx1", "Eomes", "Fabp7", "Fezf2", "Foxp1", "Lhx2", "Mef2c", "Meis2", "Nr2e1", "Nrp1", "Ntrk2", "Pax6", "Pex5l", "Pou3f2", "Pou3f3", "Ptprd", "Ptprk", "Slc1a3", "Sox2", "Sox5", "Tafa1", "Tafa2", "Tbr1", "Tle4", "Zbtb20", "Satb2", "Cadm2", "Dlgap1", "Frmd3", "Grin2a", "Etv6", "L3mbtl4", "Map9", "Nlgn1", "Nrg1", "Sh3bgrl2", "Slc12a7", "Rasgrf2", "Slc8a1", "Unc5d", "Wwtr1", "Neurod1", "Rgs20", "Notch1", "Notch2", "Nrxn1", "Sv2b", "Slc24a2", "Robo1", "Robo2", "Plxna4", "Slc1a3", "Apoe", "Gpr85", "Id2")

ht <- SCP::GroupHeatmap(seurat_object, features = GZ_CP$gene, feature_split = GZ_CP$Clusters, slot = "counts", assay = "SCT", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palette = "viridis", nlabel = 0, show_row_names = F, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE), features_label = key_markers)
ht$plot

pdf("GZ_CP_RNA_Clusters.pdf", width = 8, height = 12)
ht$plot
dev.off()

ht <- SCP::GroupHeatmap(seurat_object, features = GZ_CP$gene, feature_split = GZ_CP$Clusters, slot = "counts", assay = "SCT", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c("#0099CC", "white", "#CC0033"), nlabel = 0, show_row_names = F, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE), features_label = key_markers)
ht$plot

pdf("GZ_CP_RNA_Clusters1.pdf", width = 8, height = 12)
ht$plot
dev.off()

GZ_CP$index <- ht$feature_metadata$index[match(GZ_CP$gene, ht$feature_metadata$features)]
GZ_CP <- GZ_CP[order(GZ_CP$index),]
zscore_matrix <- ht$matrix_list$Region
colnames(zscore_matrix) <- c("E14_GZ", "E16_GZ", "E18_GZ", "E14_CP", "E16_CP", "E18_CP")
zscore_matrix <- zscore_matrix[GZ_CP$gene,]
all(rownames(zscore_matrix) == GZ_CP$gene)
GZ_CP <- cbind(GZ_CP, zscore_matrix)

write.csv(GZ_CP, "GZ_CP_meta_data.csv", quote = F, row.names = F)
#!---------------------------------------------------------------------------------------------------------------------------------------------
dir.create("Cortex_GZ_CP_Region_features")
for(i in unique(GZ_CP$Clusters)){
        gene_list <- as.character(GZ_CP[which(GZ_CP$Clusters == i),'gene'])
        enrichRes <- enrichAnalysis(genelist = gene_list, geneType = "SYMBOL", species = "mmu", database = c("go", "kegg", "MSigDb"), GO.model = "BP", MSigDb.signature = "H", sampleName = paste0("enrich_", i), minGSSize = 5, outpath = "Cortex_GZ_CP_Region_features")
}

cell_cluster_colors <- c(C1 = "#C22B86", C2 = "#F2639E", C3 = "#FBB5BB", C4 = "#5A5A89", C5 = "#1F73A7", C6 = "#68B2D3", C7 = "#4F928A", C8 = "#1C7B44", C9 = "#74C174", C10 = "#DCCF79", C11 = "#E77306", C12 = "#F28F2C", C13 = "#FFBB61")

for(i in unique(GZ_CP$Clusters)){
	enrich_go_terms <- read.csv(paste0("Cortex_GZ_CP_Region_features/enrich_", i, "_enrich_GeneOntology_BP_terms_selected.csv"))
	enrich_go_terms$logFDR <- -log10(enrich_go_terms$p.adjust)
	enrich_go_terms$Description <- factor(enrich_go_terms$Description, levels = rev(enrich_go_terms$Description))
	enrich_go_terms$Clusters <- i
	if(nrow(enrich_go_terms) > 10){
		height = 5
		width = 7
	}else if(nrow(enrich_go_terms) > 8){
		height = 4
		width = 7
	}else if(nrow(enrich_go_terms) > 6){
		height = 3
		width = 5
	}else if(nrow(enrich_go_terms) > 4){
		height = 3
		width = 5
	}else{
		height = 2
		width = 5
	}

	p0 <- ggplot(enrich_go_terms, aes_string(x = "logFDR", y = "Description", fill = "Clusters")) + theme_classic() + geom_col() + scale_fill_manual(values = cell_cluster_colors) + NoLegend() +
		geom_text(data = enrich_go_terms, aes(x = 0.05, y = Description, label = Description), size = 4, hjust = 0) +
	  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line = element_line(color = 'black', size = 1.1), axis.text = element_text(size = 12, color = "black"), plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_text(size = 16, color = "black", face = "bold")) + xlab("-log10(p.adjust)") + labs(title = i)
     print(p0)
     ggsave(paste0("Cortex_GZ_CP_Region_features/enrich_", i, "_enrich_GeneOntology_BP_terms_selected.pdf"), width = width, height = height)
}

cluster_genes <- list()
for(i in levels(GZ_CP$Clusters)){
	cluster_genes[[i]] <- GZ_CP$gene[which(GZ_CP$Clusters == i)]
}

DefaultAssay(hic_object) <- "scAB250k"
for(i in names(cluster_genes)){
	cluster_genes[[i]] <- intersect(rownames(hic_object), unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name %in% cluster_genes[[i]])]))
	hic_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(hic_object, assay = "scAB250k", slot = "data")[cluster_genes[[i]],]))
}

hic_object$C1_normalized <- normalize_to_01(hic_object$C1)
hic_object$C2_normalized <- normalize_to_01(hic_object$C2)
hic_object$C3_normalized <- normalize_to_01(hic_object$C3)
hic_object$C4_normalized <- normalize_to_01(hic_object$C4)
hic_object$C5_normalized <- normalize_to_01(hic_object$C5)
hic_object$C6_normalized <- normalize_to_01(hic_object$C6)
hic_object$C7_normalized <- normalize_to_01(hic_object$C7)
hic_object$C8_normalized <- normalize_to_01(hic_object$C8)
hic_object$C9_normalized <- normalize_to_01(hic_object$C9)
hic_object$C10_normalized <- normalize_to_01(hic_object$C10)
hic_object$C11_normalized <- normalize_to_01(hic_object$C11)
hic_object$C12_normalized <- normalize_to_01(hic_object$C12)
hic_object$C13_normalized <- normalize_to_01(hic_object$C13)

VlnPlot(hic_object, features = paste0(setdiff(names(cluster_genes), "C7"), "_normalized"), pt.size = 0, group.by = "orig.ident", split.by = "Region", cols = c("#66C2A5", "#FC8D62"), ncol = 4) & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.3, position = position_dodge(0.9)) & NoLegend() & ylab("normalized mean scAB values") & xlab("")
VlnPlot(hic_object, features = setdiff(names(cluster_genes), "C7"), pt.size = 0, group.by = "orig.ident", split.by = "Region", cols = c("#66C2A5", "#FC8D62"), ncol = 4) & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.3, position = position_dodge(0.9)) & NoLegend() & ylab("mean scAB values") & xlab("")
pdf("GZ_CP_RNA_Clusters_HiC_mean_scAB_values.pdf", width = 4)
for(i in names(cluster_genes)){
	p0 <- VlnPlot(hic_object, features = paste0(i, "_normalized"), pt.size = 0, group.by = "orig.ident", split.by = "Region", cols = c("#66C2A5", "#FC8D62")) & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.3, position = position_dodge(0.9)) & NoLegend() & ylab("normalized mean scAB values") & xlab("")
	p1 <- VlnPlot(hic_object, features = i, pt.size = 0, group.by = "orig.ident", split.by = "Region", cols = c("#66C2A5", "#FC8D62")) & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.3, position = position_dodge(0.9)) & NoLegend() & ylab("normalized mean scAB values") & xlab("")
	print(p0 + p1)
}
dev.off()

P_values <- c()
for(i in names(cluster_genes)){
	for(j in levels(hic_object$orig.ident)){
		GZ <- hic_object@meta.data[which(hic_object$orig.ident == j & hic_object$Region == "GZ"),paste0(i)]
		CP <- hic_object@meta.data[which(hic_object$orig.ident == j & hic_object$Region == "CP"),paste0(i)]
		P_valves_t <- t.test(GZ, CP)[3]
		P_valves_w <- wilcox.test(GZ, CP)[3]
		P_values <- rbind(P_values, c(i, j, P_valves_t, P_valves_w))
	}
}
colnames(P_values) <- c("Cluster", "group", "t_test_P_value", "wilcox_test_P_value")
write.csv(P_values, "GZ_CP_RNA_Clusters_HiC_mean_scAB_values_P_values.csv", quote = F, row.names = F)

VlnPlot(hic_object, features = setdiff(names(cluster_genes), "C7"), pt.size = 0, group.by = "banksy_enrich", split.by = "orig.ident", cols = c("#34D916", "#00D4E6", "#1E90FF")) & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.3, position = position_dodge(0.9)) & NoLegend()
VlnPlot(hic_object, features = paste0(names(core_genes), "_normalized"), pt.size = 0, group.by = "orig.ident", split.by = "banksy_cluster", cols = hic_object@misc$banksy_cluster_colors[match(levels(hic_object$banksy_cluster), names(hic_object@misc$banksy_cluster_colors))]) & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.3, position = position_dodge(0.9), show.legend = F)
VlnPlot(hic_object, features = "MigN_normalized", pt.size = 0, group.by = "orig.ident", split.by = "banksy_cluster", cols = hic_object@misc$banksy_cluster_colors[match(levels(hic_object$banksy_cluster), names(hic_object@misc$banksy_cluster_colors))]) & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.3, position = position_dodge(0.9), show.legend = F)

mean_values <- c()
mean_values <- c(mean_values, mean(seurat_object$MigN_normalized[which(seurat_object$orig.ident == "E14_5")]))
mean_values <- c(mean_values, mean(seurat_object$MigN_normalized[which(seurat_object$orig.ident == "E16_5")]))
mean_values <- c(mean_values, mean(seurat_object$MigN_normalized[which(seurat_object$orig.ident == "E18_5")]))
mean_values <- c(mean_values, mean(hic_object$MigN_normalized[which(hic_object$orig.ident == "E14_5")]))
mean_values <- c(mean_values, mean(hic_object$MigN_normalized[which(hic_object$orig.ident == "E16_5")]))
mean_values <- c(mean_values, mean(hic_object$MigN_normalized[which(hic_object$orig.ident == "E18_5")]))

meta_data <- data.frame(class = rep(c("RNA", "HiC"), each = 3), time = rep(c("E14_5", "E16_5", "E18_5"), 2), mean_values = mean_values)

ggplot(meta_data, aes(x = time, y = mean_values, group = class)) + geom_line(aes(color = class), size = 1, show.legend = T) + geom_point(aes(color = class), shape = 19, size = 3, show.legend = F) +
                        scale_color_manual(values = c("darkred", "darkgreen")) + ggtitle("Identified genes/region scores during cell transition") + xlab("time") + ylab("normalized gene expression scores or mean scAB values") + ylim(c(0, 1))
#!------------------------------------------------------------------------------------------------------------------------------------------------
DefaultAssay(seurat_object) <- "spARC_SCT"
DefaultAssay(hic_object) <- "spARC_scAB250kb_scale"
dir.create("Cortex_GZ_CP_Region_features/features")
for(i in GZ_CP$gene){
	try({
		p0 <- Seurat::SpatialPlot(seurat_object, features = i, shape = 22, crop = F) & scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")))
		p1 <- Seurat::SpatialPlot(hic_object, shape = 22, crop = F, features = unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name == i & gene_infor_250k$is_promoter == "promoter")])) & scale_fill_gradientn(values = seq(0, 1, 0.01), colours = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"))
		print(p0/p1)
		ggsave(paste0("Cortex_GZ_CP_Region_features/features/", i, ".png"), width = 24, height = 24)
		p0 <- FeaturePlot(seurat_object, reduction = "spatial_umap", features = i) & scale_colour_gradientn(colours = paletteContinuous(set = "solarExtra"))
		p1 <- FeaturePlot(hic_object, features = unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name == i & gene_infor_250k$is_promoter == "promoter")]), reduction = "spatial_umap") & scale_colour_gradientn(colours = paletteContinuous(set = "solarExtra"))
		print(p0/p1)
		ggsave(paste0("Cortex_GZ_CP_Region_features/features/", i, "_0.png"), width = 24, height = 24)
	})
}

DefaultAssay(seurat_object) <- "SCT"
DefaultAssay(hic_object) <- "scAB250kb_scale"
for(i in GZ_CP$gene){
	try({
		p0 <- VlnPlot(seurat_object, features = i, pt.size = 0, group.by = "orig.ident", split.by = "Region", cols = c("#00D4E6", "#00998F"))
		p1 <- VlnPlot(hic_object, pt.size = 0, features = unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name == i & gene_infor_250k$is_promoter == "promoter")]), group.by = "orig.ident", split.by = "Region", cols = c("#00D4E6", "#00998F")) + ylab("scAB values") & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.3, position = position_dodge(0.9), show.legend = F)
		print(p0 + p1)
		ggsave(paste0("Cortex_GZ_CP_Region_features/features/", i, "_vlnplot.png"), width = 12, height = 14)
	})
}
#!------------------------------------------------------------------------------------------------------------------------------------------------
DefaultAssay(hic_object) <- "scAB250kb_scale"
Idents(hic_object) <- "time_Region"

scAB_matrix <- GetAssayData(hic_object, assay = "scAB250kb", slot = "counts")
scAB_matrix <- t(scAB_matrix)
rownames(scAB_matrix) <- paste0(hic_object$time_Region)
scAB_matrix <- limma::avereps(scAB_matrix)
tmp_object <- CreateSeuratObject(counts = NULL, data = t(scAB_matrix))
tmp_object <- FindVariableFeatures(tmp_object, assay = "RNA", selection.method = "mvp")
scAB_matrix <- t(scAB_matrix)[VariableFeatures(tmp_object),]
#scAB_matrix <- scale(scAB_matrix)
dist_mat <- as.dist(1 - cor(t(scAB_matrix)))
hc_result <- hclust(dist_mat, method = "average")
hc_order <- order.dendrogram(as.dendrogram(hc_result))
scAB_matrix <- scAB_matrix[hc_order,]
scAB_matrix <- scale(scAB_matrix)
ComplexHeatmap::Heatmap(scale(scAB_matrix), cluster_rows = F, show_row_names = F)

E14_5_GZ_CP <- FindMarkers(hic_object, ident.1 = "E14_5_CP", ident.2 = "E14_5_GZ", min.pct = 0.1, test.use = 'LR', fc.name = "avg_diff")
E16_5_GZ_CP <- FindMarkers(hic_object, ident.1 = "E16_5_CP", ident.2 = "E16_5_GZ", min.pct = 0.1, test.use = 'LR', fc.name = "avg_diff")
E18_5_GZ_CP <- FindMarkers(hic_object, ident.1 = "E18_5_CP", ident.2 = "E18_5_GZ", min.pct = 0.1, test.use = 'LR', fc.name = "avg_diff")

E14_5_GZ_CP <- E14_5_GZ_CP %>% filter(abs(avg_diff) >= 0.3, p_val_adj < 0.05)
E16_5_GZ_CP <- E16_5_GZ_CP %>% filter(abs(avg_diff) >= 0.3, p_val_adj < 0.05)
E18_5_GZ_CP <- E18_5_GZ_CP %>% filter(abs(avg_diff) >= 0.3, p_val_adj < 0.05)

features <- unique(c(rownames(E14_5_GZ_CP), rownames(E16_5_GZ_CP), rownames(E18_5_GZ_CP)))

scAB_matrix <- GetAssayData(hic_object, assay = "scAB250kb", slot = "counts")[features,]
scAB_matrix <- t(scAB_matrix)
rownames(scAB_matrix) <- paste0(hic_object$time_Region)
scAB_matrix <- t(limma::avereps(scAB_matrix))
mean_E14 <- rowMeans(scAB_matrix[,1:2])
mean_E16 <- rowMeans(scAB_matrix[,3:4])
mean_E18 <- rowMeans(scAB_matrix[,5:6])
scAB_matrix_center <- cbind(scAB_matrix[,1:2] - mean_E14, scAB_matrix[,3:4] - mean_E16, scAB_matrix[,5:6] - mean_E18)
tmp_object <- CreateSeuratObject(counts = scAB_matrix, data = scAB_matrix_center)
tmp_object[["scAB250kb"]] <- CreateAssayObject(data = scAB_matrix_center)
tmp_object$Region <- c("GZ", "CP", "CP", "GZ", "CP", "GZ")
tmp_object$orig.ident <- c("E14_5", "E14_5", "E16_5", "E16_5", "E18_5", "E18_5")
tmp_object$Region <- factor(tmp_object$Region, levels = c("GZ", "CP"))
tmp_object$orig.ident <- factor(tmp_object$orig.ident, levels = c("E14_5", "E16_5", "E18_5"))

ht <- SCP::GroupHeatmap(tmp_object, features = rownames(tmp_object), n_split = 2, assay = "scAB250kb", slot = "data", exp_method = "raw", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(tmp_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot

GZ_CP_HiC <- ht$feature_metadata
GZ_CP_HiC <- GZ_CP_HiC[order(GZ_CP_HiC$index),]
scAB_matrix_center <- scAB_matrix_center[GZ_CP_HiC$features,]
colnames(scAB_matrix_center) <- c("E14_GZ", "E14_CP", "E16_CP", "E16_GZ", "E18_CP", "E18_GZ")

GZ_CP_HiC$final_cluster <- NA
bin = 0.015
for(i in 1:nrow(GZ_CP_HiC)){
	tmp <- max(scAB_matrix_center[i,])
	tmp_names <- colnames(scAB_matrix_center)[which(scAB_matrix_center[i,] == tmp)]
	tmp_sec <- max(scAB_matrix_center[i,setdiff(colnames(scAB_matrix_center), tmp_names)])
	tmp_sec_names <- colnames(scAB_matrix_center)[which(scAB_matrix_center[i,] == tmp_sec)]
	if(tmp > (tmp_sec + bin)){
		GZ_CP_HiC$final_cluster[i] <- tmp_names
	}else if(all(c(tmp_names, tmp_sec_names) %in% c("E14_GZ", "E16_GZ"))){
		if(tmp > (tmp_sec  + scAB_matrix_center[i,"E18_GZ"] + bin)){
			GZ_CP_HiC$final_cluster[i] <- tmp_names
		}else if(tmp_sec > (scAB_matrix_center[i,"E18_GZ"] + bin)){
			GZ_CP_HiC$final_cluster[i] <- "E14_E16_GZ"
		}else{
			GZ_CP_HiC$final_cluster[i] <- "E14_E16_E18_GZ"
		}
	}else if(all(c(tmp_names, tmp_sec_names) %in% c("E16_GZ", "E18_GZ"))){
		if(tmp > (tmp_sec  + scAB_matrix_center[i,"E14_GZ"] + bin)){
			GZ_CP_HiC$final_cluster[i] <- tmp_names
		}else if(tmp_sec > (scAB_matrix_center[i,"E14_GZ"] + bin)){
			GZ_CP_HiC$final_cluster[i] <- "E16_E18_GZ"
		}else{
			GZ_CP_HiC$final_cluster[i] <- "E14_E16_E18_GZ"
		}
	}else if(all(c(tmp_names, tmp_sec_names) %in% c("E14_GZ", "E18_GZ"))){
			if(tmp > (tmp_sec  + scAB_matrix_center[i,"E16_GZ"] + bin)){
				GZ_CP_HiC$final_cluster[i] <- tmp_names
			}else{
				GZ_CP_HiC$final_cluster[i] <- "E14_E16_E18_GZ"
			}
	}else if(all(c(tmp_names, tmp_sec_names) %in% c("E14_CP", "E16_CP"))){
		if(tmp > (tmp_sec  + scAB_matrix_center[i,"E18_CP"] + bin)){
			GZ_CP_HiC$final_cluster[i] <- tmp_names
		}else if(tmp_sec > (scAB_matrix_center[i,"E18_CP"] + bin)){
			GZ_CP_HiC$final_cluster[i] <- "E14_E16_CP"
		}else{
			GZ_CP_HiC$final_cluster[i] <- "E14_E16_E18_CP"
		}
	}else if(all(c(tmp_names, tmp_sec_names) %in% c("E16_CP", "E18_CP"))){
		if(tmp > (tmp_sec  + scAB_matrix_center[i,"E14_CP"] + bin)){
			GZ_CP_HiC$final_cluster[i] <- tmp_names
		}else if(tmp_sec > (scAB_matrix_center[i,"E14_CP"] + bin)){
			GZ_CP_HiC$final_cluster[i] <- "E16_E18_CP"
		}else{
			GZ_CP_HiC$final_cluster[i] <- "E14_E16_E18_CP"
		}
	}else if(all(c(tmp_names, tmp_sec_names) %in% c("E14_CP", "E18_CP"))){
			if(tmp > (tmp_sec  + scAB_matrix_center[i,"E16_CP"] + bin)){
				GZ_CP_HiC$final_cluster[i] <- tmp_names
			}else{
				GZ_CP_HiC$final_cluster[i] <- "E14_E16_E18_CP"
			}
	}else{
		GZ_CP_HiC$final_cluster[i] <- "E14_E16_E18"
	}
}

GZ_CP_HiC$final_cluster <- factor(GZ_CP_HiC$final_cluster, levels = c("E14_GZ", "E14_E16_GZ", "E16_GZ", "E16_E18_GZ", "E18_GZ", "E14_E16_E18_GZ", "E14_E16_E18", "E14_CP", "E14_E16_CP", "E16_CP", "E16_E18_CP", "E18_CP", "E14_E16_E18_CP"))

ht <- SCP::GroupHeatmap(tmp_object, features = GZ_CP_HiC$features, feature_split = GZ_CP_HiC$final_cluster, assay = "scAB250kb", slot = "data", exp_method = "raw", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(tmp_object$Region), names(cluster_colors))]), cell_split_palcolor = list(paletteContinuous(set = "sambaNight")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot

GZ_CP_HiC$Clusters <- as.character(GZ_CP_HiC$final_cluster)
GZ_CP_HiC$Clusters[which(GZ_CP_HiC$final_cluster == "E14_GZ")] <- "C1"
GZ_CP_HiC$Clusters[which(GZ_CP_HiC$final_cluster == "E14_E16_GZ")] <- "C2"
GZ_CP_HiC$Clusters[which(GZ_CP_HiC$final_cluster == "E16_GZ")] <- "C3"
GZ_CP_HiC$Clusters[which(GZ_CP_HiC$final_cluster == "E16_E18_GZ")] <- "C4"
GZ_CP_HiC$Clusters[which(GZ_CP_HiC$final_cluster == "E18_GZ")] <- "C5"
GZ_CP_HiC$Clusters[which(GZ_CP_HiC$final_cluster == "E14_E16_E18_GZ")] <- "C6"
GZ_CP_HiC$Clusters[which(GZ_CP_HiC$final_cluster == "E14_E16_E18")] <- "C7"
GZ_CP_HiC$Clusters[which(GZ_CP_HiC$final_cluster == "E14_CP")] <- "C8"
GZ_CP_HiC$Clusters[which(GZ_CP_HiC$final_cluster == "E14_E16_CP")] <- "C9"
GZ_CP_HiC$Clusters[which(GZ_CP_HiC$final_cluster == "E16_CP")] <- "C10"
GZ_CP_HiC$Clusters[which(GZ_CP_HiC$final_cluster == "E16_E18_CP")] <- "C11"
GZ_CP_HiC$Clusters[which(GZ_CP_HiC$final_cluster == "E18_CP")] <- "C12"
GZ_CP_HiC$Clusters[which(GZ_CP_HiC$final_cluster == "E14_E16_E18_CP")] <- "C13"

GZ_CP_HiC$Clusters <- factor(GZ_CP_HiC$Clusters, levels = paste0("C", 1:13))

ht <- SCP::GroupHeatmap(tmp_object, features = GZ_CP_HiC$features, feature_split = GZ_CP_HiC$Clusters, assay = "scAB250kb", slot = "data", exp_method = "raw", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(tmp_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot

pdf("GZ_CP_HiC_Clusters.pdf", width = 8, height = 12)
ht$plot
dev.off()

ht <- SCP::GroupHeatmap(tmp_object, features = GZ_CP_HiC$features, feature_split = GZ_CP_HiC$Clusters, assay = "scAB250kb", slot = "data", exp_method = "raw", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(tmp_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c(rev(col_scale_mako)[4:10], col_scale_acton[1:7]), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot

pdf("GZ_CP_HiC_Clusters1.pdf", width = 8, height = 12)
ht$plot
dev.off()

GZ_CP_HiC$index <- ht$feature_metadata$index[match(GZ_CP_HiC$features, ht$feature_metadata$features)]
GZ_CP_HiC <- GZ_CP_HiC[order(GZ_CP_HiC$index),]
scAB_matrix_center <- scAB_matrix_center[GZ_CP_HiC$features,]
colnames(scAB_matrix_center) <- c("E14_GZ", "E14_CP", "E16_CP", "E16_GZ", "E18_CP", "E18_GZ")
colnames(scAB_matrix_center) <- paste0(colnames(scAB_matrix_center), "_center")
scAB_matrix <- scAB_matrix[GZ_CP_HiC$features,]
GZ_CP_HiC <- cbind(GZ_CP_HiC, scAB_matrix_center, scAB_matrix)

write.csv(GZ_CP_HiC, "GZ_CP_HiC_meta_data.csv", quote = F, row.names = F)

dir.create("Cortex_HiC_GZ_CP_Region_features")
HiC_cluster_genes <- list()
for(i in levels(GZ_CP_HiC$Clusters)){
	HiC_cluster_genes[[i]] <- GZ_CP_HiC$features[which(GZ_CP_HiC$Clusters == i)]
	HiC_cluster_genes[[i]] <- unique(gene_infor_250k$gene_name[which(gene_infor_250k$gene_region %in% HiC_cluster_genes[[i]] & gene_infor_250k$gene_type == "protein_coding")])
	if(length(grep("Gm[0-9]|Rik$|Rik[0-9]$", HiC_cluster_genes[[i]])) > 0){
		HiC_cluster_genes[[i]] <- HiC_cluster_genes[[i]][-grep("Gm[0-9]|Rik$|Rik[0-9]$", HiC_cluster_genes[[i]])]
	}
}
for(i in names(HiC_cluster_genes)){
        gene_list <- HiC_cluster_genes[[i]]
        enrichRes <- enrichAnalysis(genelist = gene_list, geneType = "SYMBOL", species = "mmu", database = c("go", "kegg", "MSigDb"), GO.model = "BP", MSigDb.signature = "H", sampleName = paste0("enrich_", i), minGSSize = 5, outpath = "Cortex_HiC_GZ_CP_Region_features")
}
hic_rna_features <- c()
for(i in names(HiC_cluster_genes)){
        gene_list <- data.frame(gene = HiC_cluster_genes[[i]], Clusters = i)
        hic_rna_features <- rbind(hic_rna_features, gene_list)
}

hic_rna_features$Clusters <- factor(hic_rna_features$Clusters, levels = paste0("C", 1:13))

key_markers <- c("Bcl11b", "Celf4", "Cux1", "Cux2", "Egfr", "Emx1", "Eomes", "Fabp7", "Fezf2", "Foxp1", "Lhx2", "Mef2c", "Meis2", "Nr2e1", "Nrp1", "Ntrk2", "Pax6", "Pex5l", "Pou3f2", "Pou3f3", "Ptprd", "Ptprk", "Slc1a3", "Sox2", "Sox5", "Tafa1", "Tafa2", "Tbr1", "Tle4", "Zbtb20", "Satb2", "Cadm2", "Dlgap1", "Frmd3", "Grin2a", "Etv6", "L3mbtl4", "Map9", "Nlgn1", "Nrg1", "Sh3bgrl2", "Slc12a7", "Rasgrf2", "Slc8a1", "Unc5d", "Wwtr1", "Neurod1", "Rgs20", "Notch1", "Notch2", "Nrxn1", "Sv2b", "Slc24a2", "Robo1", "Robo2", "Plxna4", "Slc1a3", "Apoe", "Gpr85", "Id2")

ht <- SCP::GroupHeatmap(seurat_object, features = hic_rna_features$gene, feature_split = hic_rna_features$Clusters, slot = "counts", assay = "SCT", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c("#0099CC", "white", "#CC0033"), nlabel = 0, show_row_names = F, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE), features_label = key_markers)
ht$plot

overlap_rna_features <- intersect(hic_rna_features$gene, GZ_CP$gene)

overlap_hic_features <- intersect(features, gene_infor_250k$gene_region[which(gene_infor_250k$gene_name %in% overlap_rna_features)])
GZ_CP_HiC_overlap <- GZ_CP_HiC[which(GZ_CP_HiC$features %in% overlap_hic_features),]

ht <- SCP::GroupHeatmap(tmp_object, features = GZ_CP_HiC_overlap$features, feature_split = GZ_CP_HiC_overlap$Clusters, assay = "scAB250kb", slot = "data", exp_method = "raw", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(tmp_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c(rev(col_scale_mako)[4:10], col_scale_acton[1:7]), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot

meta_data <- gene_infor_250k[which(gene_infor_250k$gene_region %in% overlap_hic_features & gene_infor_250k$gene_name %in% overlap_rna_features & gene_infor_250k$gene_type == "protein_coding"),]
meta_data <- unique(meta_data[,-c(1, which(colnames(meta_data) == "gene_id"))])
meta_data$region_gene <- paste0(meta_data$gene_region, "-", meta_data$gene_name)
meta_data$idx <- paste0(meta_data$gene_name, "-", 1:nrow(meta_data))
meta_data$Clusters <- GZ_CP_HiC$Clusters[match(meta_data$gene_region, GZ_CP_HiC$features)]

key_markers <- c("Adamts3", "Adamtsl1", "Ccser1", "Cdh20", "Cdh9", "Dach1", "Dlgap1", "Epb41l2", "Egfr", "Gpr85", "Grin3a", "Kcnj3", "Frmd3", "Meis1", "Meis2", "Neto1", "Nlgn1", "Notch2", "Nrg1", "Pax6", "Pcdh15", "Pcdh9", "Pou3f2", "Ptprd", "Ptprk", "Rad21", "Slc12a7", "Slc2a3", "Unc5d", "Yap1", "Zfhx4")

meta_data$label <- NA
meta_data$label[match(key_markers, meta_data$gene_name)] <- "T"

overlap_scAB_matrix_center <- scAB_matrix_center[meta_data$gene_region,]
rownames(overlap_scAB_matrix_center) <- meta_data$idx
colnames(overlap_scAB_matrix_center) <- c("E14_5_GZ", "E14_5_CP", "E16_5_CP", "E16_5_GZ", "E18_5_CP", "E18_5_GZ")
tmp_object[["overlap"]] <- CreateAssayObject(overlap_scAB_matrix_center)

rna_matrix <- GetAssayData(seurat_object, assay = "SCT", slot = "counts")
rna_matrix <- rna_matrix[meta_data$gene_name,]
rownames(rna_matrix) <- meta_data$idx
seurat_object[["overlap"]] <- CreateAssayObject(rna_matrix)

ht1 <- SCP::GroupHeatmap(tmp_object, features = meta_data$idx, feature_split = meta_data$Clusters, assay = "overlap", slot = "data", exp_method = "raw", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(hic_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c(rev(col_scale_mako)[4:10], col_scale_acton[1:7]), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE), features_label = meta_data$idx[which(meta_data$label == "T")])
ht1$plot

pdf("HiC_RNA_overlap_Clusters_hiC.pdf", width = 8, height = 12)
ht1$plot
dev.off()

ht2 <- SCP::GroupHeatmap(seurat_object, features = ht1$feature_metadata$features[order(ht1$feature_metadata$index)], feature_split = ht1$feature_metadata$feature_split[order(ht1$feature_metadata$index)], assay = "overlap", slot = "counts", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(hic_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = T, cluster_columns = FALSE, heatmap_palcolor = c("#0099CC", "white", "#CC0033"), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE), features_label = meta_data$idx[which(meta_data$label == "T")])
ht2$plot

pdf("HiC_RNA_overlap_Clusters_rna.pdf", width = 8, height = 12)
ht2$plot
dev.off()
#!---------------------------------------------------------------------------------------------------------------------------------------------
GZ_CP <- read.csv("GZ_CP_meta_data.csv")
gene_info <- gene_infor_100k[which(gene_infor_100k$gene_name %in% GZ_CP$gene & gene_infor_100k$is_promoter == "promoter"),]
GZ_CP <- GZ_CP[match(gene_info$gene_name, GZ_CP$gene),]
gene_info <- cbind(gene_info, GZ_CP)

#scAB_matrix <- GetAssayData(hic_object, assay = "scAB250kb", slot = "counts")
#scAB_matrix <- t(scAB_matrix)
#rownames(scAB_matrix) <- paste0(hic_object$time_Region)
#scAB_matrix <- t(limma::avereps(scAB_matrix))
#scAB_matrix <- scAB_matrix[which(rownames(scAB_matrix) %in% gene_info$gene_region),c("E14_5_GZ", "E16_5_GZ", "E18_5_GZ", "E14_5_CP", "E16_5_CP", "E18_5_CP")]
#gene_info <- gene_info[which(gene_info$gene_region %in% rownames(scAB_matrix)),]
#scAB_matrix <- scAB_matrix[gene_info$gene_region,]

scAB_matrix <- c()
for(i in c("E14_GZ", "E16_GZ", "E18_GZ", "E14_CP", "E16_CP", "E18_CP")){
	bulkAB_files <- read.table(paste0("Track/", i, ".bedGraph"), sep = "\t", header = F)
	bulkAB_files$gene_region <- paste0(bulkAB_files[,1], ":", bulkAB_files[,2], "-", bulkAB_files[,3])
	bulkAB_files <- bulkAB_files[,5:4]
	colnames(bulkAB_files) <- c("gene_region", paste0(i, "_hic"))
	if(i == "E14_GZ"){
		scAB_matrix <- bulkAB_files
	}else{
		scAB_matrix <- cbind(scAB_matrix, bulkAB_files)
	}
}
scAB_matrix <- scAB_matrix[,c(1, 2, 4, 6, 8, 10, 12)]
rownames(scAB_matrix) <- scAB_matrix$gene_region
scAB_matrix <- scAB_matrix[,2:7]
tmp <- intersect(gene_info$gene_region, rownames(scAB_matrix))

gene_info <- gene_info[which(gene_info$gene_region %in% tmp),]
scAB_matrix <- scAB_matrix[gene_info$gene_region,]

temp_object <- AverageExpression(seurat_object, assays = "SCT", features = unique(gene_info$gene), group.by = "time_Region", slot = "counts", return.seurat = T)
rna_matrix <- GetAssayData(temp_object, slot = "data")
rna_matrix <- rna_matrix[gene_info$gene,]
colnames(rna_matrix) <- paste0(c("E14_GZ", "E16_GZ", "E18_GZ", "E14_CP", "E16_CP", "E18_CP"), "_rna")

gene_info <- cbind(gene_info, scAB_matrix, rna_matrix)
gene_info <- gene_info[,c(1:13, 20:31)]

gene_info$cor_pearson <- NA
gene_info$p_value_pearson <- NA
gene_info$cor_spearman <- NA
gene_info$p_value_spearman <- NA

gene_info <- gene_info[order(gene_info$index),]
for(i in 1:nrow(gene_info)){
	corr <- cor.test(x = as.numeric(gene_info[i,14:19]), y = as.numeric(gene_info[i,20:25]), method = "pearson")
	gene_info$cor_pearson[i] <- corr[[4]]
	gene_info$p_value_pearson[i] <- corr[[3]]
	corr <- cor.test(x = as.numeric(gene_info[i,14:19]), y = as.numeric(gene_info[i,20:25]), method = "spearman")
	gene_info$cor_spearman[i] <- corr[[4]]
	gene_info$p_value_spearman[i] <- corr[[3]]
}
write.csv(gene_info, "GZ_CP_meta_data_new_100k.csv", quote = F, row.names = F)

gene_info <- gene_info[order(gene_info$cor_pearson, decreasing = T),]
gene_info$index_cor <- 1:nrow(gene_info)
key_markers <- c("Foxp1", "Ptprd", "Dlgap1", "Pcdh15", "Pax6", "Tafa2", "Meis2", "Sv2b", "Satb2", "Grin2a", "Lhx2", "Tafa1", "Notch2", "Prdm5", "Neurod1")
gene_info$is_label <- NA
gene_info$is_label[match(key_markers, gene_info$gene)] <- "label"

tmp_info <- gene_info[which(gene_info$is_label == "label"),]
gene_info <- gene_info[-which(gene_info$is_label == "label"),]
gene_info <- rbind(gene_info, tmp_info)

gene_info <- gene_info[which(gene_info$cor_pearson > 0.7),]

gene_info$Clusters <- factor(gene_info$Clusters, levels = paste0("C", 1:13))

key_markers <- c("Bcl11b", "Celf4", "Cux1", "Cux2", "Egfr", "Emx1", "Eomes", "Fabp7", "Fezf2", "Foxp1", "Lhx2", "Mef2c", "Meis2", "Nr2e1", "Nrp1", "Ntrk2", "Pax6", "Pex5l", "Pou3f2", "Pou3f3", "Ptprd", "Ptprk", "Slc1a3", "Sox2", "Sox5", "Tafa1", "Tafa2", "Tbr1", "Tle4", "Zbtb20", "Satb2", "Cadm2", "Dlgap1", "Frmd3", "Grin2a", "Etv6", "L3mbtl4", "Map9", "Nlgn1", "Nrg1", "Sh3bgrl2", "Slc12a7", "Rasgrf2", "Slc8a1", "Unc5d", "Wwtr1", "Neurod1", "Rgs20", "Notch1", "Notch2", "Nrxn1", "Sv2b", "Slc24a2", "Robo1", "Robo2", "Plxna4", "Slc1a3", "Apoe", "Gpr85", "Id2")

ht <- SCP::GroupHeatmap(seurat_object, features = gene_info$gene, feature_split = gene_info$Clusters, slot = "counts", assay = "SCT", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c("#0099CC", "white", "#CC0033"), nlabel = 0, show_row_names = F, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE), features_label = key_markers)
pdf("GZ_CP_RNA_Clusters_cor.pdf", width = 8, height = 12)
ht$plot
dev.off()
#!---------------------------------------
# useless
gene_info <- gene_info[which(gene_info$),]

gene_info <- gene_info[-which(gene_info$gene_region %in% c("chr1:172200000-172300000", "chr10:42500000-42600000", "chr14:73300000-73400000", "chr9:121200000-121300000")),]

scAB_matrix <- as.matrix(gene_info[, 14:19])
colnames(scAB_matrix) <- c("E14_GZ", "E16_GZ", "E18_GZ", "E14_CP", "E16_CP", "E18_CP")
rownames(scAB_matrix) <- gene_info$gene_region

mean_E14 <- rowMeans(scAB_matrix[,c(1, 4)])
mean_E16 <- rowMeans(scAB_matrix[,c(2, 5)])
mean_E18 <- rowMeans(scAB_matrix[,c(3, 6)])
scAB_matrix_center <- cbind(scAB_matrix[,c(1, 4)] - mean_E14, scAB_matrix[,c(2, 5)] - mean_E16, scAB_matrix[,c(3, 6)] - mean_E18)
colnames(scAB_matrix_center) <- c("E14_GZ", "E16_GZ", "E18_GZ", "E14_CP", "E16_CP", "E18_CP")
tmp_object <- CreateSeuratObject(counts = scAB_matrix, data = scAB_matrix_center)
tmp_object[["AB100kb"]] <- CreateAssayObject(data = scAB_matrix_center)
tmp_object$Region <- c("GZ", "GZ", "GZ", "CP", "CP", "CP")
tmp_object$orig.ident <- c("E14_5", "E16_5", "E18_5", "E14_5", "E16_5", "E18_5")
tmp_object$Region <- factor(tmp_object$Region, levels = c("GZ", "CP"))
tmp_object$orig.ident <- factor(tmp_object$orig.ident, levels = c("E14_5", "E16_5", "E18_5"))

ht <- SCP::GroupHeatmap(tmp_object, features = gene_info$gene_region, feature_split = gene_info$Clusters, assay = "AB100kb", slot = "data", exp_method = "raw", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(tmp_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c(rev(col_scale_mako)[4:10], col_scale_acton[1:7]), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot

pdf("GZ_CP_HiC_Clusters1.pdf", width = 8, height = 12)
ht$plot
dev.off()
#!---------------------------------------
tmp_rna <- normalize_to_01(as.vector(as.matrix(gene_info[,20:25])))
tmp_rna <- matrix(tmp_rna, nrow = nrow(gene_info))
colnames(tmp_rna) <- colnames(gene_info)[20:25]
gene_info[,20:25] <- tmp_rna

tmp_hic <- normalize_to_01(as.vector(as.matrix(gene_info[,14:19])))
tmp_hic <- matrix(tmp_hic, nrow = nrow(gene_info))
colnames(tmp_hic) <- colnames(gene_info)[14:19]
gene_info[,14:19] <- tmp_hic

tmp_info <- data.frame(gene = rep(gene_info$gene, each = 12), Cluster = rep(gene_info$Clusters, each = 12), method = rep(c(rep(c("HiC", "RNA"), each = 6)), nrow(gene_info)), time_Region = rep(colnames(gene_info)[14:25], nrow(gene_info)), values = as.vector(unlist(t(gene_info[,14:25]))))
tmp_info$time_Region <- gsub(tmp_info$time_Region, pattern = "_hic|_rna", replacement = "")
tmp_info$time_Region <- factor(tmp_info$time_Region, levels = c("E14_GZ", "E16_GZ", "E18_GZ", "E14_CP", "E16_CP", "E18_CP"))

tmp_info$gene_method <- paste0(tmp_info$gene, "_", tmp_info$method)

pdf("correlation_RNA_HiC_GZ_vs_CP_per_gene.pdf", width = 12)
for(i in unique(tmp_info$gene)){
	p0 <- ggplot(tmp_info[which(tmp_info$gene == i),], aes(x = time_Region, y = values, group = gene_method)) + 
			geom_line(aes(color = method), size = 1, show.legend = T) + geom_point(aes(color = method), shape = 19, size = 3, show.legend = F) + 
			scale_color_manual(values = c("darkred", "darkgreen")) + 
			theme_classic(base_size = 10) + xlab("time Region") + ylab("normalized expression/mean scAB values") + ylim(c(0, 1)) + labs(title = i) +
			theme(panel.grid = element_blank(), legend.position = 'right', legend.justification = c(0, 1), plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_text(size = 12, color = "black", face = "bold")) + RotatedAxis()
	print(p0)
}
dev.off()

dir.create("Cortex_GZ_CP_Region_correlation_features")
for(i in unique(gene_info$Clusters)){
        gene_list <- as.character(gene_info[which(gene_info$Clusters == i),'gene'])
        enrichRes <- enrichAnalysis(genelist = gene_list, geneType = "SYMBOL", species = "mmu", database = c("go", "kegg", "MSigDb"), GO.model = "BP", MSigDb.signature = "H", sampleName = paste0("enrich_", i), minGSSize = 5, outpath = "Cortex_GZ_CP_Region_correlation_features")
}
cell_cluster_colors <- c(C1 = "#C22B86", C2 = "#F2639E", C3 = "#FBB5BB", C4 = "#5A5A89", C5 = "#1F73A7", C6 = "#68B2D3", C7 = "#4F928A", C8 = "#1C7B44", C9 = "#74C174", C10 = "#DCCF79", C11 = "#E77306", C12 = "#F28F2C", C13 = "#FFBB61")
for(i in c("C1", "C3", "C4", "C6", "C11", "C12", "C13")){
	enrich_go_terms <- read.csv(paste0("Cortex_GZ_CP_Region_correlation_features/enrich_", i, "_enrich_GeneOntology_BP_terms_selected.csv"))
	enrich_go_terms$logFDR <- -log10(enrich_go_terms$p.adjust)
	enrich_go_terms$Description <- factor(enrich_go_terms$Description, levels = rev(enrich_go_terms$Description))
	enrich_go_terms$Clusters <- i
	p0 <- ggplot(enrich_go_terms, aes_string(x = "logFDR", y = "Description", fill = "Clusters")) + theme_classic() + geom_col() + scale_fill_manual(values = cell_cluster_colors) + NoLegend() +
		geom_text(data = enrich_go_terms, aes(x = 0.05, y = Description, label = Description), size = 4, hjust = 0) +
	  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line = element_line(color = 'black', size = 1.1), axis.text = element_text(size = 12, color = "black"), plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_text(size = 16, color = "black", face = "bold")) + xlab("-log10(p.adjust)") + labs(title = i)
     print(p0)
     ggsave(paste0("Cortex_GZ_CP_Region_correlation_features/enrich_", i, "_enrich_GeneOntology_BP_terms_selected.pdf"), width = 6, height = 3.5)
}

library(ggforce)
library(ggwordcloud)

gene_info <- gene_info[which(gene_info$p_value_pearson < 0.05),]
pdf("correlation_RNA_HiC_GZ_vs_CP_wordcloud.pdf")
for(i in levels(gene_info$Clusters)){
	try({
		tmp_info <- gene_info[which(gene_info$Clusters == i),]
		tmp_info$angle <- 90 * sample(c(0, 1), nrow(tmp_info), replace = TRUE, prob = c(60, 40))
		tmp_info$p_value_pearson <- max(tmp_info$p_value_pearson) - tmp_info$p_value_pearson
		colors <- SCP::palette_scp(tmp_info$p_value_pearson, type = "continuous", palette = "Spectral", palcolor = NULL, matched = FALSE)
		colors_value <- seq(min(tmp_info$p_value_pearson, na.rm = TRUE), quantile(tmp_info$p_value_pearson, 0.99, na.rm = TRUE) + 0.001, length.out = 100)
			  
		p0 <- ggplot(tmp_info, aes(label = gene, size = cor_pearson, color = p_value_pearson, angle = angle)) + 
				geom_text_wordcloud(rm_outside = TRUE, eccentricity = 1, shape = "square", show.legend = TRUE, grid_margin = 3) + # star, cardioid, diamond, square, triangle-forward, triangle-upright, pentagon, circle
				scale_color_gradientn(name = "P values:", colours = colors, values = scales::rescale(colors_value), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barheight = 4, barwidth = 1, title.hjust = 0)) +
				scale_size(name = "Count", range = c(2, 8), breaks = ceiling(seq(min(tmp_info$cor_pearson * 10, na.rm = TRUE), max(tmp_info$cor_pearson * 10, na.rm = TRUE), length.out = 3))) +
				guides(size = guide_legend(override.aes = list(colour = "black", label = "G"), order = 1)) + labs(title = i) + SCP::theme_scp() +
				theme(aspect.ratio = 1, legend.box.margin = margin(0, 0, 0, 0), legend.margin = margin(0, 0, 0, 0), legend.position = "right", legend.direction = "vertical", plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_text(size = 16, color = "black", face = "bold"))
		print(p0)
	})
}
dev.off()
#!------------------------------------------------------------------------------------------------------------------------
gene_info <- read.csv("GZ_CP_meta_data_new_100k.csv")
Idents(seurat_object) <- "time_Region"
rna_FC_E14 <- FoldChange(seurat_object, features = gene_info$gene, ident.1 = "E14_5_CP", ident.2 = "E14_5_GZ")
rna_FC_E16 <- FoldChange(seurat_object, features = gene_info$gene, ident.1 = "E16_5_CP", ident.2 = "E16_5_GZ")
rna_FC_E18 <- FoldChange(seurat_object, features = gene_info$gene, ident.1 = "E18_5_CP", ident.2 = "E18_5_GZ")
colnames(rna_FC_E14) <- paste0("E14_rna_", colnames(rna_FC_E14))
colnames(rna_FC_E16) <- paste0("E16_rna_", colnames(rna_FC_E16))
colnames(rna_FC_E18) <- paste0("E18_rna_", colnames(rna_FC_E18))

Idents(seurat_object) <- "Region"
rna_FC <- FoldChange(seurat_object, features = gene_info$gene, ident.1 = "CP", ident.2 = "GZ")

gene_info <- cbind(gene_info, rna_FC_E14$E14_rna_avg_log2FC, rna_FC_E16$E16_rna_avg_log2FC, rna_FC_E18$E18_rna_avg_log2FC)
colnames(gene_info)[30:32] <- c("E14_avg_log2FC", "E16_avg_log2FC", "E18_avg_log2FC")

gene_info$E14_AB_change <- gene_info$E14_CP_hic - gene_info$E14_GZ_hic
gene_info$E16_AB_change <- gene_info$E16_CP_hic - gene_info$E16_GZ_hic
gene_info$E18_AB_change <- gene_info$E18_CP_hic - gene_info$E18_GZ_hic

gene_info$rna_FC <- rna_FC$avg_log2FC
gene_info$hic_AB_change <- rowMeans(gene_info[,c("E14_AB_change", "E16_AB_change", "E18_AB_change")])

write.csv(gene_info, "GZ_CP_meta_data_new_100k.csv", row.names = F, quote = F)

key_markers <- c("Pax6", "Tafa1", "Tafa2", "Dlgap1", "Grin2a", "Sv2b", "Notch2", "Zbtb20")
gene_info$is_label <- NA
gene_info$is_label[match(key_markers, gene_info$gene)] <- "label"
gene_info <- gene_info[order(abs(gene_info$cor_pearson)),]

tmp_info <- gene_info[which(gene_info$is_label == "label"),]
gene_info <- gene_info[-which(gene_info$is_label == "label"),]
gene_info <- rbind(gene_info, tmp_info)

p <- ggplot() +
		geom_point(data = gene_info, aes(x = rna_FC, y = hic_AB_change, color = cor_pearson, size = cor_pearson), alpha = 0.8) +
		scale_size(name = "cor", range = c(0.5, 2), breaks = ceiling(seq(min(abs(gene_info$cor_pearson) * 10, na.rm = TRUE), max(gene_info$cor_pearson * 10, na.rm = TRUE), length.out = 3))) +
		geom_hline(yintercept = 0, color = "black", linetype = 2) + 
		geom_vline(xintercept = 0, color = "grey", linetype = 2) +
		geom_text_repel(data = gene_info[which(gene_info$is_label == "label"),], 
						aes(x = rna_FC, y = hic_AB_change, label = gene), min.segment.length = 0, max.overlaps = 100, segment.colour = "grey40", color = "black", bg.color = "white", bg.r = 0.1, size = 4, force = 20) +
		labs(x = "Spatial RNA avg logFC", y = "Spatial HiC scAB change values") +
		scale_color_gradientn(
        name = "pearson cor", colors = SCP::palette_scp(gene_info$cor_pearson, type = "continuous", palette = "Spectral", palcolor = NULL, matched = FALSE),
        values = scales::rescale(seq(min(gene_info$cor_pearson, na.rm = TRUE), quantile(gene_info$cor_pearson, 0.99, na.rm = TRUE) + 0.001, length.out = 100)),
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barheight = 4, barwidth = 1, title.hjust = 0, order = 1)) +
		theme(aspect.ratio = NULL)
print(p)
ggsave("GZ_CP_FC_pointplot_all_FC_time_new.pdf")

p1 <- ggplot() +
		geom_point(data = gene_info, aes(x = E14_avg_log2FC, y = E14_AB_change, color = cor_pearson, size = cor_pearson), alpha = 0.8) +
		scale_size(name = "cor", range = c(0.5, 2), breaks = ceiling(seq(min(abs(gene_info$cor_pearson) * 10, na.rm = TRUE), max(gene_info$cor_pearson * 10, na.rm = TRUE), length.out = 3))) +
		geom_hline(yintercept = 0, color = "black", linetype = 2) +
		geom_vline(xintercept = 0, color = "grey", linetype = 2) +
		geom_text_repel(data = gene_info[which(gene_info$is_label == "label"),], 
						aes(x = E14_avg_log2FC, y = E14_AB_change, label = gene), min.segment.length = 0, max.overlaps = 100, segment.colour = "grey40", color = "black", bg.color = "white", bg.r = 0.1, size = 4, force = 20) +
		labs(x = "Spatial RNA avg logFC", y = "Spatial HiC scAB change values", title = "E14 GZ vs CP") +
		scale_color_gradientn(
        name = "pearson cor", colors = SCP::palette_scp(gene_info$cor_pearson, type = "continuous", palette = "Spectral", palcolor = NULL, matched = FALSE),
        values = scales::rescale(seq(min(gene_info$cor_pearson, na.rm = TRUE), quantile(gene_info$cor_pearson, 0.99, na.rm = TRUE) + 0.001, length.out = 100)),
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barheight = 4, barwidth = 1, title.hjust = 0, order = 1)) +
		theme(aspect.ratio = NULL)
p2 <- ggplot() +
		geom_point(data = gene_info, aes(x = E16_avg_log2FC, y = E16_AB_change, color = cor_pearson, size = cor_pearson), alpha = 0.8) +
		scale_size(name = "cor", range = c(0.5, 2), breaks = ceiling(seq(min(abs(gene_info$cor_pearson) * 10, na.rm = TRUE), max(gene_info$cor_pearson * 10, na.rm = TRUE), length.out = 3))) +
		geom_hline(yintercept = 0, color = "black", linetype = 2) +
		geom_vline(xintercept = 0, color = "grey", linetype = 2) +
		geom_text_repel(data = gene_info[which(gene_info$is_label == "label"),], 
						aes(x = E16_avg_log2FC, y = E16_AB_change, label = gene), min.segment.length = 0, max.overlaps = 100, segment.colour = "grey40", color = "black", bg.color = "white", bg.r = 0.1, size = 4, force = 20) +
		labs(x = "Spatial RNA avg logFC", y = "Spatial HiC scAB change values", title = "E16 GZ vs CP") +
		scale_color_gradientn(
        name = "pearson cor", colors = SCP::palette_scp(gene_info$cor_pearson, type = "continuous", palette = "Spectral", palcolor = NULL, matched = FALSE),
        values = scales::rescale(seq(min(gene_info$cor_pearson, na.rm = TRUE), quantile(gene_info$cor_pearson, 0.99, na.rm = TRUE) + 0.001, length.out = 100)),
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barheight = 4, barwidth = 1, title.hjust = 0, order = 1)) +
		theme(aspect.ratio = NULL)
p3 <- ggplot() +
		geom_point(data = gene_info, aes(x = E18_avg_log2FC, y = E18_AB_change, color = cor_pearson, size = cor_pearson), alpha = 0.8) +
		scale_size(name = "cor", range = c(0.5, 2), breaks = ceiling(seq(min(abs(gene_info$cor_pearson) * 10, na.rm = TRUE), max(gene_info$cor_pearson * 10, na.rm = TRUE), length.out = 3))) +
		geom_hline(yintercept = 0, color = "black", linetype = 2) +
		geom_vline(xintercept = 0, color = "grey", linetype = 2) +
		geom_text_repel(data = gene_info[which(gene_info$is_label == "label"),], 
						aes(x = E18_avg_log2FC, y = E18_AB_change, label = gene), min.segment.length = 0, max.overlaps = 100, segment.colour = "grey40", color = "black", bg.color = "white", bg.r = 0.1, size = 4, force = 20) +
		labs(x = "Spatial RNA avg logFC", y = "Spatial HiC scAB change values", title = "E18 GZ vs CP") +
		scale_color_gradientn(
        name = "pearson cor", colors = SCP::palette_scp(gene_info$cor_pearson, type = "continuous", palette = "Spectral", palcolor = NULL, matched = FALSE),
        values = scales::rescale(seq(min(gene_info$cor_pearson, na.rm = TRUE), quantile(gene_info$cor_pearson, 0.99, na.rm = TRUE) + 0.001, length.out = 100)),
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barheight = 4, barwidth = 1, title.hjust = 0, order = 1)) +
		theme(aspect.ratio = NULL)
print(p1 + p2 + p3)
ggsave("GZ_CP_FC_pointplot_per_FC_time.pdf", width = 21, height = 6.5)

p4 <- ggplot() +
		geom_point(data = gene_info[which(is.na(gene_info$is_label)),], aes(x = E14_avg_log2FC, y = E14_AB_change, color = cor_pearson, size = cor_pearson), alpha = 0.6, shape = 17) +
		geom_point(data = gene_info[which(is.na(gene_info$is_label)),], aes(x = E16_avg_log2FC, y = E16_AB_change, color = cor_pearson, size = cor_pearson), alpha = 0.6, shape = 18) +
		geom_point(data = gene_info[which(is.na(gene_info$is_label)),], aes(x = E18_avg_log2FC, y = E18_AB_change, color = cor_pearson, size = cor_pearson), alpha = 0.6, shape = 16) +
		geom_point(data = gene_info[which(gene_info$is_label == "label"),], aes(x = E14_avg_log2FC, y = E14_AB_change, color = cor_pearson, size = cor_pearson), alpha = 1, shape = 17) +
		geom_point(data = gene_info[which(gene_info$is_label == "label"),], aes(x = E16_avg_log2FC, y = E16_AB_change, color = cor_pearson, size = cor_pearson), alpha = 1, shape = 18) +
		geom_point(data = gene_info[which(gene_info$is_label == "label"),], aes(x = E18_avg_log2FC, y = E18_AB_change, color = cor_pearson, size = cor_pearson), alpha = 1, shape = 16) +
		scale_size(name = "cor", range = c(0.5, 2), breaks = ceiling(seq(min(abs(gene_info$cor_pearson) * 10, na.rm = TRUE), max(gene_info$cor_pearson * 10, na.rm = TRUE), length.out = 3))) +
		geom_hline(yintercept = 0, color = "black", linetype = 2) +
		geom_vline(xintercept = 0, color = "grey", linetype = 2) +
		geom_text_repel(data = gene_info[which(gene_info$is_label == "label"),], 
						aes(x = E14_avg_log2FC, y = E14_AB_change, label = gene), min.segment.length = 0, max.overlaps = 100, segment.colour = "grey40", color = "black", bg.color = "white", bg.r = 0.1, size = 4, force = 20) +
		geom_text_repel(data = gene_info[which(gene_info$is_label == "label"),], 
						aes(x = E16_avg_log2FC, y = E16_AB_change, label = gene), min.segment.length = 0, max.overlaps = 100, segment.colour = "grey40", color = "black", bg.color = "white", bg.r = 0.1, size = 4, force = 20) +
		geom_text_repel(data = gene_info[which(gene_info$is_label == "label"),], 
						aes(x = E18_avg_log2FC, y = E18_AB_change, label = gene), min.segment.length = 0, max.overlaps = 100, segment.colour = "grey40", color = "black", bg.color = "white", bg.r = 0.1, size = 4, force = 20) +
		labs(x = "Spatial RNA avg logFC", y = "Spatial HiC scAB change values", title = "GZ vs CP cross cortex development") +
		scale_color_gradientn(
        name = "pearson cor", colors = SCP::palette_scp(gene_info$cor_pearson, type = "continuous", palette = "Spectral", palcolor = NULL, matched = FALSE),
        values = scales::rescale(seq(min(gene_info$cor_pearson, na.rm = TRUE), quantile(gene_info$cor_pearson, 0.99, na.rm = TRUE) + 0.001, length.out = 100)),
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barheight = 4, barwidth = 1, title.hjust = 0, order = 1)) +
		theme(aspect.ratio = NULL, plot.title = element_text(hjust = 0.5, size = 18), axis.title = element_text(size = 16, color = "black", face = "bold"))
print(p4)
ggsave("GZ_CP_FC_pointplot_per_FC_time_combined.pdf", width = 8, height = 8)
gene_info <- gene_info[which(gene_info$cor_pearson > 0.7),]

gene_info$Clusters <- factor(gene_info$Clusters, levels = paste0("C", 1:13))

key_markers <- c("Bcl11b", "Celf4", "Cux1", "Cux2", "Egfr", "Emx1", "Eomes", "Fabp7", "Fezf2", "Foxp1", "Lhx2", "Mef2c", "Meis2", "Nr2e1", "Nrp1", "Ntrk2", "Pax6", "Pex5l", "Pou3f2", "Pou3f3", "Ptprd", "Ptprk", "Slc1a3", "Sox2", "Sox5", "Tafa1", "Tafa2", "Tbr1", "Tle4", "Zbtb20", "Satb2", "Cadm2", "Dlgap1", "Frmd3", "Grin2a", "Etv6", "L3mbtl4", "Map9", "Nlgn1", "Nrg1", "Sh3bgrl2", "Slc12a7", "Rasgrf2", "Slc8a1", "Unc5d", "Wwtr1", "Neurod1", "Rgs20", "Notch1", "Notch2", "Nrxn1", "Sv2b", "Slc24a2", "Robo1", "Robo2", "Plxna4", "Slc1a3", "Apoe", "Gpr85", "Id2")

ht <- SCP::GroupHeatmap(seurat_object, features = gene_info$gene, feature_split = gene_info$Clusters, slot = "counts", assay = "SCT", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c("#0099CC", "white", "#CC0033"), nlabel = 0, show_row_names = F, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE), features_label = key_markers)
pdf("GZ_CP_RNA_Clusters_cor_0.7.pdf", width = 8, height = 12)
ht$plot
dev.off()
#!---------------------------------------
# useless
gene_info <- gene_info[-which(gene_info$gene_region %in% c("chr1:172200000-172300000", "chr10:42500000-42600000", "chr14:73300000-73400000", "chr9:121200000-121300000")),]

scAB_matrix <- as.matrix(gene_info[, 14:19])
colnames(scAB_matrix) <- c("E14_GZ", "E16_GZ", "E18_GZ", "E14_CP", "E16_CP", "E18_CP")
rownames(scAB_matrix) <- gene_info$gene_region

tmp_object <- CreateSeuratObject(counts = scAB_matrix)
tmp_object[["AB100k"]] <- CreateAssayObject(data = scAB_matrix)
tmp_object$Region <- c("GZ", "GZ", "GZ", "CP", "CP", "CP")
tmp_object$orig.ident <- c("E14_5", "E16_5", "E18_5", "E14_5", "E16_5", "E18_5")
tmp_object$Region <- factor(tmp_object$Region, levels = c("GZ", "CP"))
tmp_object$orig.ident <- factor(tmp_object$orig.ident, levels = c("E14_5", "E16_5", "E18_5"))

gene_info <- gene_info[which(gene_info$Clusters != "C7"),]

ht <- SCP::GroupHeatmap(tmp_object, features = gene_info$gene_region, feature_split = gene_info$Clusters, assay = "AB100k", slot = "data", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(tmp_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c(rev(col_scale_mako)[4:10], col_scale_acton[1:7]), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot

pdf("GZ_CP_HiC_Clusters_cor_0.7.pdf", width = 8, height = 12)
ht$plot
dev.off()

gene_info$final_cluster1 <- as.character(gene_info$Clusters)
gene_info$final_cluster1[which(gene_info$final_cluster1 %in% paste0("C", 1:7))] <- "C1"
gene_info$final_cluster1[which(gene_info$final_cluster1 %in% paste0("C", 8:13))] <- "C2"
gene_info$final_cluster1 <- factor(gene_info$final_cluster1, levels = c("C1", "C2"))

scAB_matrix <- scAB_matrix[gene_info$gene_region,]
rownames(scAB_matrix) <- gene_info$gene
tmp_object[["overlap"]] <- CreateAssayObject(data = scAB_matrix)

ht2 <- SCP::GroupHeatmap(seurat_object, features = gene_info$gene, feature_split = gene_info$Clusters, assay = "SCT", slot = "counts", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(hic_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = T, cluster_columns = FALSE, heatmap_palcolor = c("#0099CC", "white", "#CC0033"), nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE), features_label = key_markers)
ht2$plot

pdf("GZ_CP_HiC_Clusters_cor_0.7_RNA_cluster.pdf", width = 8, height = 12)
ht2$plot
dev.off()

ht1 <- SCP::GroupHeatmap(tmp_object, features = ht2$feature_metadata$features[order(ht2$feature_metadata$index)], feature_split = ht2$feature_metadata$feature_split[order(ht2$feature_metadata$index)], slot = "data", assay = "overlap", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(tmp_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = F, cluster_columns = FALSE, heatmap_palcolor = c(rev(col_scale_mako)[4:10], col_scale_acton[1:7]), nlabel = 0, show_row_names = F, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE), features_label = key_markers)
ht1$plot

pdf("GZ_CP_HiC_Clusters_cor_0.7_RNA_order.pdf", width = 8, height = 12)
ht1$plot
dev.off()

ht1 <- SCP::GroupHeatmap(tmp_object, features = gene_info$gene, feature_split = gene_info$Clusters, slot = "data", assay = "overlap", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(tmp_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c(rev(col_scale_mako)[4:10], col_scale_acton[1:7]), nlabel = 0, show_row_names = F, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE), features_label = key_markers)
ht1$plot

pdf("GZ_CP_HiC_Clusters_cor_0.7.pdf", width = 8, height = 12)
ht1$plot
dev.off()
#!------------------------------------------------------------------------------------------------------------------------
DefaultAssay(seurat_object) <- "spARC_SCT"
DefaultAssay(hic_object) <- "spARC_scAB250kb_scale"
seurat_object$shape = "RNA"
hic_object$shape = "HiC"

i = "Tafa1"
i = "Nr4a2"

p0 <- Seurat::SpatialPlot(seurat_object, features = i, shape = 22, crop = F) & scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")))
p1 <- Seurat::SpatialPlot(hic_object, shape = 22, crop = F, features = unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name == i & gene_infor_250k$is_promoter == "promoter")])) & scale_fill_gradientn(values = seq(0, 1, 0.01), colours = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"))
print(p0/p1)

p1 <- Seurat::SpatialPlot(hic_object, shape = 22, crop = F, features = unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name == i & gene_infor_250k$is_promoter == "promoter")]), images = "E14_5") & scale_fill_gradientn(colours = paletteContinuous(set = "solarExtra", n = 240)[1:80])
p2 <- Seurat::SpatialPlot(hic_object, shape = 22, crop = F, features = unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name == i & gene_infor_250k$is_promoter == "promoter")]), images = "E16_5") & scale_fill_gradientn(colours = paletteContinuous(set = "solarExtra", n = 240)[80:150])
p3 <- Seurat::SpatialPlot(hic_object, shape = 22, crop = F, features = unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name == i & gene_infor_250k$is_promoter == "promoter")]), images = "E18_5") & scale_fill_gradientn(colours = paletteContinuous(set = "solarExtra", n = 240)[90:240])
print(p1 + p2 + p3)
ggsave(paste0(i, "_Spatial.pdf"))
p0 <- FeaturePlot(seurat_object, reduction = "spatial_umap", features = i, shape.by = "shape", pt.size = 0.8) & scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) & scale_shape_manual(values = 15) & ylim(c(0, 1080)) & xlim(c(0, 5500))
p1 <- FeaturePlot(hic_object, features = unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name == i & gene_infor_250k$is_promoter == "promoter")]), reduction = "spatial_umap", shape.by = "shape", pt.size = 0.8) & scale_colour_gradientn(values = seq(0, 1, 0.01), colours = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#F8BB3A", "#FCB31C", "#F99F1C","#FF3606","#E30407")) & scale_shape_manual(values = 15) & ylim(c(0, 1080)) & xlim(c(0, 5500))
print(p0/p1)
ggsave(paste0(i, "_Spatial0.pdf"))

DefaultAssay(seurat_object) <- "SCT"
DefaultAssay(hic_object) <- "scAB250kb_scale"

p0 <- VlnPlot(seurat_object, features = i, pt.size = 0, group.by = "orig.ident", split.by = "Region", cols = c("#66C2A5", "#FC8D62"))
p1 <- VlnPlot(hic_object, pt.size = 0, features = unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name == i & gene_infor_250k$is_promoter == "promoter")]), group.by = "orig.ident", split.by = "Region", cols = c("#66C2A5", "#FC8D62")) + ylab("scAB values") & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.3, position = position_dodge(0.9), show.legend = F)
print(p0 + p1)
ggsave(paste0(i, "_vlnplot.pdf"))

P_values <- c()
for(j in levels(hic_object$orig.ident)){
	GZ <- GetAssayData(hic_object, assay = "scAB250kb_scale", layer = "data")[unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name == i & gene_infor_250k$is_promoter == "promoter")]), which(hic_object$orig.ident == j & hic_object$Region == "GZ")]
	CP <- GetAssayData(hic_object, assay = "scAB250kb_scale", layer = "data")[unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name == i & gene_infor_250k$is_promoter == "promoter")]), which(hic_object$orig.ident == j & hic_object$Region == "CP")]
	P_valves_t <- t.test(GZ, CP)[3]
	P_valves_w <- wilcox.test(GZ, CP)[3]
	P_values <- rbind(P_values, c(i, j, P_valves_t, P_valves_w))
}
colnames(P_values) <- c("feature", "group", "t_test_P_value", "wilcox_test_P_value")
write.csv(P_values, paste0(i, "_HiC_mean_scAB_values_P_values.csv"), quote = F, row.names = F)

VlnPlot(seurat_object, features = i, pt.size = 0, group.by = "orig.ident", split.by = "Cluster", cols = ClusterIDColors)

VlnPlot(seurat_object, features = i, pt.size = 0, group.by = "Cluster", cols = ClusterIDColors) + NoLegend()
ggsave(paste0(i, "_SCT_Cluster_vlnplot.pdf"))
VlnPlot(hic_object, pt.size = 0, features = unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name == i & gene_infor_250k$is_promoter == "promoter")]), group.by = "banksy_cluster", cols = hic_object@misc$banksy_cluster_colors) + ylab("scAB values") + NoLegend() & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.3, show.legend = F)
ggsave(paste0(i, "_scAB_banksy_cluster_vlnplot.pdf"))

temp_object <- readRDS("output/development_brain_mouse_scRNA_2021_Nature.rds")
DimPlot(temp_object, reduction = "cca_umap", cols = selected_colors) + SetAxes()
ggsave("scRNA_2021_Nature_E14_16_18_cell_identity.pdf")
FeaturePlot(temp_object, reduction = "cca_umap", features = "Fam19a1", order = T, min.cutoff = "q10", max.cutoff = "q99", cols = c("lightgrey", "red")) + SetAxes()
ggsave("scRNA_2021_Nature_E14_16_18_tafa1_expression.pdf")
VlnPlot(temp_object, features = "Fam19a1", cols = selected_colors) + NoLegend()
ggsave("scRNA_2021_Nature_E14_16_18_tafa1_expression_Vlnplot.pdf")
VlnPlot(temp_object, features = "Fam19a1", idents = "SCPN", split.by = "orig.ident")
ggsave("scRNA_2021_Nature_E14_16_18_tafa1_SCPN_expression_Vlnplot.pdf")
#!------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LGE CPU LS features
conserved_features <- read.csv("Region_all_conserved_markers_cell_type.csv")
de_features <- read.csv("Time_Region_ST.csv")

TF_file_list <- list.files("Trajectory_analysis/moscot/tutorials/figures_MSN_Spatial/")
TF_file_list <- TF_file_list[grep("*push.csv", TF_file_list)]
TF_file_list <- gsub(TF_file_list, pattern = ".csv", replacement = "")
core_list <- gsub(TF_file_list, pattern = "_push", replacement = "")
core_list <- gsub(core_list, pattern = "_drivers_tf", replacement = "")
core_genes <- list()
for(i in 1:length(TF_file_list)){
	tmp_data <- read.csv(paste0("Trajectory_analysis/moscot/tutorials/figures_MSN_Spatial/", TF_file_list[i], ".csv"))
	colnames(tmp_data) <- c("genes", "corr", "pval", "qval", "ci_low", "ci_high")
	tmp_data$FDR <- log10(tmp_data$qval)
	tmp_data <- tmp_data[-grep("[0-9]Rik$", tmp_data$genes),]
	core_genes[[core_list[i]]] <- tmp_data$genes[which(tmp_data$corr >= 0.4)]
}

Region_features <- unique(c(unlist(core_genes), conserved_features$gene[which(conserved_features$final_cluster %in% c("LGE", "CPU", "LS"))], de_features$gene[which(de_features$final_cluster %in% c("LGE", "CPU", "LS") & de_features$pct.2 < 0.3)]))

DefaultAssay(seurat_object) <- "spARC_SCT"
DefaultAssay(hic_object) <- "spARC_scAB250kb_scale"
dir.create("Region_features")
for(i in Region_features){
	try({
		p0 <- Seurat::SpatialPlot(seurat_object, features = i, shape = 22, crop = F) & scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")))
		p1 <- Seurat::SpatialPlot(hic_object, shape = 22, crop = F, features = unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name == i & gene_infor_250k$is_promoter == "promoter")])) & scale_fill_gradientn(values = seq(0, 1, 0.01), colours = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"))
		print(p0/p1)
		ggsave(paste0("Region_features/", i, ".png"), width = 24, height = 24)
		p0 <- FeaturePlot(seurat_object, reduction = "spatial_umap", features = i) & scale_colour_gradientn(colours = paletteContinuous(set = "solarExtra"))
		p1 <- FeaturePlot(hic_object, features = unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name == i & gene_infor_250k$is_promoter == "promoter")]), reduction = "spatial_umap") & scale_colour_gradientn(colours = paletteContinuous(set = "solarExtra"))
		print(p0/p1)
		ggsave(paste0("Region_features/", i, "_0.png"), width = 20, height = 10)
	})
}

DefaultAssay(seurat_object) <- "SCT"
DefaultAssay(hic_object) <- "scAB250kb_scale"
for(i in Region_features){
	try({
		p0 <- VlnPlot(seurat_object, features = i, pt.size = 0, group.by = "orig.ident", split.by = "Region", cols = cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))])
		p1 <- VlnPlot(hic_object, pt.size = 0, features = unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name == i & gene_infor_250k$is_promoter == "promoter")]), group.by = "orig.ident", split.by = "Region", cols = cluster_colors[match(levels(hic_object$Region), names(cluster_colors))]) + ylab("scAB values") & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.3, position = position_dodge(0.9), show.legend = F)
		print(p0 + p1)
		ggsave(paste0("Region_features/", i, "_vlnplot.png"), width = 16, height = 14)
	})
}
#!----------------------------------------------------------------------------
hic_object$Short_Long_Ratio[which(hic_object$Short_Long_Ratio < 1)] <- 1
hic_object$Short_Long_Ratio[which(hic_object$Short_Long_Ratio > 3)] <- 3

Seurat::SpatialPlot(hic_object, shape = 22, crop = F, features = "Short_Long_Ratio", max.cutoff = "q98") & scale_fill_gradientn(colours = paletteContinuous(set = "sambaNight"))
ggsave("Short_Long_Ratio.pdf")

Seurat::SpatialPlot(seurat_object, shape = 22, crop = F, features = "CytoTRACE2_Relative") & scale_fill_gradientn(colours = rev(paletteContinuous(set = "beach")))
ggsave("CytoTRACE2_score_all_spots.pdf")
seurat_object <- subset(seurat_object, subset = Region != "Others")
source("~/scripts/seurat2scanpy/shiny_st.R")
options(browser = "/usr/bin/firefox")
seurat_object <- shiny_st(seurat = seurat_object, isVisium = F, assay = "SCT", image = "E14_5")
seurat_object <- subset(seurat_object, subset = CytoTRACE2_Relative != "removed_cells")
Seurat::SpatialPlot(seurat_object, shape = 22, crop = F, features = "CytoTRACE2_Relative") & scale_fill_gradientn(colours = paletteContinuous(set = "beach"))
ggsave("CytoTRACE2_score.pdf")














