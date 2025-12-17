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

cluster_colors <- c("E14_5" = "#EF5EB4", "E16_5" = "#8A2BD2", "E18_5" = "#0047AB", "Adult15" = "#B312A6",
"VZ" = "#61E2A4", "Cortex" = "#59B375", "MGE" = "#F97B72", "LGE" = "#89288F", "CPU" = "#CC61B0", "ChP" = "#FF26A8", "LCS" = "#C1FF73", "CLA" = "#33FF00", "LS" = "#E68316", "Epd" = "#94FFB5", "PIR" = "#B5EFB5", "AEP" = "#653EB3", "others" = "#FEE52C",
"AP" = "#03FFF4", "Cortex_DP_AP" = "#03FFF4", "Cortex_MP_AP" = "#036DF4", "IPC" = "#0BD3B1", "Cortex_IPC" = "#0BD3B1", "MigN" = "#62CFE8", "Cortex_DP_MigN" = "#62CFE8", "Cortex_MP_MigN" = "#2F7DD1",
"Cortex_Layer_6b" = "#99FFBC", "CThPN" = "#7EC136", "Cortex_CThPN" = "#7EC136", "Cortex_DP_CThPN" = "#00CC14", "Cortex_MP_CThPN" = "#2D739E", "SCPN" = "#34A047", "Cortex_SCPN" = "#34A047", "Cortex_DP_SCPN" = "#00991F", "Cortex_MP_SCPN" = "#5C5CA6", "UL_CPN" = "#01545a", "Cortex_UL_CPN" = "#01545a",
"CTGL_ITL6GL" = "#66C5CC", "PTGL" = "#C9DB74", "ITL5GL" = "#87C55F", "ITL4GL" = "#54990F", "ITL23GL" = "#017351", "ITL1" = "#3DCCB1", "Astrocytes" = "#16F2F2", "CR" = "#F2F318", "Cortex_CR" = "#F2F318", "MGE_derived_InN" = "#D4E727", "Cortex_MGE_derived_InN" = "#D4E727",
"ODC1" = "#9EFF99", "ODC2" = "#64C2FC", "SSTGA" = "#CC79A7", "PVGA" = "#FF26CB", "VIPGA" = "#267DFF", 
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
#!-----------------------------------------------------------------------------------------------------------
meta_data <- read.table("../../data/scRNA_Data/metadata/metaData_scDevSC.txt", header = T, row.names = 1, sep = "\t")
data_matrix <- readMM("../../data/scRNA_Data/expression/gene_sorted-matrix.mtx")

barcodes <- read.table("../../data/scRNA_Data/expression/barcodes.tsv")[,1]
gene_names <- read.table("../../data/scRNA_Data/expression/genes.tsv")
rownames(data_matrix) <- gene_names$V2
colnames(data_matrix) <- barcodes

count_matrix <- round(expm1(data_matrix))

seurat_object <- CreateSeuratObject(counts = count_matrix, data = data_matrix, meta.data = meta_data)

umap <- read.table("../../data/scRNA_Data/cluster/cluster_scDevSC.merged.umap.txt", header = T, row.names = 1)
seurat_object[['umap']] <- CreateDimReducObject(embeddings = as.matrix(umap))

seurat_object$New_cellType <- gsub(seurat_object$New_cellType, pattern = " ", replacement = "_")
seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[-which(seurat_object$New_cellType %in% c("Doublet", "Red_blood_cells", "Low_quality_cells"))])
seurat_object$New_cellType <- factor(seurat_object$New_cellType, levels = c("Apical_progenitors", "Intermediate_progenitors", "Migrating_neurons", "Immature_neurons", "Cajal_Retzius_cells", "CThPN", "SCPN", "NP", "Layer_6b", "Layer_4", "DL_CPN", "DL_CPN_1", "DL_CPN_2", "UL_CPN", "Interneurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Cycling_glial_cells", "Ependymocytes", "Endothelial_cells", "VLMC", "Pericytes"))

colnames(seurat_object@meta.data)[17] <- "cell_identity"
Idents(seurat_object) <- "cell_identity"

saveRDS(seurat_object, "output/brain_mouse_scRNA_2021_Nature.rds")
dir.create("cell2location_Res")
dir.create("cell2location_Res/anndata")

seurat2scanpy(seurat_object, savefile = "cell2location_Res/anndata/scRNA_ref.h5ad")
#!-----------------------------------------------------------------------------------------------------------
gene.names <- read.table("../../data/Mus_musculus.GRCm38.102.gene.id_name.tsv", sep = "\t", header = T)

# remove unnecessary pseudogene, rRNA, snoRNA, tRNA, ribozyme
removed.genes <- gene.names$gene.name[-grep("protein_coding", gene.names$gene.bio_type)]

meta_data <- read.csv("../../data/scRNA_Data/GSE224560_FC_Droplet_Paired-Tag_metadata/meta_data.csv", row.names = 1)
meta_data <- meta_data[-which(meta_data$Anno0.8_yel %in% c('Doublet (MGL / ITLGL)', 'Doublet (OGC / OBGA)', 'Unknown (MGL / OGC)')),]
meta_data <- meta_data[,c(4:7, 17:20, 22, 23)]
colnames(meta_data)[5] <- "cell_identity"

# deal Droplet-Paired-Tag FC data
sample_list <- list.files("../../data/scRNA_Data/GSE224560_RAW/")
samples <- gsub(pattern = 'GSM[0-9]*_|_raw_feature_bc_matrix.h5', replacement = '', x = sample_list)

seurat_list <- list()
count = 1

for(i in 1:length(sample_list)){
	count.matrix <- Read10X_h5(paste0("../../data/scRNA_Data/GSE224560_RAW/", sample_list[i]), unique.features = FALSE)
	#rownames(count.matrix) <- gene.names$V2[match(rownames(count.matrix), gene.names[,1])]
	#count.matrix <- avereps(count.matrix) # average of the same gene name
	colnames(count.matrix) <- paste0(samples[i], "-", colnames(count.matrix))
	print(dim(count.matrix))
	count.matrix <- count.matrix[,rownames(meta_data)[which(meta_data$sample == samples[i])]]
	ID <- as.character(rownames(count.matrix))
	d <- duplicated(ID)
	ID <- factor(ID, levels = unique(ID))
	count.matrix <- rowsum(as.matrix(count.matrix), ID, reorder = FALSE, na.rm = TRUE)
	print(dim(count.matrix))
	removed.gene <- intersect(rownames(count.matrix), removed.genes)
	count.matrix <- count.matrix[-match(removed.gene, rownames(count.matrix)),]
	print(dim(count.matrix))

	temp_object <- CreateSeuratObject(counts = count.matrix, project = samples[i], min.cells = 1, min.features = 1)
	temp_object[["percent.mt"]] <- PercentageFeatureSet(object = temp_object, pattern = "^mt")
	temp_object[["percent.ribo"]] <- PercentageFeatureSet(object = temp_object, pattern = "^Rp[sl]")
	seurat_list[count] <- temp_object
	count <- count + 1
}
names(seurat_list) <- samples
seurat_object <- merge(seurat_list[[1]], seurat_list[2:length(seurat_list)])
all(colnames(seurat_object) == rownames(meta_data))
seurat_object@meta.data <- cbind(seurat_object@meta.data, meta_data)
seurat_object <- JoinLayers(seurat_object)

seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object)
seurat_object <- RunUMAP(seurat_object, dims = 1:30, return.model = T)

library(SeuratWrappers)
library(harmony)
seurat_object <- RunHarmony(seurat_object, group.by.vars = "orig.ident",  reduction.name = 'pca', reduction.save = 'harmony')
seurat_object <- RunUMAP(seurat_object, reduction = "harmony", dims = 1:30, reduction.name = "harmony_umap", return.model = T)

seurat_object$cell_identity <- factor(seurat_object$cell_identity, levels = c("ASC", "OPC", "OGC", "VLMC", "MGL", "VIPGA", "PVGA", "SSTGA", "OBGA", "D12MSN", "STRGA", "OBGL", "CLAGL", "ITL23GL", "PTGL", "ITL45GL", "ITL5GL", "ITL6GL", "CTGL", "NPGL"))

p0 <- DimPlot(seurat_object, reduction = "harmony_umap", group.by = c("cell_identity"), label = T, cols = unique(seurat_object$Subclass.color))
p1 <- DimPlot(seurat_object, reduction = "umap", group.by = c("cell_identity"), label = T, cols = unique(seurat_object$Subclass.color))
p0 + p1

saveRDS(seurat_object, "output/FC_mouse_snRNA_2023_NSMB.rds")
seurat2scanpy(seurat_object, savefile = "cell2location_Res/anndata/scRNA_Adult_FC_ref.h5ad")
#!-----------------------------------------------------------------------------------------------------------
#! in python
# import scanpy as sc
# import numpy as np
# import pandas as pd
# from scipy.sparse import csr_matrix, csc_matrix
# import scipy.io as sio
# adata = sc.read_h5ad("/data/test_data/Adult_brain_sc/WMB-10Xv3-CB-raw.h5ad")
# adata.X = adata.X.astype(np.int32)
# adata.X = adata.X.tocsc()
# sio.mmwrite('WMB-10Xv3-CB-raw.mtx', adata.X)
#
#samples = ['WMB-10Xv3-CB-raw.h5ad', 'WMB-10Xv3-CTXsp-raw.h5ad', 'WMB-10Xv3-HPF-raw.h5ad', 'WMB-10Xv3-HY-raw.h5ad', 'WMB-10Xv3-Isocortex-1-raw.h5ad', 'WMB-10Xv3-Isocortex-2-raw.h5ad', 'WMB-10Xv3-MB-raw.h5ad', 'WMB-10Xv3-MY-raw.h5ad', 'WMB-10Xv3-OLF-raw.h5ad', 'WMB-10Xv3-P-raw.h5ad', 'WMB-10Xv3-PAL-raw.h5ad', 'WMB-10Xv3-STR-raw.h5ad', 'WMB-10Xv3-TH-raw.h5ad']
#for i in samples:
#	adata = sc.read_h5ad("/data/test_data/Adult_brain_sc/" + i)
#	adata.X = adata.X.astype(np.int32)
#	adata.X = adata.X.tocsc()
#	sio.mmwrite( + '.mtx', adata.X)
#	adata.obs.to_csv("/data/zhanglin/Projects/brain_project/" + i + ".barcode")
#	adata.var.to_csv("/data/zhanglin/Projects/brain_project/" + i + ".features")

meta_data <- read.csv("/data/test_data/Adult_brain_sc/cell_type_metadata/cell_metadata.csv")
meta_data <- meta_data[which(meta_data$library_method == "10Xv3"),]

h5ad_files <- list.files("/data/test_data/Adult_brain_sc")
h5ad_files <- h5ad_files[grep(".h5ad", h5ad_files)]
mat_files <- gsub(h5ad_files, pattern = ".h5ad", replacement = '')

genes_info <- read.csv("/data/test_data/Adult_brain_sc/cell_type_metadata/gene.csv")
genes_info <- genes_info[-which(genes_info$mapped_ncbi_identifier == ''),]

seurat_list <- list()
for(i in 2:length(h5ad_files)){
	count_matrix <- t(Matrix::readMM(paste0(mat_files[i], ".mtx")))
	barcodes <- read.csv(paste0(mat_files[i], ".barcode"))
	features <- read.csv(paste0(mat_files[i], ".features"))
	rownames(count_matrix) <- features$gene_symbol
	colnames(count_matrix) <- barcodes$cell_label
	count_matrix = as(count_matrix, "CsparseMatrix")
	count_matrix <- count_matrix[rownames(count_matrix) %in% genes_info$gene_symbol,]
	
	tmp_names <- as.data.frame(table(rownames(count_matrix)))
	unique_genes <- as.character(tmp_names$Var1[which(tmp_names$Freq == 1)])
	tmp_matrix <- count_matrix[-match(unique_genes, rownames(count_matrix)),]
	count_matrix <- count_matrix[match(unique_genes, rownames(count_matrix)),]
	ID <- as.character(rownames(tmp_matrix))
	ID <- factor(ID, levels = unique(ID))
	tmp_matrix <- as(rowsum(as.matrix(tmp_matrix), ID, reorder = FALSE, na.rm = TRUE), "dgCMatrix")
	count_matrix <- rbind(count_matrix, tmp_matrix)
	#ID <- as.character(rownames(count_matrix))
	#ID <- factor(ID, levels = unique(ID))
	#count_matrix <- rowsum(as.matrix(count_matrix), ID, reorder = FALSE, na.rm = TRUE)
	print(dim(count_matrix))
	temp_object <- CreateSeuratObject(counts = count_matrix)
	seurat_list[[i]] <- temp_object
}
seurat_object <- merge(seurat_list[[1]], seurat_list[2:length(seurat_list)])
tmp <- intersect(colnames(seurat_object), meta_data$cell_label)
meta_data <- meta_data[match(tmp, meta_data$cell_label),]

seurat_object <- subset(seurat_object, cells = tmp)
seurat_object@meta.data <- cbind(seurat_object@meta.data, meta_data)
umap_matrix <- as.matrix(seurat_object@meta.data[,c("x", "y")])
colnames(umap_matrix) <- c("UMAP_1", "UMAP_2")
seurat_object[['umap']] <- CreateDimReducObject(embeddings = umap_matrix)

cluster_meta_data <- read.table("cluster_meta_data", sep = "\t", header = T, check.names = F)
tmp <-cluster_meta_data[match(tmp$cluster_alias, cluster_meta_data$cl),]
all(tmp$cl == seurat_object$cluster_alias)
seurat_object@meta.data <- cbind(seurat_object@meta.data, tmp[,2:47])

qsave(seurat_object, "seurat_object.qs")

selected_samples <- c("CTXsp", "HPF", "Isocortex-1", "Isocortex-2", "STR", "PAL")
selected_cells <- c()
for(i in paste0("WMB-10Xv3-", selected_samples, "-raw")){
	barcodes <- read.csv(paste0(i, ".barcode"))
	selected_cells <- c(selected_cells, barcodes$cell_label)
}
selected_cells <- intersect(selected_cells, colnames(seurat_object))
seurat_object <- subset(seurat_object, cells = selected_cells)

tmp <- as.data.frame(table(seurat_object$library_label))
tmp <- tmp[order(tmp$Freq),]

tmp <- as.data.frame(table(seurat_object$subclass_label))
tmp <- tmp[order(tmp$Freq),]

selected_subclass_label <- tmp$Var1[which(tmp$Freq >= 100)]
seurat_object <- subset(seurat_object, subset = subclass_label %in% selected_subclass_label)

selected_cells <- WhichCells(seurat_object, downsample = 1000)
seurat_object <- subset(seurat_object, cells = selected_cells)

library(reticulate)
use_python("/home/chenlab/softwares/miniconda3/envs/Seurat5/bin")

sc <- import("scanpy")
sq <- import("squidpy")
np <- import("numpy")
pd <- import("pandas")
ad <- import("anndata")

meta_data <- seurat_object@meta.data
for(i in 1:ncol(meta_data)){
	if(class(meta_data[,i]) == "factor"){
		meta_data[,i] <- as.character(meta_data[,i])
	}
}

data_matrix <- GetAssayData(seurat_object, assay = "RNA", slot = "counts")
adata = ad$AnnData(X = r_to_py(t(data_matrix)))
adata$obs = r_to_py(seurat_object@meta.data)
rownames(adata$var) = rownames(data_matrix)
adata$obsm[paste0("X_umap")] = r_to_py(np$array(seurat_object@reductions[["umap"]]@cell.embeddings))
adata$write(savefile)
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

seurat_list[[1]] <- subset(seurat_list[[1]], subset = nCount_RNA < 12000 & nFeature_RNA < 5000 & percent.mt < 5 & percent.ribo < 3)
seurat_list[[2]] <- subset(seurat_list[[2]], subset = nCount_RNA < 25000 & nFeature_RNA < 7000 & percent.mt < 5 & percent.ribo < 3)
seurat_list[[3]] <- subset(seurat_list[[3]], subset = nCount_RNA < 20000 & nFeature_RNA < 6000 & percent.mt < 5 & percent.ribo < 3)
seurat_list[[4]] <- subset(seurat_list[[4]], subset = nCount_RNA < 10000 & nFeature_RNA < 5000 & percent.mt < 5 & percent.ribo < 3)
seurat_list[[5]] <- subset(seurat_list[[5]], subset = nCount_RNA < 10000 & nFeature_RNA < 5000 & percent.mt < 5 & percent.ribo < 3)

saveRDS(seurat_list, "output/seurat_list.rds")

dir.create("SCT_vs_banksy_GraphPCA_SpatialPCA_Res")
setwd("SCT_vs_banksy_GraphPCA_SpatialPCA_Res")

library(SpatialPCA)
library(SeuratWrappers)
library(Banksy)

library(reticulate)
use_python("/home/yiyelinfeng/softwares/miniconda3/envs/R43/bin")

sc <- import("scanpy")
sq <- import("squidpy")
np <- import("numpy")
pd <- import("pandas")
sg <- import("GraphPCA")
ad <- import("anndata")

cc.genes <- cc.genes.updated.2019
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
gene.mouse <- read.csv("/media/yiyelinfeng/data/Repository/baseFiles/ortholog_mouse.csv", header = T)

for(i in 1:length(seurat_list)){
	tmp_object <- seurat_list[[i]]
	DefaultAssay(tmp_object) <- "RNA"
	tmp_object <- SCTransform(tmp_object, vst.flavor = "v2", method = "glmGamPoi")

	s.gene <- intersect(rownames(tmp_object), s.genes)
	g2m.gene <- intersect(rownames(tmp_object), g2m.genes)
	s.gene <- intersect(rownames(tmp_object), gene.mouse$Mouse.gene.name[match(s.genes, gene.mouse$Gene.name)])
	g2m.gene <- intersect(rownames(tmp_object), gene.mouse$Mouse.gene.name[match(g2m.genes, gene.mouse$Gene.name)])
	
	tmp_object <- CellCycleScoring(tmp_object, s.features = s.gene, g2m.features = g2m.gene, set.ident = FALSE, nbin = 12)
	tmp_object$CC.Difference <- tmp_object$S.Score - tmp_object$G2M.Score

	tmp_object <- RunPCA(tmp_object, assay = "SCT", reduction.name = "pca", npcs = 30)
	tmp_object <- RunUMAP(tmp_object, reduction = "pca", dims = 1:30, min.dist = 0.3, reduction.name = "SCT_umap", return.model = T)
	tmp_object <- FindNeighbors(tmp_object, reduction = "pca", dims = 1:30, annoy.metric = "cosine", graph.name = c('SCT_nn', 'SCT_snn'))
	tmp_object <- FindClusters(tmp_object, cluster.name = "SCT_cluster", resolution = 0.5, graph.name = "SCT_snn")
	
	X <- GetAssayData(tmp_object, assay = "SCT", slot = "scale.data")
	locs <- Seurat::GetTissueCoordinates(tmp_object)[,seq_len(2)]
	colnames(locs) <- c("array_row", "array_col")
	tmp_object@meta.data <- cbind(tmp_object@meta.data, locs)

	sp_count <- GetAssayData(tmp_object, assay = "RNA", layer = "counts")
	sparkversion <- ifelse(ncol(sp_count) > 3000, "sparkx", "spark")
	temp_object <- CreateSpatialPCAObject(counts = sp_count, location = as.matrix(locs), project = "SpatialPCA", gene.type = "spatial", sparkversion = sparkversion, numCores_spark = 5, customGenelist = NULL, min.loctions = 20, min.features = 20)
	bandwidthtype <- ifelse(ncol(sp_count) > 3000, "Silverman", "SJ")
	temp_object <- SpatialPCA_buildKernel(temp_object, kerneltype = "gaussian", bandwidthtype = bandwidthtype)
	temp_object <- SpatialPCA_EstimateLoading(temp_object, fast = FALSE, SpatialPCnum = 30)
	temp_object <- SpatialPCA_SpatialPCs(temp_object, fast = FALSE)
	tmp_object[["spPCA"]] <- CreateAssay5Object(data = temp_object@normalized_expr)
	SpatialPCA <- t(temp_object@SpatialPCs)
	colnames(SpatialPCA) <- paste0("spPCA", 1:30)
	tmp_object[['SpatialPCA']] <- CreateDimReducObject(embeddings = SpatialPCA, assay = "spPCA", key = "spPCA_")
	
	tmp_object <- RunUMAP(tmp_object, reduction = "SpatialPCA", dims = 1:30, min.dist = 0.3, reduction.name = "spPCA_umap", return.model = T)
	tmp_object <- FindNeighbors(tmp_object, dims = 1:30, reduction = 'SpatialPCA', annoy.metric = "cosine", graph.name = c('spPCA_nn', 'spPCA_snn'))
	tmp_object <- FindClusters(tmp_object, resolution = 0.5, cluster.name = "spPCA_cluster", graph.name = 'spPCA_snn')
	
	adata = ad$AnnData(X = t(X))
	adata$obs = tmp_object@meta.data
	adata$obs$index = colnames(tmp_object)
	adata$var$index = rownames(X)
	adata$var["highly_variable"] = rownames(X)
	locs <- Seurat::GetTissueCoordinates(tmp_object)[,seq_len(2)]
	location <- as.matrix(t(locs))
	GraphPCA_Res <- sg$Run_GPCA(adata, location = locs, n_components = as.integer(30), method = "knn", platform = "Visium", save_reconstruction = TRUE)

	GraphPCA <- GraphPCA_Res[[1]]
	colnames(GraphPCA) <- paste0("GraphPCA", 1:30)
	rownames(GraphPCA) <- colnames(tmp_object)
	tmp_object[['GraphPCA']] <- CreateDimReducObject(embeddings = GraphPCA, assay = "RNA", key = "GraphPCA_")
	tmp_object <- RunUMAP(tmp_object, reduction = "GraphPCA", dims = 1:30, min.dist = 0.3, reduction.name = "GraphPCA_umap", return.model = T)
	tmp_object <- FindNeighbors(tmp_object, reduction = "GraphPCA", dims = 1:30, annoy.metric = "cosine", graph.name = c('GraphPCA_nn', 'GraphPCA_snn'))
	tmp_object <- FindClusters(tmp_object, cluster.name = "GraphPCA_cluster", resolution = 0.3, graph.name = "GraphPCA_snn")

	tmp_object <- RunBanksy(tmp_object, lambda = 0.2, assay = 'SCT', slot = 'data', dimx = "array_row", dimy = "array_col", features = 'variable', k_geom = 15)
	DefaultAssay(tmp_object) <- "BANKSY"
	tmp_object <- RunPCA(tmp_object, assay = "BANKSY", reduction.name = "banksy_pca", features = rownames(tmp_object), npcs = 30)
	tmp_object <- RunUMAP(tmp_object, reduction = "banksy_pca", dims = 1:30, min.dist = 0.3, reduction.name = "banksy_umap", return.model = T)
	tmp_object <- FindNeighbors(tmp_object, reduction = "banksy_pca", dims = 1:30, annoy.metric = "cosine", graph.name = c('banksy_nn', 'banksy_snn'))
	tmp_object <- FindClusters(tmp_object, cluster.name = "banksy_cluster", resolution = 0.3, graph.name = "banksy_snn")
		
	p1 <- DimPlot(tmp_object, reduction = "SCT_umap", group.by = "SCT_cluster", cols = selected_colors, label = T) + NoLegend()
	p2 <- DimPlot(tmp_object, reduction = "banksy_umap", group.by = "banksy_cluster", cols = selected_colors, label = T) + NoLegend()
	p3 <- DimPlot(tmp_object, reduction = "GraphPCA_umap", group.by = "GraphPCA_cluster", cols = selected_colors, label = T) + NoLegend()
	p4 <- DimPlot(tmp_object, reduction = "spPCA_umap", group.by = "spPCA_cluster", cols = selected_colors, label = T) + NoLegend()
	print((p1 + p2)/(p3 + p4))
	ggsave(paste0("SCT_vs_banksy_vs_GraphPCA_SpatialPCA_for_", names(seurat_list)[i], "_umap.png"), width = 14, height = 14)
	p1 <- Seurat::SpatialPlot(tmp_object, image.scale = "hires", shape = 22, group.by = "SCT_cluster", cols = selected_colors)
	p2 <- Seurat::SpatialPlot(tmp_object, image.scale = "hires", shape = 22, group.by = "banksy_cluster", cols = selected_colors)
	p3 <- Seurat::SpatialPlot(tmp_object, image.scale = "hires", shape = 22, group.by = "GraphPCA_cluster", cols = selected_colors)
	p4 <- Seurat::SpatialPlot(tmp_object, image.scale = "hires", shape = 22, group.by = "spPCA_cluster", cols = selected_colors)
	print((p1 + p2)/(p3 + p4))
	ggsave(paste0("SCT_vs_banksy_vs_GraphPCA_SpatialPCA_for_", names(seurat_list)[i], "_Spatial.png"), width = 16, height = 14)
	seurat_list[[i]] <- tmp_object
}
saveRDS(seurat_list, "../output/seurat_list_SCT_banksy_GraphPCA_spatialPCA.rds")
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

tmp_object$cell_identity <- as.character(seurat_object$cell_identity[match(paste0("E14_5_", colnames(tmp_object)), colnames(seurat_object))])
DimPlot(tmp_object, reduction = "SCT_umap", group.by = "cell_identity", cols = cluster_colors, label = T)
tmp_object <- FindClusters(tmp_object, cluster.name = "SCT_cluster", resolution = 1.7, graph.name = "SCT_snn")

tmp_object$cell_identity <- as.character(tmp_object$SCT_cluster)
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(2))] <- "Cortex_DP_AP"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(14))] <- "Cortex_MP_AP"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(8))] <- "Cortex_IPC"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(6))] <- "Cortex_DP_MigN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(22))] <- "Cortex_MP_MigN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(21))] <- "Cortex_CThPN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(12, 16))] <- "Cortex_SCPN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(18))] <- "Cortex_CR"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(4))] <- "LGE_MGE_AP"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(23))] <- "MGE_IPC"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(9))] <- "MGE_InN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(13))] <- "LGE_IPC"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(3))] <- "LGE_InN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(5))] <- "CPU_pre_MSN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(0, 1))] <- "CPU_MSN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(24))] <- "LCS_IMN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(20))] <- "Epd_ExN" # "Tfap2d", "Tafa1", "Lhfpl3", "Epha3", "Epha6", "Ppp1r14c", "Thsd7a", "Fam135b", "Pcdh19", "Syn3", "Ndst3", "Elmod1", "Tafa5", "Il1rapl1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(26))] <- "CLA_ExN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(15))] <- "PIR_ExN" # "Tafa1", "Unc13c", "Adgrl3", "Tafa2", "Fstl4", "Sphkap", "Epha7", "Csmd1", "Cntn1", "Npas3", "Plppr5", "Lhfpl3", "Egfem1", "Dner", "Ptpre", "Htr2a", "Slc1a2"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(10))] <- "LS_AP"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(25))] <- "AEP_AP" # "Lrp2", "Bcan", "Pdlim3", "Slc6a11", "Tll2", "Shh", "Itga2", "Cdh19", "Sulf1", "Mamdc2", "Ednrb", "Fabp7", "Cybrd1", "A2m", "Mfge8", "Lhcgr"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(7))] <- "Endo"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(11, 19, 27))] <- "Fibrob"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(17))] <- "Stromal"

tmp_object$cell_identity <- factor(tmp_object$cell_identity, levels = c("Cortex_DP_AP", "Cortex_MP_AP", "Cortex_IPC", "Cortex_DP_MigN", "Cortex_MP_MigN", "Cortex_CThPN", "Cortex_SCPN", "Cortex_CR", "LGE_MGE_AP", "MGE_IPC", "MGE_InN", "LGE_IPC", "LGE_InN", "CPU_pre_MSN", "CPU_MSN", "LCS_IMN", "Epd_ExN", "PIR_ExN", "CLA_ExN", "LS_AP", "AEP_AP", "Endo", "Fibrob", "Stromal"))
Seurat::SpatialPlot(tmp_object, group.by = "cell_identity", shape = 22, cols = cluster_colors, crop = F)

tmp_object$cell_type <- as.character(tmp_object$cell_identity)
tmp_object$cell_type <- gsub(tmp_object$cell_type, pattern = "Cortex_", replacement = "")
tmp_object$cell_type <- gsub(tmp_object$cell_type, pattern = "DP_", replacement = "")
tmp_object$cell_type <- gsub(tmp_object$cell_type, pattern = "MP_", replacement = "")

tmp_object$cell_type <- factor(tmp_object$cell_type, levels = c("AP", "IPC", "MigN", "CThPN", "SCPN", "CR", "LGE_MGE_AP", "MGE_IPC", "MGE_InN", "LGE_IPC", "LGE_InN", "CPU_pre_MSN", "CPU_MSN", "LCS_IMN", "Epd_ExN", "PIR_ExN", "CLA_ExN", "LS_AP", "AEP_AP", "Endo", "Fibrob", "Stromal"))
Seurat::SpatialPlot(tmp_object, group.by = "cell_type", shape = 22, cols = cluster_colors, crop = F)
#!------------------------------
Idents(tmp_object) <- "cell_identity"
all_markers <- FindAllMarkers(tmp_object, min.pct = 0.1, only.pos = T)
de_markers <- all_markers[which(all_markers$avg_log2FC > 0.5 & all_markers$p_val_adj < 0.01 & all_markers$pct.1 > 0.3),]
write.csv(de_markers, "E14_5_ST_cell_identity_DE_markers.csv", quote = F)

key_markers <- c("Satb2", "Top2a", "Mki67", "Neurog2", "Eomes", "Hmga2", "Pax6", "Dct", "Adamts19", "Emid1", "Slc17a6", "Nrn1", "Neurod1", "Adamtsl3", "Cdh12", "Neurod2", "Neurod6", "Cntn2", "Abcc8", "Nfe2l3", "Trpc3", "Rgs6", "Dync1i1", "Opcml", "Nr4a3", "Fezf2", "Reln", "Dpp10", "Trp73",  "Zhx2","Vit", "Ednrb", "Ildr2", "Sp9", "Dlx1", "Dlx6", "Stk33", "Nkx2-1", "Lhx6", "Lhx8", "Rai2", "Esrrg", "Rffl", "Slc18a2", "Dlx2", "Gucy1a1", "Zfp503", "Isl1", "Rarb", "Oprm1", "Frmd6", "Ikzf1", "Grid1", "Gad2", "Gad1", "Tfap2d", "Nr2f1", "Tafa2", "Tafa1", "Lhfpl3", "Unc13c", "Fstl4", "Sphkap", "Plppr5", "Nr4a2", "Galnt14", "Zic4", "Otx2", "Fgf17", "Lrp2", "Bcan", "Slc6a11", "Kdr", "Flt1", "Ptprb", "Ahnak", "Col5a1", "Dcn", "Bmp6", "Ranbp3l", "Itih5")

DotPlot(tmp_object, features = key_markers, col.min = 0.1, dot.min = 0.1) + RotatedAxis() + labs(y = "cell identity", x = "Features") + scale_color_gradientn(colours = paletteContinuous(set = "whitePurple"))
ggsave("Development_QC_cell_identity/E14_5_cell_identity_key_markers.pdf")
#!------------------------------
tmp_object$Cluster <- names(ClusterNames)[match(as.character(tmp_object$cell_identity), ClusterNames)]
tmp_object$Cluster <- factor(tmp_object$Cluster, levels = intersect(ClusterID, unique(tmp_object$Cluster)))
Seurat::SpatialPlot(tmp_object, group.by = "Cluster", shape = 22, cols = ClusterIDColors)

tmp_object <- FindClusters(tmp_object, resolution = seq(0.1, 1, 0.1), graph.name = "banksy_snn")

wrap_plots(map(seq(0.1, 1, 0.1), function(x) DimPlot(tmp_object, reduction = "banksy_umap", group.by = paste0("banksy_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)
Seurat::SpatialPlot(tmp_object, shape = 22, cols = selected_colors, group.by = paste0("banksy_snn_res.", seq(0.1, 1, 0.1)), ncol = 5) & NoLegend()
tmp_object$banksy_cluster <- as.character(tmp_object$banksy_snn_res.0.2)
tmp_object$banksy_cluster[which(tmp_object$banksy_snn_res.0.4 == 13)] <- 8
tmp_object$banksy_cluster <- paste0("D", tmp_object$banksy_cluster)

tmp_object$banksy_domains <- as.character(tmp_object$banksy_cluster)
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D4"))] <- "Cortex_1"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D1"))] <- "Cortex_2"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D2"))] <- "Cortex_3"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D3"))] <- "Cortex_5"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D8"))] <- "LGE0"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D6"))] <- "LGE"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D0"))] <- "CPU"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D5"))] <- "MGE"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D7"))] <- "LS"
tmp_object$banksy_domains <- factor(tmp_object$banksy_domains, levels = c("Cortex_1", "Cortex_2", "Cortex_3", "Cortex_5", "LGE0", "LGE", "CPU", "MGE", "LS"))
tmp_object$banksy_cluster <- domainID[match(as.character(tmp_object$banksy_domains), domainNames)]
tmp_object$banksy_cluster <- factor(tmp_object$banksy_cluster, levels = intersect(domainID, unique(tmp_object$banksy_cluster)))
Seurat::SpatialPlot(tmp_object, group.by = "banksy_cluster", shape = 22, cols = domainIDColors, crop = FALSE)

# cell2location analysis
cell2location <- read.csv("cell2location_Res/development_FC_cell2location/sinlge_timepoint/ST_E14_5_sp_cell2location.csv", row.names = 1)
all(rownames(cell2location) == colnames(tmp_object))
colnames(cell2location)[18:27] <- paste0("E14_5_", colnames(cell2location)[18:27])
tmp_object@meta.data <- cbind(tmp_object@meta.data, cell2location[,18:27])
cell2location <- read.csv("cell2location_Res/development_FC_cell2location/timepoints/ST_E14_5_sp_cell2location.csv", row.names = 1, check.names = F)
all(rownames(cell2location) == colnames(tmp_object))
tmp_object@meta.data <- cbind(tmp_object@meta.data, cell2location[,18:36])
#!-----------------------------------------------------------------
temp_object <- readRDS("output/E14_brain_mouse_scRNA_2021_Nature.rds")

DefaultAssay(tmp_object) <- "RNA"
tmp_object <- NormalizeData(tmp_object)
anchors <- FindTransferAnchors(reference = temp_object, query = tmp_object, dims = 1:30, reference.reduction = "pca", reference.assay= "RNA", query.assay = "RNA", reduction = "pcaproject")
tmp_object <- MapQuery(anchorset = anchors, reference = temp_object, query = tmp_object, refdata = list(celltype = "cell_identity"), reference.reduction = "pca", new.reduction.name = "ref.pca", reduction.model = "umap", projectumap.args = list(reduction.name = "ref.umap"))
colnames(tmp_object@meta.data)[81] <- "predicted_E14_5_celltype_2021_Nature.score"
colnames(tmp_object@meta.data)[82] <- "predicted_E14_5_celltype_2021_Nature"

tmp_object$predicted_E14_5_celltype_2021_Nature <- factor(tmp_object$predicted_E14_5_celltype_2021_Nature, levels = intersect(levels(temp_object$cell_identity), unique(tmp_object$predicted_E14_5_celltype_2021_Nature)))

predictions <- table(tmp_object$Cluster, tmp_object$predicted_E14_5_celltype_2021_Nature)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)
predictions$Var1 <- paste0(predictions$Var1, "_", ClusterNames[match(predictions$Var1, ClusterID)])
predictions$Var1 <- factor(predictions$Var1, levels = intersect(paste0(ClusterID, "_", ClusterNames), unique(predictions$Var1)))

p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
    low = "#ffffc8", high = "#7d0025") + xlab("SCT clusters") + ylab("Predicted cell type label with 2021 Nature E14.5") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
print(p1)
ggsave("predicted_E14_5_celltype_2021_Nature_Clusters_all.pdf")

predictions <- predictions[grep("Cortex_|C30_Endo|C31_Fibrob|C33_Stromal|MGE_|LGE_IPC", predictions$Var1),]

p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
    low = "#ffffc8", high = "#7d0025") + xlab("SCT clusters") + ylab("Predicted cell type label with 2021 Nature E14.5") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
print(p1)
ggsave("predicted_E14_5_celltype_2021_Nature_Clusters.pdf")

DefaultAssay(tmp_object) <- "SCT"
Idents(tmp_object) <- "cell_identity"

Seurat::SpatialPlot(tmp_object, shape = 22, cols = ClusterIDColors, group.by = "Cluster", crop = FALSE)
Seurat::SpatialPlot(tmp_object, shape = 22, cols = domainIDColors, group.by = "banksy_cluster")
Seurat::SpatialPlot(tmp_object, shape = 22, cols = cluster_colors)
DimPlot(tmp_object, cols = ClusterIDColors, group.by = "Cluster", label = T, repel = T) + SetAxes() + NoLegend()
DimPlot(tmp_object, cols = cluster_colors, label = T, repel = T) + SetAxes()
DimPlot(tmp_object, reduction = "banksy_umap", cols = domainIDColors, group.by = "banksy_cluster", label = T, repel = T) + SetAxes()

library(ggsankey)

meta_data <- tmp_object@meta.data %>% make_long(Cluster, banksy_cluster)
meta_data$node <- factor(meta_data$node, levels = c(levels(tmp_object$Cluster), intersect(domainID, unique(meta_data$node))))
meta_data$next_node <- factor(meta_data$next_node, levels = c(levels(tmp_object$Cluster), intersect(domainID, unique(meta_data$next_node))))

ggplot(meta_data, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) + geom_sankey(color = "black", lwd = 0.2, flow.alpha = 0.5, show.legend = FALSE, na.rm = FALSE) + scale_x_discrete(expand = c(0, 0.8)) + 
		scale_fill_manual(values = c(ClusterIDColors, domainIDColors), drop = FALSE) + geom_sankey_text(size = 3.2, color = "black") + theme_void() + theme(axis.text.x = element_text())
ggsave("Development_QC_cell_identity/E14_5_SCT_spot_cluster_banksy_domian_sankey_plot.pdf")
#!---------------------------------------------------------------------------------------------------------------
library(UCell)
source("/home/yiyelinfeng/scripts/Rscripts/lung_project/IPF/spatial-lung-fibrosis/scripts/custom_colors.R")

selected_clusters <- paste0("E14_5_", levels(temp_object$cell_identity))

Idents(tmp_object) <- "Cluster"
DefaultAssay(tmp_object) <- "SCT"
all_markers <- FindAllMarkers(tmp_object, only.pos = T, min.pct = 0.1)

de_markers <- all_markers[which(all_markers$avg_log2FC > 0.5 & all_markers$p_val_adj < 0.01 & all_markers$pct.1 > 0.3),]
marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:30]
	count = count + 1
}
names(marker_genes) <- unique(de_markers$cluster)

tmp_object <- UCell::AddModuleScore_UCell(tmp_object, features = marker_genes, ncores = 1, name = "")

cluster_c2l_density <- cor(tmp_object@meta.data[,which(colnames(tmp_object@meta.data) %in% c(levels(tmp_object$Cluster), selected_clusters))])
diag(cluster_c2l_density) <- 0
cluster_c2l_density <- as.data.frame(cluster_c2l_density)
cluster_c2l_density <- cluster_c2l_density[setdiff(rownames(cluster_c2l_density), levels(tmp_object$Cluster)) , levels(tmp_object$Cluster)]
colnames(cluster_c2l_density) <- paste0(colnames(cluster_c2l_density), "_", ClusterNames[match(colnames(cluster_c2l_density), ClusterID)])
cluster_c2l_density <- cluster_c2l_density[selected_clusters, paste0(levels(tmp_object$Cluster), "_", levels(tmp_object$cell_identity))]

pheatmap::pheatmap(cluster_c2l_density, cellwidth = 15, cellheight = 15, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, treeheight_col = 20, treeheight_row = 20, color = col_scale_div_custom2, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_Nat_2021_", names(seurat_list)[1], "_all.pdf"), width = 10, height = 10)

cluster_c2l_density <- cluster_c2l_density[selected_clusters, grep("Cortex_|C33_Endo|C34_Fibrob|C36_Stromal|MGE_|LGE_", colnames(cluster_c2l_density))]

pheatmap::pheatmap(cluster_c2l_density, cellwidth = 15, cellheight = 15, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, treeheight_col = 20, treeheight_row = 20, color = col_scale_div_custom2, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_Nat_2021_", names(seurat_list)[1], ".pdf"), width = 10, height = 10)

tmp_meta <- as.matrix(tmp_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(tmp_object$Cluster, "_", tmp_object$cell_identity)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)
tmp_meta <- tmp_meta[paste0(levels(tmp_object$Cluster), "_", levels(tmp_object$cell_identity)), selected_clusters]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[1], "_zscore_all.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[1], "_zscore_all1.pdf"), width = 10, height = 10)

tmp_meta <- as.matrix(tmp_object@meta.data[,selected_clusters[1:7]])
rownames(tmp_meta) <- paste0(tmp_object$Cluster, "_", tmp_object$cell_identity)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta <- tmp_meta[grep("Cortex_", rownames(tmp_meta)),]
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)
tmp_meta0 <- tmp_meta0[intersect(paste0(levels(tmp_object$Cluster), "_", levels(tmp_object$cell_identity)), rownames(tmp_meta0)),]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[1], "_zscore0.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[1], "_zscore01.pdf"), width = 10, height = 10)

selected_clusters <- c("Apical_progenitors", "Intermediate_progenitors", "Migrating_neurons", "Immature_neurons", "CThPN", "SCPN", "NP", "Layer_6b", "Layer_4", "DL_CPN", "UL_CPN", "Cajal_Retzius_cells", "Interneurons", "Astrocytes", "Oligodendrocytes", "Cycling_glial_cells", "Endothelial_cells", "Pericytes", "Microglia")

tmp_meta <- as.matrix(tmp_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(tmp_object$Cluster, "_", tmp_object$cell_identity)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)
tmp_meta0 <- tmp_meta0[paste0(levels(tmp_object$Cluster), "_", levels(tmp_object$cell_identity)), selected_clusters]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[1], "_zscore_all_times.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[1], "_zscore_all1_times.pdf"), width = 10, height = 10)

tmp_meta <- as.matrix(tmp_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(tmp_object$Cluster, "_", tmp_object$cell_identity)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta <- tmp_meta[grep("Cortex_|C33_Endo|C34_Fibrob|C36_Stromal|MGE_|LGE_", rownames(tmp_meta)),]
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)
tmp_meta0 <- tmp_meta0[intersect(paste0(levels(tmp_object$Cluster), "_", levels(tmp_object$cell_identity)), rownames(tmp_meta0)), selected_clusters]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[1], "_zscore_times.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[1], "_zscore1_times.pdf"), width = 10, height = 10)

Idents(tmp_object) <- "cell_identity"
seurat_list$E14_5 <- tmp_object
#!---------------------------------------------------------------------------------------------------------------
DimPlot(temp_object, group.by = "cell_identity", cols = ArchRPalettes$kelly, label = T, repel = T) + SetAxes()

tmp_object@assays <- tmp_object@assays[1]
tmp_object@meta.data <- tmp_object@meta.data[,c(1:5, 8:15, 18, 19, 40:43, 83, 84)]

coembed <- merge(x = temp_object, y = tmp_object)
coembed <- JoinLayers(coembed)
coembed@assays$RNA@layers <- coembed@assays$RNA@layers[1:2]

embedding <- tmp_object@reductions$ref.umap@cell.embeddings
embeddings <- rbind(temp_object@reductions$umap@cell.embeddings, embedding)
all(rownames(embeddings) == colnames(coembed))

coembed[["ref.umap"]] <- CreateDimReducObject(embeddings = embeddings, key = "ref.umap_", assay = 'RNA')

coembed$Cluster <- NA
coembed$Cluster[match(colnames(tmp_object), colnames(coembed))] <- tmp_object$Cluster
coembed$cell_identity0 <- NA
coembed$cell_identity0[match(colnames(tmp_object), colnames(coembed))] <- as.character(tmp_object$cell_identity)
DimPlot(coembed, group.by = "cell_identity0", cols = cluster_colors) + SetAxes()

qsave(coembed, "output/E14_brain_mouse_scRNA_2021_Nature_coembed.qs")

coembed <- subset(coembed, cells = c(colnames(temp_object), colnames(tmp_object)[grep("Cortex_|Endo|Fibrob|Stromal|MGE_|LGE_IPC", tmp_object$cell_identity)]))
DimPlot(coembed, group.by = "cell_identity0", cols = cluster_colors) + SetAxes()
#!---------------------------------------------------------------------------------------------------------------
# deal E16_5 spatial RNA
tmp_object <- seurat_list$E16_5
tmp_object <- FindClusters(tmp_object, resolution = seq(0.1, 2, 0.1), graph.name = "SCT_snn")
wrap_plots(map(seq(0.6, 2, 0.1), function(x) DimPlot(tmp_object, reduction = "SCT_umap", group.by = paste0("SCT_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)
Seurat::SpatialPlot(tmp_object, shape = 22, cols = selected_colors, group.by = paste0("SCT_snn_res.", seq(1.1, 2, 0.1)), ncol = 5) & NoLegend()

tmp_object$cell_identity <- as.character(seurat_object$cell_identity[match(paste0("E16_5_", colnames(tmp_object)), colnames(seurat_object))])
DimPlot(tmp_object, reduction = "SCT_umap", group.by = "cell_identity", cols = cluster_colors, label = T)
tmp_object <- FindClusters(tmp_object, cluster.name = "SCT_cluster", resolution = 1.8, graph.name = "SCT_snn")

tmp_object$cell_identity <- as.character(tmp_object$SCT_cluster)
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(12))] <- "Cortex_DP_AP"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(27))] <- "Cortex_MP_AP"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(9, 13))] <- "Cortex_IPC"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(3))] <- "Cortex_DP_MigN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(16, 24))] <- "Cortex_MP_MigN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(15))] <- "Cortex_Layer_6b"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(0))] <- "Cortex_DP_CThPN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(7))] <- "Cortex_MP_CThPN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(5))] <- "Cortex_DP_SCPN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(4))] <- "Cortex_MP_SCPN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(6))] <- "Cortex_UL_CPN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(14))] <- "Cortex_MGE_derived_InN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(20))] <- "LGE_AP"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(17))] <- "LGE_IPC"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(2))] <- "LGE_InN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(8))] <- "CPU_pre_MSN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(1, 11))] <- "CPU_MSN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(25))] <- "CPU_Endo"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(26))] <- "ChP_MCC"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(23))] <- "LCS_IMN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(22))] <- "CLA_ExN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(28))] <- "LS_AP"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(18, 19))] <- "LS_IPC"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(10))] <- "Endo"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(21))] <- "Fibrob"

tmp_object$cell_identity <- factor(tmp_object$cell_identity, levels = c("Cortex_DP_AP", "Cortex_MP_AP", "Cortex_IPC", "Cortex_DP_MigN", "Cortex_MP_MigN", "Cortex_Layer_6b", "Cortex_DP_CThPN", "Cortex_MP_CThPN", "Cortex_DP_SCPN", "Cortex_MP_SCPN", "Cortex_UL_CPN", "Cortex_MGE_derived_InN", "LGE_AP", "LGE_IPC", "LGE_InN", "CPU_pre_MSN", "CPU_MSN", "CPU_Endo", "ChP_MCC", "LCS_IMN", "CLA_ExN", "LS_AP", "LS_IPC", "Endo", "Fibrob"))
Seurat::SpatialPlot(tmp_object, group.by = "cell_identity", shape = 22, cols = cluster_colors, crop = F)
#!------------------------------
Idents(tmp_object) <- "cell_identity"
all_markers <- FindAllMarkers(tmp_object, min.pct = 0.1, only.pos = T)
de_markers <- all_markers[which(all_markers$avg_log2FC > 0.5 & all_markers$p_val_adj < 0.01 & all_markers$pct.1 > 0.3),]
write.csv(de_markers, "E16_5_ST_cell_identity_DE_markers.csv", quote = F)

key_markers <- c("Satb2", "Top2a", "Mki67", "Eomes", "Pax6", "Veph1", "Ccdc80", "Aldoc", "Neurog2", "Hmga2", "Dct", "Pter", "Hopx", "Ltbp1", "Slc17a6", "Nrn1", "Adamts2","Calb2", "Neurod1", "Cemip", "Nr3c2", "Abca8a", "Tmem132d", "Pappa2", "Nfe2l3", "Klhl1", "Npy", "Npr3",  "Nr4a3","Igfbp3", "Nwd2", "Tafa1", "Tcerg1l", "Rorb", "Ntf3", "Trmt9b", "Chst8", "Usp43", "Tafa2", "Bcl6", "Inhba", "Frem2", "Mkx", "Mafb", "Lhx6", "Npas1", "Nxph1", "Lhcgr", "Pdzph1", "Vcam1", "Sp9", "Dlx6", "Dlx1", "Dlx2", "Crb1", "Six3", "Drd2", "Dchs2", "Zfp503", "Ikzf1", "Cntnap3", "Dclk3", "Oprm1", "Ano3", "Rxrg", "Myo3b", "Kdr", "Zfp366", "Fli1", "Cldn5", "Tmem72", "Folr1", "Clic6", "Gmnc", "Tfap2d", "Lypd6", "Nr4a2", "Prss12", "Runx2", "Daw1", "Hepacam", "Mlc1", "Glis3", "Zic4", "Zic1", "Zic2", "Epas1", "Igfbp7", "Ndnf", "Col3a1", "Col1a1")

DotPlot(tmp_object, features = key_markers, col.min = 0.1, dot.min = 0.1) + RotatedAxis() + labs(y = "cell identity", x = "Features") + scale_color_gradientn(colours = paletteContinuous(set = "whitePurple"))
ggsave("Development_QC_cell_identity/E16_5_cell_identity_key_markers.pdf")
#!----------------------------
tmp_object$cell_type <- as.character(tmp_object$cell_identity)
tmp_object$cell_type <- gsub(tmp_object$cell_type, pattern = "Cortex_", replacement = "")
tmp_object$cell_type <- gsub(tmp_object$cell_type, pattern = "DP_", replacement = "")
tmp_object$cell_type <- gsub(tmp_object$cell_type, pattern = "MP_", replacement = "")

tmp_object$cell_type <- factor(tmp_object$cell_type, levels = c("AP", "IPC", "MigN", "Layer_6b", "CThPN", "SCPN", "UL_CPN", "MGE_derived_InN", "LGE_AP", "LGE_IPC", "LGE_InN", "CPU_pre_MSN", "CPU_MSN", "CPU_Endo", "ChP_MCC", "LCS_IMN", "CLA_ExN", "LS_AP", "LS_IPC", "Endo", "Fibrob"))
Seurat::SpatialPlot(tmp_object, group.by = "cell_type", shape = 22, cols = cluster_colors, crop = FALSE)

tmp_object$Cluster <- ClusterID[match(as.character(tmp_object$cell_identity), ClusterNames)]
tmp_object$Cluster <- factor(tmp_object$Cluster, levels = intersect(ClusterID, unique(tmp_object$Cluster)))
Seurat::SpatialPlot(tmp_object, group.by = "Cluster", shape = 22, cols = ClusterIDColors, crop = F)

tmp_object <- FindClusters(tmp_object, resolution = seq(0.1, 1, 0.1), graph.name = "banksy_snn")

wrap_plots(map(seq(0.1, 1, 0.1), function(x) DimPlot(tmp_object, reduction = "banksy_umap", group.by = paste0("banksy_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)
Seurat::SpatialPlot(tmp_object, shape = 22, cols = selected_colors, group.by = paste0("banksy_snn_res.", seq(0.1, 1, 0.1)), ncol = 5) & NoLegend()
tmp_object$banksy_cluster <- as.character(tmp_object$banksy_snn_res.0.2)
tmp_object$banksy_cluster[which(tmp_object$banksy_snn_res.0.5 == 14)] <- 11
tmp_object$banksy_cluster <- paste0("D", tmp_object$banksy_cluster)

tmp_object$banksy_domains <- as.character(tmp_object$banksy_cluster)
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D4"))] <- "Cortex_1"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D3"))] <- "Cortex_2"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D0"))] <- "Cortex_3"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D2"))] <- "Cortex_4"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D6"))] <- "Cortex_5"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D9"))] <- "RigN"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D11"))] <- "LGE0"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D5"))] <- "LGE"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D1"))] <- "CPU"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D10"))] <- "ChP"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D8"))] <- "CLA"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D7"))] <- "LS"
tmp_object$banksy_domains <- factor(tmp_object$banksy_domains, levels = c("Cortex_1", "Cortex_2", "Cortex_3", "Cortex_4", "Cortex_5", "MigN", "LGE0", "LGE", "CPU", "ChP", "CLA", "LS"))
tmp_object$banksy_cluster <- domainID[match(as.character(tmp_object$banksy_domains), domainNames)]
tmp_object$banksy_cluster <- factor(tmp_object$banksy_cluster, levels = intersect(domainID, unique(tmp_object$banksy_cluster)))
Seurat::SpatialPlot(tmp_object, group.by = "banksy_cluster", shape = 22, cols = domainIDColors, crop = FALSE)

# cell2location analysis
cell2location <- read.csv("cell2location_Res/development_FC_cell2location/sinlge_timepoint/ST_E16_5_sp_cell2location.csv", row.names = 1)
all(rownames(cell2location) == colnames(tmp_object))
colnames(cell2location)[18:27] <- paste0("E16_5_", colnames(cell2location)[18:27])
tmp_object@meta.data <- cbind(tmp_object@meta.data, cell2location[,18:27])
cell2location <- read.csv("cell2location_Res/development_FC_cell2location/timepoints/ST_E16_5_sp_cell2location.csv", row.names = 1, check.names = F)
all(rownames(cell2location) == colnames(tmp_object))
tmp_object@meta.data <- cbind(tmp_object@meta.data, cell2location[,18:36])
#!-----------------------------------------------------------------
temp_object <- readRDS("output/E16_brain_mouse_scRNA_2021_Nature.rds")

DefaultAssay(tmp_object) <- "RNA"
tmp_object <- NormalizeData(tmp_object)
anchors <- FindTransferAnchors(reference = temp_object, query = tmp_object, dims = 1:30, reference.reduction = "pca", reference.assay= "RNA", query.assay = "RNA", reduction = "pcaproject")
tmp_object <- MapQuery(anchorset = anchors, reference = temp_object, query = tmp_object, refdata = list(celltype = "cell_identity"), reference.reduction = "pca", new.reduction.name = "ref.pca", reduction.model = "umap", projectumap.args = list(reduction.name = "ref.umap"))
colnames(tmp_object@meta.data)[82] <- "predicted_E16_5_celltype_2021_Nature.score"
colnames(tmp_object@meta.data)[83] <- "predicted_E16_5_celltype_2021_Nature"

tmp_object$predicted_E16_5_celltype_2021_Nature <- factor(tmp_object$predicted_E16_5_celltype_2021_Nature, levels = intersect(levels(temp_object$cell_identity), unique(tmp_object$predicted_E16_5_celltype_2021_Nature)))

predictions <- table(tmp_object$Cluster, tmp_object$predicted_E16_5_celltype_2021_Nature)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)
predictions$Var1 <- paste0(predictions$Var1, "_", ClusterNames[match(predictions$Var1, ClusterID)])
predictions$Var1 <- factor(predictions$Var1, levels = intersect(paste0(ClusterID, "_", ClusterNames), unique(predictions$Var1)))

p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
    low = "#ffffc8", high = "#7d0025") + xlab("SCT clusters") + ylab("Predicted cell type label with 2021 Nature E16.5") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
print(p1)
ggsave("predicted_E16_5_celltype_2021_Nature_Clusters_all.pdf")

predictions <- predictions[grep("Cortex_|C30_Endo|C31_Fibrob|C33_Stromal|LGE_AP|LGE_IPC", predictions$Var1),]

p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
    low = "#ffffc8", high = "#7d0025") + xlab("SCT clusters") + ylab("Predicted cell type label with 2021 Nature E16.5") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
print(p1)
ggsave("predicted_E16_5_celltype_2021_Nature_Clusters.pdf")

DefaultAssay(tmp_object) <- "SCT"
Idents(tmp_object) <- "cell_identity"

Seurat::SpatialPlot(tmp_object, shape = 22, cols = ClusterIDColors, group.by = "Cluster", crop = FALSE)
Seurat::SpatialPlot(tmp_object, shape = 22, cols = domainIDColors, group.by = "banksy_cluster", crop = FALSE)
Seurat::SpatialPlot(tmp_object, shape = 22, cols = cluster_colors, crop = FALSE)
DimPlot(tmp_object, cols = ClusterIDColors, group.by = "Cluster", label = T, repel = T) + SetAxes() + NoLegend()
DimPlot(tmp_object, cols = cluster_colors, label = F) + SetAxes()
DimPlot(tmp_object, reduction = "banksy_umap", cols = domainIDColors, group.by = "banksy_cluster", label = T, repel = T) + SetAxes()

library(ggsankey)

meta_data <- tmp_object@meta.data %>% make_long(Cluster, banksy_cluster)
meta_data$node <- factor(meta_data$node, levels = c(levels(tmp_object$Cluster), intersect(domainID, unique(meta_data$node))))
meta_data$next_node <- factor(meta_data$next_node, levels = c(levels(tmp_object$Cluster), intersect(domainID, unique(meta_data$next_node))))

ggplot(meta_data, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) + geom_sankey(color = "black", lwd = 0.2, flow.alpha = 0.5, show.legend = FALSE, na.rm = FALSE) + scale_x_discrete(expand = c(0, 0.8)) + 
		scale_fill_manual(values = c(ClusterIDColors, domainIDColors), drop = FALSE) + geom_sankey_text(size = 3.2, color = "black") + theme_void() + theme(axis.text.x = element_text())
ggsave("Development_QC_cell_identity/E16_5_SCT_spot_cluster_banksy_domian_sankey_plot.pdf")
#!--------------------------------------------------------------------------------
library(UCell)
source("/home/yiyelinfeng/scripts/Rscripts/lung_project/IPF/spatial-lung-fibrosis/scripts/custom_colors.R")

selected_clusters <- paste0("E16_5_", levels(temp_object$cell_identity))

Idents(tmp_object) <- "Cluster"
DefaultAssay(tmp_object) <- "SCT"
all_markers <- FindAllMarkers(tmp_object, only.pos = T, min.pct = 0.1)

de_markers <- all_markers[which(all_markers$avg_log2FC > 0.5 & all_markers$p_val_adj < 0.01 & all_markers$pct.1 > 0.3),]
marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:30]
	count = count + 1
}
names(marker_genes) <- unique(de_markers$cluster)

tmp_object <- UCell::AddModuleScore_UCell(tmp_object, features = marker_genes, ncores = 1, name = "")

cluster_c2l_density <- cor(tmp_object@meta.data[,which(colnames(tmp_object@meta.data) %in% c(levels(tmp_object$Cluster), selected_clusters))])
diag(cluster_c2l_density) <- 0
cluster_c2l_density <- as.data.frame(cluster_c2l_density)
cluster_c2l_density <- cluster_c2l_density[setdiff(rownames(cluster_c2l_density), levels(tmp_object$Cluster)) , levels(tmp_object$Cluster)]
colnames(cluster_c2l_density) <- paste0(colnames(cluster_c2l_density), "_", ClusterNames[match(colnames(cluster_c2l_density), ClusterID)])

pheatmap::pheatmap(cluster_c2l_density, cellwidth = 15, cellheight = 15, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_div_custom2, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_Nat_2021_", names(seurat_list)[2], "_all.pdf"), width = 10, height = 10)

cluster_c2l_density <- cluster_c2l_density[selected_clusters, grep("Cortex_|C33_Endo|C34_Fibrob|LGE_AP|LGE_", colnames(cluster_c2l_density))]

pheatmap::pheatmap(cluster_c2l_density, cellwidth = 15, cellheight = 15, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_div_custom2, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_Nat_2021_", names(seurat_list)[2], ".pdf"), width = 10, height = 10)

tmp_meta <- as.matrix(tmp_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(tmp_object$Cluster, "_", tmp_object$cell_identity)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)
tmp_meta0 <- tmp_meta0[paste0(levels(tmp_object$Cluster), "_", ClusterNames[match(levels(tmp_object$Cluster), names(ClusterNames))]), selected_clusters]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[2], "_zscore_all.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[2], "_zscore_all1.pdf"), width = 10, height = 10)

tmp_meta <- as.matrix(tmp_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(tmp_object$Cluster, "_", tmp_object$cell_identity)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta <- tmp_meta[grep("Cortex_|C33_Endo|C34_Fibrob|LGE_AP|LGE_", rownames(tmp_meta)),]
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)
tmp_meta0 <- tmp_meta0[intersect(paste0(levels(tmp_object$Cluster), "_", ClusterNames[match(levels(tmp_object$Cluster), names(ClusterNames))]), rownames(tmp_meta0)), selected_clusters]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[2], "_zscore.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[2], "_zscore1.pdf"), width = 10, height = 10)

tmp_meta <- as.matrix(tmp_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(tmp_object$Cluster, "_", tmp_object$cell_identity)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta <- tmp_meta[grep("Cortex_", rownames(tmp_meta)),]
tmp_meta <- tmp_meta[-grep("MGE_",rownames(tmp_meta)),selected_clusters[1:7]]
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)
tmp_meta0 <- tmp_meta0[intersect(paste0(levels(tmp_object$Cluster), "_", ClusterNames[match(levels(tmp_object$Cluster), names(ClusterNames))]), rownames(tmp_meta0)),]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[2], "_zscore0.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[2], "_zscore01.pdf"), width = 10, height = 10)

selected_clusters <- c("Apical_progenitors", "Intermediate_progenitors", "Migrating_neurons", "Immature_neurons", "Layer_6b", "NP", "CThPN", "Layer_4", "SCPN", "DL_CPN", "UL_CPN", "Cajal_Retzius_cells", "Interneurons", "Astrocytes", "Oligodendrocytes", "Cycling_glial_cells", "Endothelial_cells", "Pericytes", "Microglia")

tmp_meta <- as.matrix(tmp_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(tmp_object$Cluster, "_", tmp_object$cell_identity)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)
tmp_meta0 <- tmp_meta0[paste0(levels(tmp_object$Cluster), "_", ClusterNames[match(levels(tmp_object$Cluster), names(ClusterNames))]), selected_clusters]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[2], "_zscore_all_times.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[2], "_zscore_all1_times.pdf"), width = 10, height = 10)

tmp_meta <- as.matrix(tmp_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(tmp_object$Cluster, "_", tmp_object$cell_identity)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta <- tmp_meta[grep("Cortex_|C33_Endo|C34_Fibrob|LGE_AP|LGE_", rownames(tmp_meta)),]
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)
tmp_meta0 <- tmp_meta0[intersect(paste0(levels(tmp_object$Cluster), "_", ClusterNames[match(levels(tmp_object$Cluster), names(ClusterNames))]), rownames(tmp_meta0)), selected_clusters]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[2], "_zscore_times.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[2], "_zscore1_times.pdf"), width = 10, height = 10)

Idents(tmp_object) <- "cell_identity"
seurat_list$E16_5 <- tmp_object
#!-------------------------------------------------------------------------------------------------
DimPlot(temp_object, group.by = "cell_identity", cols = ArchRPalettes$kelly, label = T, repel = T) + SetAxes()

tmp_object@assays <- tmp_object@assays[1]
tmp_object@meta.data <- tmp_object@meta.data[,c(1:5, 8:15, 39:41, 81, 82)]

coembed <- merge(x = temp_object, y = tmp_object)
coembed <- JoinLayers(coembed)
coembed@assays$RNA@layers <- coembed@assays$RNA@layers[1:2]

embedding <- tmp_object@reductions$ref.umap@cell.embeddings
embeddings <- rbind(temp_object@reductions$umap@cell.embeddings, embedding)
all(rownames(embeddings) == colnames(coembed))

coembed[["ref.umap"]] <- CreateDimReducObject(embeddings = embeddings, key = "ref.umap_", assay = 'RNA')

coembed$Cluster <- NA
coembed$Cluster[match(colnames(tmp_object), colnames(coembed))] <- tmp_object$Cluster
coembed$cell_identity0 <- NA
coembed$cell_identity0[match(colnames(tmp_object), colnames(coembed))] <- as.character(tmp_object$cell_identity)

DimPlot(coembed, group.by = "cell_identity0", cols = cluster_colors) + SetAxes()

qsave(coembed, "output/E16_brain_mouse_scRNA_2021_Nature_coembed.qs")
coembed <- subset(coembed, cells = c(colnames(temp_object), colnames(tmp_object)[grep("Cortex_|$Endo|$Fibrob|LGE_IPC|LGE_AP", tmp_object$cell_identity)]))
DimPlot(coembed, group.by = "cell_identity0", cols = cluster_colors) + SetAxes()
#!---------------------------------------------------------------------------------------------------------------
# deal E18_5 spatial RNA
tmp_object <- seurat_list$E18_5
tmp_object <- FindClusters(tmp_object, resolution = seq(0.1, 2, 0.1), graph.name = "SCT_snn")
wrap_plots(map(seq(0.6, 2, 0.1), function(x) DimPlot(tmp_object, reduction = "SCT_umap", group.by = paste0("SCT_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)
Seurat::SpatialPlot(tmp_object, shape = 22, cols = selected_colors, group.by = paste0("SCT_snn_res.", seq(1.1, 2, 0.1)), ncol = 5) & NoLegend()

tmp_object$cell_identity <- as.character(seurat_object$cell_identity[match(paste0("E18_5_", colnames(tmp_object)), colnames(seurat_object))])
DimPlot(tmp_object, reduction = "SCT_umap", group.by = "cell_identity", cols = cluster_colors, label = T)
tmp_object <- FindClusters(tmp_object, cluster.name = "SCT_cluster", resolution = 1.3, graph.name = "SCT_snn")

tmp_object$cell_identity <- as.character(tmp_object$SCT_cluster)
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(9))] <- "Cortex_DP_AP"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(29))] <- "Cortex_MP_AP"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(7))] <- "Cortex_IPC"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(1, 17))] <- "Cortex_DP_MigN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(16))] <- "Cortex_MP_MigN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(14))] <- "Cortex_Layer_6b"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(2, 3))] <- "Cortex_DP_CThPN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(6))] <- "Cortex_MP_CThPN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(8))] <- "Cortex_DP_SCPN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(19))] <- "Cortex_MP_SCPN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(5, 12))] <- "Cortex_UL_CPN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(15))] <- "Cortex_MGE_derived_InN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(24))] <- "LGE_AP"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(22))] <- "LGE_IPC"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(10))] <- "LGE_InN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(4))] <- "CPU_pre_MSN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(0, 21))] <- "CPU_MSN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(23))] <- "CPU_Endo"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(25))] <- "ChP_MCC"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(20))] <- "LCS_IMN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(28))] <- "Cortex_DP_CThPN"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(26))] <- "LS_AP"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(11, 18))] <- "LS_IPC"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(13))] <- "Endo"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(27))] <- "MG"

tmp_object$cell_identity <- factor(tmp_object$cell_identity, levels = c("Cortex_DP_AP", "Cortex_MP_AP", "Cortex_IPC", "Cortex_DP_MigN", "Cortex_MP_MigN", "Cortex_Layer_6b", "Cortex_DP_CThPN", "Cortex_MP_CThPN", "Cortex_DP_SCPN", "Cortex_MP_SCPN", "Cortex_UL_CPN", "Cortex_MGE_derived_InN", "LGE_AP", "LGE_IPC", "LGE_InN", "CPU_pre_MSN", "CPU_MSN", "CPU_Endo", "ChP_MCC", "LCS_IMN", "LS_AP", "LS_IPC", "Endo", "MG"))
Seurat::SpatialPlot(tmp_object, group.by = "cell_identity", shape = 22, cols = cluster_colors, crop = F)
#!------------------------------
Idents(tmp_object) <- "cell_identity"
all_markers <- FindAllMarkers(tmp_object, min.pct = 0.1, only.pos = T)
de_markers <- all_markers[which(all_markers$avg_log2FC > 0.5 & all_markers$p_val_adj < 0.01 & all_markers$pct.1 > 0.3),]
write.csv(de_markers, "E18_5_ST_cell_identity_DE_markers.csv", quote = F)

key_markers <- c("Satb2", "Top2a", "Mki67", "Eomes", "Pax6", "Veph1", "Ccdc80", "Aldoc", "Neurog2", "Dct", "Pter", "Hopx", "Ltbp1", "Sema3c", "Slc17a6", "Pou3f2", "Nrn1", "Prox1", "Ndst3", "Lypd6", "Adamts2", "Fst", "Calb2", "Cemip", "Nr3c2", "Abca8a", "Tmem132d", "Pappa2", "Nfe2l3", "Lmo7", "Necab1", "Cdh9", "Sv2b", "Klhl1", "Npy", "Rspo3", "Npr3", "Pde1a", "Nr4a3", "Etv5", "Nwd2", "Tafa1", "Tcerg1l", "Adamts20", "Rorb", "L3mbtl4", "Rasgrf2", "Arhgap28", "Mlip", "Mkx", "Ntf3", "Usp43", "Inhba", "Tafa2", "Bcl6", "Pdzrn4", "Mafb", "Nxph2",  "Npas1", "Nxph1", "Lhx6", "Dgkg", "Lhcgr", "Vcam1", "Pdzph1", "Hepacam", "Mlc1", "Dlx2", "Sp8", "Akna", "Smoc1", "St18", "Cdca7", "Sp9", "Dlx6", "Dlx1", "Mob3b", "Six3", "Drd2", "Dchs2", "Ikzf1", "Cntnap3", "Dclk3", "Ano3", "Rxrg", "Oprm1", "Myo3b", "Tgfa", "Ngef", "Gpr88", "Cldn5", "Tmem72", "Folr1", "Clic6", "Gmnc", "Tfap2d", "Pdgfra", "Cspg4", "Pstpip2", "Prkcq", "Shc4", "Dynlrb2", "Daw1", "Crocc2", "Deup1", "Wdr63", "Zic4", "Zic1", "Snhg11", "Fgd5", "Flt1", "Eng", "Epas1", "Igfbp7", "Adgrl4", "C1qb", "Cd180", "Dock2", "Spp1")

DotPlot(tmp_object, features = key_markers, col.min = 0.1, dot.min = 0.1) + RotatedAxis() + labs(y = "cell identity", x = "Features") + scale_color_gradientn(colours = paletteContinuous(set = "whitePurple"))
ggsave("Development_QC_cell_identity/E18_5_cell_identity_key_markers.pdf")
#!----------------------------
tmp_object$cell_type <- as.character(tmp_object$cell_identity)
tmp_object$cell_type <- gsub(tmp_object$cell_type, pattern = "Cortex_", replacement = "")
tmp_object$cell_type <- gsub(tmp_object$cell_type, pattern = "DP_", replacement = "")
tmp_object$cell_type <- gsub(tmp_object$cell_type, pattern = "MP_", replacement = "")

tmp_object$cell_type <- factor(tmp_object$cell_type, levels = c("AP", "IPC", "MigN", "Layer_6b", "CThPN", "SCPN", "UL_CPN", "MGE_derived_InN", "LGE_AP", "LGE_IPC", "LGE_InN", "CPU_pre_MSN", "CPU_MSN", "CPU_Endo", "ChP_MCC", "LCS_IMN", "LS_AP", "LS_IPC", "Endo", "MG"))
Seurat::SpatialPlot(tmp_object, group.by = "cell_type", shape = 22, cols = cluster_colors)

tmp_object$Cluster <- ClusterID[match(as.character(tmp_object$cell_identity), ClusterNames)]
tmp_object$Cluster <- factor(tmp_object$Cluster, levels = intersect(ClusterID, unique(tmp_object$Cluster)))
Seurat::SpatialPlot(tmp_object, group.by = "Cluster", shape = 22, cols = ClusterIDColors, crop = F)

tmp_object <- FindClusters(tmp_object, resolution = seq(0.1, 1, 0.1), graph.name = "banksy_snn")

wrap_plots(map(seq(0.1, 1, 0.1), function(x) DimPlot(tmp_object, reduction = "banksy_umap", group.by = paste0("banksy_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)
Seurat::SpatialPlot(tmp_object, shape = 22, cols = selected_colors, group.by = paste0("banksy_snn_res.", seq(0.1, 1, 0.1)), ncol = 5) & NoLegend()
tmp_object$banksy_cluster <- as.character(tmp_object$banksy_snn_res.0.2)
tmp_object$banksy_cluster[which(tmp_object$banksy_snn_res.0.6 == 16)] <- 10
tmp_object$banksy_cluster[which(tmp_object$banksy_snn_res.0.6 == 14)] <- 11
tmp_object$banksy_cluster[which(tmp_object$banksy_snn_res.0.6 == 9)] <- 12
tmp_object$banksy_cluster <- paste0("D", tmp_object$banksy_cluster)

tmp_object$banksy_domains <- as.character(tmp_object$banksy_cluster)
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D4"))] <- "Cortex_1"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D3"))] <- "Cortex_2"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D0", "D7"))] <- "Cortex_3"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D2"))] <- "Cortex_4"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D12"))] <- "MigN"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D11"))] <- "LGE0"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D6"))] <- "LGE"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D1"))] <- "CPU"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D10"))] <- "ChP"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D5"))] <- "LS"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c("D9"))] <- "End"
tmp_object$banksy_domains <- factor(tmp_object$banksy_domains, levels = c("Cortex_1", "Cortex_2", "Cortex_3", "Cortex_4", "MigN", "LGE0", "LGE", "CPU", "ChP", "LS", "End"))
tmp_object$banksy_cluster <- domainID[match(as.character(tmp_object$banksy_domains), domainNames)]
tmp_object$banksy_cluster <- factor(tmp_object$banksy_cluster, levels = intersect(domainID, unique(tmp_object$banksy_cluster)))
Seurat::SpatialPlot(tmp_object, group.by = "banksy_cluster", shape = 22, cols = domainIDColors, crop = F)

# cell2location analysis
cell2location <- read.csv("cell2location_Res/development_FC_cell2location/sinlge_timepoint/ST_E18_5_sp_cell2location.csv", row.names = 1)
all(rownames(cell2location) == colnames(tmp_object))
colnames(cell2location)[18:34] <- paste0("E18_5_", colnames(cell2location)[18:34])
tmp_object@meta.data <- cbind(tmp_object@meta.data, cell2location[,18:34])
cell2location <- read.csv("cell2location_Res/development_FC_cell2location/timepoints/ST_E18_5_sp_cell2location.csv", row.names = 1, check.names = F)
all(rownames(cell2location) == colnames(tmp_object))
tmp_object@meta.data <- cbind(tmp_object@meta.data, cell2location[,18:36])
#!-----------------------------------------------------------------
temp_object <- readRDS("output/E18_brain_mouse_scRNA_2021_Nature.rds")

DefaultAssay(tmp_object) <- "RNA"
tmp_object <- NormalizeData(tmp_object)
anchors <- FindTransferAnchors(reference = temp_object, query = tmp_object, dims = 1:30, reference.reduction = "cca", reference.assay= "RNA", query.assay = "RNA", reduction = "pcaproject")
tmp_object <- MapQuery(anchorset = anchors, reference = temp_object, query = tmp_object, refdata = list(celltype = "cell_identity"), reference.reduction = "cca", new.reduction.name = "ref.pca", reduction.model = "cca_umap", projectumap.args = list(reduction.name = "ref.umap"))
colnames(tmp_object@meta.data)[89] <- "predicted_E18_5_celltype_2021_Nature.score"
colnames(tmp_object@meta.data)[90] <- "predicted_E18_5_celltype_2021_Nature"
tmp_object$predicted_E18_5_celltype_2021_Nature <- factor(tmp_object$predicted_E18_5_celltype_2021_Nature, levels = intersect(levels(temp_object$cell_identity), unique(tmp_object$predicted_E18_5_celltype_2021_Nature)))

predictions <- table(tmp_object$Cluster, tmp_object$predicted_E18_5_celltype_2021_Nature)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)
predictions$Var1 <- paste0(predictions$Var1, "_", ClusterNames[match(predictions$Var1, ClusterID)])
predictions$Var1 <- factor(predictions$Var1, levels = intersect(paste0(ClusterID, "_", ClusterNames), unique(predictions$Var1)))

p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
    low = "#ffffc8", high = "#7d0025") + xlab("SCT clusters") + ylab("Predicted cell type label with 2021 Nature E18.5") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
print(p1)
ggsave("predicted_E18_5_celltype_2021_Nature_Clusters_all.pdf")

predictions <- predictions[grep("Cortex_|C30_Endo|C31_Fibrob|C33_Stromal|LGE_AP|LGE_IPC", predictions$Var1),]

p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
    low = "#ffffc8", high = "#7d0025") + xlab("SCT clusters") + ylab("Predicted cell type label with 2021 Nature E18.5") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
print(p1)
ggsave("predicted_E18_5_celltype_2021_Nature_Clusters.pdf")

DefaultAssay(tmp_object) <- "SCT"
Idents(tmp_object) <- "cell_identity"

Seurat::SpatialPlot(tmp_object, shape = 22, cols = ClusterIDColors, group.by = "Cluster", crop = F)
Seurat::SpatialPlot(tmp_object, shape = 22, cols = domainIDColors, group.by = "banksy_cluster", crop = F)
Seurat::SpatialPlot(tmp_object, shape = 22, cols = cluster_colors, crop = F)
DimPlot(tmp_object, cols = ClusterIDColors, group.by = "Cluster", label = T, repel = T) + SetAxes() + NoLegend()
DimPlot(tmp_object, cols = cluster_colors) + SetAxes()
DimPlot(tmp_object, reduction = "banksy_umap", cols = domainIDColors, group.by = "banksy_cluster", label = T, repel = T) + SetAxes()

library(ggsankey)

meta_data <- tmp_object@meta.data %>% make_long(Cluster, banksy_cluster)
meta_data$node <- factor(meta_data$node, levels = c(levels(tmp_object$Cluster), intersect(domainID, unique(meta_data$node))))
meta_data$next_node <- factor(meta_data$next_node, levels = c(levels(tmp_object$Cluster), intersect(domainID, unique(meta_data$next_node))))

ggplot(meta_data, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) + geom_sankey(color = "black", lwd = 0.2, flow.alpha = 0.5, show.legend = FALSE, na.rm = FALSE) + scale_x_discrete(expand = c(0, 0.8)) + 
		scale_fill_manual(values = c(ClusterIDColors, domainIDColors), drop = FALSE) + geom_sankey_text(size = 3.2, color = "black") + theme_void() + theme(axis.text.x = element_text())
ggsave("Development_QC_cell_identity/E18_5_SCT_spot_cluster_banksy_domian_sankey_plot.pdf")
#!--------------------------------------------------------------------------------
library(UCell)
source("/home/yiyelinfeng/scripts/Rscripts/lung_project/IPF/spatial-lung-fibrosis/scripts/custom_colors.R")

selected_clusters <- paste0("E18_5_", levels(temp_object$cell_identity))

Idents(tmp_object) <- "Cluster"
DefaultAssay(tmp_object) <- "SCT"
all_markers <- FindAllMarkers(tmp_object, only.pos = T, min.pct = 0.1)

de_markers <- all_markers[which(all_markers$avg_log2FC > 0.5 & all_markers$p_val_adj < 0.01 & all_markers$pct.1 > 0.3),]
marker_genes <- list()
count = 1
for(i in unique(de_markers$cluster)){
	marker_genes[[count]] <- de_markers$gene[which(de_markers$cluster == i)][1:30]
	count = count + 1
}
names(marker_genes) <- unique(de_markers$cluster)

tmp_object <- UCell::AddModuleScore_UCell(tmp_object, features = marker_genes, ncores = 1, name = "")

cluster_c2l_density <- cor(tmp_object@meta.data[,which(colnames(tmp_object@meta.data) %in% c(levels(tmp_object$Cluster), selected_clusters))])
diag(cluster_c2l_density) <- 0
cluster_c2l_density <- as.data.frame(cluster_c2l_density)
cluster_c2l_density <- cluster_c2l_density[setdiff(rownames(cluster_c2l_density), levels(tmp_object$Cluster)) , levels(tmp_object$Cluster)]
colnames(cluster_c2l_density) <- paste0(colnames(cluster_c2l_density), "_", ClusterNames[match(colnames(cluster_c2l_density), ClusterID)])
cluster_c2l_density <- cluster_c2l_density[selected_clusters,]

pheatmap::pheatmap(cluster_c2l_density, cellwidth = 15, cellheight = 15, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_div_custom2, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_Nat_2021_", names(seurat_list)[3], "_all.pdf"), width = 10, height = 10)

cluster_c2l_density <- cluster_c2l_density[selected_clusters, grep("Cortex_|C33_Endo|LGE_LGE_|C35_MG", colnames(cluster_c2l_density))]

pheatmap::pheatmap(cluster_c2l_density, cellwidth = 15, cellheight = 15, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_div_custom2, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_Nat_2021_", names(seurat_list)[3], ".pdf"), width = 10, height = 10)

tmp_meta <- as.matrix(tmp_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(tmp_object$Cluster, "_", tmp_object$cell_identity)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)
tmp_meta0 <- tmp_meta0[paste0(levels(tmp_object$Cluster), "_", ClusterNames[match(levels(tmp_object$Cluster), names(ClusterNames))]), selected_clusters]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[3], "_zscore_all.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[3], "_zscore_all1.pdf"), width = 10, height = 10)

tmp_meta <- as.matrix(tmp_object@meta.data[,selected_clusters[1:9]])
rownames(tmp_meta) <- paste0(tmp_object$Cluster, "_", tmp_object$cell_identity)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta <- tmp_meta[grep("Cortex_", rownames(tmp_meta)),]
tmp_meta <- tmp_meta[-grep("MGE_", rownames(tmp_meta)),]
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)
tmp_meta0 <- tmp_meta0[intersect(paste0(levels(tmp_object$Cluster), "_", ClusterNames[match(levels(tmp_object$Cluster), names(ClusterNames))]), rownames(tmp_meta0)),]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[3], "_zscore0.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[3], "_zscore01.pdf"), width = 10, height = 10)

selected_clusters <- c("Apical_progenitors", "Intermediate_progenitors", "Migrating_neurons", "Immature_neurons", "Layer_6b", "NP", "CThPN", "Layer_4", "SCPN", "DL_CPN", "UL_CPN", "Cajal_Retzius_cells", "Interneurons", "Astrocytes", "Oligodendrocytes", "Cycling_glial_cells", "Endothelial_cells", "Pericytes", "Microglia")

tmp_meta <- as.matrix(tmp_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(tmp_object$Cluster, "_", tmp_object$cell_identity)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)
tmp_meta0 <- tmp_meta0[paste0(levels(tmp_object$Cluster), "_", ClusterNames[match(levels(tmp_object$Cluster), names(ClusterNames))]),]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[3], "_zscore_all_times.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[3], "_zscore_all1_times.pdf"), width = 10, height = 10)

tmp_meta <- as.matrix(tmp_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(tmp_object$Cluster, "_", tmp_object$cell_identity)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta <- tmp_meta[grep("Cortex_|C33_Endo|LGE_LGE_|C35_MG", rownames(tmp_meta)),]
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta0))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)
tmp_meta0 <- tmp_meta0[intersect(paste0(levels(tmp_object$Cluster), "_", ClusterNames[match(levels(tmp_object$Cluster), names(ClusterNames))]), rownames(tmp_meta0)),]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[3], "_zscore_times.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_", names(seurat_list)[3], "_zscore1_times.pdf"), width = 10, height = 10)

Idents(tmp_object) <- "cell_identity"
seurat_list$E18_5 <- tmp_object
#!-----------------------------------------------------------------------------------
DimPlot(temp_object, reduction = "cca_umap", group.by = "cell_identity", cols = ArchRPalettes$kelly, label = T, repel = T) + SetAxes()

tmp_object@assays <- tmp_object@assays[1]
tmp_object@meta.data <- tmp_object@meta.data[,c(1:5, 8:15, 18:19, 40:42, 89, 90)]

coembed <- merge(x = temp_object, y = tmp_object)
coembed <- JoinLayers(coembed)
coembed@assays$RNA@layers <- coembed@assays$RNA@layers[1:2]

embedding <- tmp_object@reductions$ref.umap@cell.embeddings
embeddings <- rbind(temp_object@reductions$cca_umap@cell.embeddings, embedding)
all(rownames(embeddings) == colnames(coembed))

coembed[["ref.umap"]] <- CreateDimReducObject(embeddings = embeddings, key = "ref.umap_", assay = 'RNA')

coembed$Cluster <- NA
coembed$Cluster[match(colnames(tmp_object), colnames(coembed))] <- tmp_object$Cluster
coembed$cell_identity0 <- NA
coembed$cell_identity0[match(colnames(tmp_object), colnames(coembed))] <- as.character(tmp_object$cell_identity)

DimPlot(coembed, group.by = "cell_identity0", cols = cluster_colors) + SetAxes()

qsave(coembed, "output/E18_brain_mouse_scRNA_2021_Nature_coembed.qs")

coembed <- subset(coembed, cells = c(colnames(temp_object), colnames(tmp_object)[grep("Cortex_|$Endo|$Fibrob|LGE_IPC|LGE_AP|MG", tmp_object$cell_identity)]))
DimPlot(coembed, group.by = "cell_identity0", cols = cluster_colors) + SetAxes()
#!-----------------------------------------------------------------------
seurat_object <- merge(seurat_list[[1]], seurat_list[2:length(seurat_list)], add.cell.ids = names(seurat_list))

DefaultAssay(seurat_object) <- "spARC_SCT"
seurat_object <- JoinLayers(seurat_object)

DefaultAssay(seurat_object) <- "RNA"
seurat_object$orig.ident <- factor(seurat_object$orig.ident, levels = names(seurat_list))
Idents(seurat_object) <- "orig.ident"

VlnPlot(seurat_object, features = c('nCount_RNA', 'nFeature_RNA', "percent.mt", "percent.ribo"), pt.size = 0, ncol = 2) + NoLegend()
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
seurat_object <- FindClusters(seurat_object, resolution = seq(0.1, 1.5, 0.1), graph.name = 'harmony_snn')
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

resolution <- seq(0.1, 1, 0.1)
DefaultAssay(seurat_object) <- "SCT"
seurat_object <- RunBanksy(seurat_object, lambda = 0.8, assay = 'SCT', slot = 'data', dimx = "sdimx", dimy = "sdimy", features = 'variable', group = 'orig.ident', split.scale = TRUE, k_geom = 15)
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
wrap_plots(map(seq(0.1, 1, 0.1), function(x) DimPlot(seurat_object, reduction = "umap_banksy_harmony", group.by = paste0("banksy_harmony_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)
#!------------------------------------------------------------------
cell2location <- read.csv("cell2location_Res/development_FC_cell2location/ST_Adult7_sp_cell2location.csv", row.names = 1)
all(rownames(cell2location) == colnames(seurat_list[[4]]))
colnames(cell2location)[18:40] <- paste0("D_", colnames(cell2location)[18:40])
seurat_list[[4]]@meta.data <- cbind(seurat_list[[4]]@meta.data, cell2location[,18:40])
cell2location <- read.csv("cell2location_Res/adult_FC_cell2location/ST_Adult7_sp_adult_FC_cell2location.csv", row.names = 1)
all(rownames(cell2location) == colnames(seurat_list[[4]]))
seurat_list[[4]]@meta.data <- cbind(seurat_list[[4]]@meta.data, cell2location[,18:37])

cell2location <- read.csv("cell2location_Res/adult_cell2location/ST_Adult7_sp_cell2location.csv", row.names = 1, check.names = F)
all(rownames(cell2location) == colnames(seurat_list[[4]]))
seurat_list[[4]]@meta.data <- cbind(seurat_list[[4]]@meta.data, cell2location[,18:138])

cell2location <- read.csv("cell2location_Res/development_FC_cell2location/ST_Adult15_sp_cell2location.csv", row.names = 1)
all(rownames(cell2location) == colnames(seurat_list[[5]]))
colnames(cell2location)[18:40] <- paste0("D_", colnames(cell2location)[18:40])
seurat_list[[5]]@meta.data <- cbind(seurat_list[[5]]@meta.data, cell2location[,18:40])

cell2location <- read.csv("cell2location_Res/adult_FC_cell2location/ST_Adult15_sp_adult_FC_cell2location.csv", row.names = 1)
all(rownames(cell2location) == colnames(seurat_list[[5]]))
seurat_list[[5]]@meta.data <- cbind(seurat_list[[5]]@meta.data, cell2location[,18:37])

cell2location <- read.csv("cell2location_Res/adult_cell2location/ST_Adult15_sp_cell2location.csv", row.names = 1, check.names = F)
all(rownames(cell2location) == colnames(seurat_list[[5]]))
seurat_list[[5]]@meta.data <- cbind(seurat_list[[5]]@meta.data, cell2location[,18:138])
#!----------------------------------------------------------------------------------------------------------------------------------------------------------
tmp_object <- readRDS("output/brain_mouse_scRNA_2021_Nature.rds")
tmp_object <- ScaleData(tmp_object)
tmp_object <- RunPCA(tmp_object, npcs = 100)
tmp_object <- RunUMAP(tmp_object, dims = 1:50, return.model = T, reduction.name = "umap1")

tmp_object@reductions$umap1@cell.embeddings[,1] <- tmp_object@reductions$umap@cell.embeddings[,1]
tmp_object@reductions$umap1@cell.embeddings[,2] <- tmp_object@reductions$umap@cell.embeddings[,2]

tmp_object@reductions$umap1@misc$model$embedding[,1] <- tmp_object@reductions$umap@cell.embeddings[,1]
tmp_object@reductions$umap1@misc$model$embedding[,2] <- tmp_object@reductions$umap@cell.embeddings[,2]

tmp_object@reductions <- tmp_object@reductions[2:3]
names(tmp_object@reductions)[2] <- "umap"

temp_object <- readRDS("output/FC_mouse_snRNA_2023_NSMB.rds")
#!-------------------------------------------------------------------------------------------------------
DefaultAssay(seurat_object) <- "RNA"
seurat_object <- NormalizeData(seurat_object)
anchors <- FindTransferAnchors(reference = tmp_object, query = seurat_object, dims = 1:30, reference.reduction = "pca", reference.assay= "RNA", query.assay = "RNA", reduction = "pcaproject")

seurat_object <- MapQuery(anchorset = anchors, reference = tmp_object, query = seurat_object, refdata = list(celltype = "cell_identity"), reference.reduction = "pca", new.reduction.name = "ref.pca", reduction.model = "umap", projectumap.args = list(reduction.name = "ref.umap"))
colnames(seurat_object@meta.data)[178] <- "predicted.celltype_2021_Nature.score"
colnames(seurat_object@meta.data)[179] <- "predicted.celltype_2021_Nature"

anchors <- FindTransferAnchors(reference = temp_object, query = seurat_object, dims = 1:30, reference.reduction = "harmony", reference.assay= "RNA", query.assay = "RNA", reduction = "pcaproject")

seurat_object <- MapQuery(anchorset = anchors, reference = temp_object, query = seurat_object, refdata = list(celltype = "cell_identity"), reference.reduction = "harmony", new.reduction.name = "ref.harmony", reduction.model = "harmony_umap", projectumap.args = list(reduction.name = "ref.harmony.umap"))
colnames(seurat_object@meta.data)[178] <- "predicted.celltype_2023_NSMB.score"
colnames(seurat_object@meta.data)[179] <- "predicted.celltype_2023_NSMB"

seurat_list[[1]] <- seurat_object
#!--------------------------------------------------------------------------------
seurat_object <- readRDS("output/brain_development_ST_final_tutorial.rds")
DefaultAssay(seurat_object) <- "RNA"
seurat_object@assays <- seurat_object@assays[1]
seurat_object@meta.data <- seurat_object@meta.data[,1:5]
seurat_object@reductions <- list()
seurat_object@images <- list()

coembed <- merge(x = seurat_object, y = tmp_object)

embeddings <- c()
for(i in 1:3){
	embedding <- seurat_list[[i]]@reductions$ref.umap@cell.embeddings
	rownames(embedding) <- paste0(names(seurat_list)[i], "_", rownames(embedding))
	embeddings <- rbind(embeddings, embedding)
}

all(rownames(embeddings) == colnames(seurat_object))

embeddings <- rbind(embeddings, tmp_object@reductions$umap@cell.embeddings)

all(rownames(embeddings) == colnames(coembed))

coembed[["ref.umap"]] <- CreateDimReducObject(embeddings = embeddings, key = "ref.umap_", assay = 'RNA')
coembed$cell_identity <- factor(coembed$cell_identity, levels = levels(tmp_object$cell_identity))

meta_data <- c()
for(i in 1:3){
	tmp_meta <- seurat_list[[i]]@meta.data[,c("SCT_cluster", "cluster_cellcharter", "predicted.celltype_2021_Nature.score", "predicted.celltype_2021_Nature")]
	rownames(tmp_meta) <- paste0(names(seurat_list)[i], "_", rownames(tmp_meta))
	meta_data <- rbind(meta_data, tmp_meta)
}
all(rownames(meta_data) == colnames(seurat_object))
coembed <- AddMetaData(coembed, metadata = meta_data)

seurat_object <- readRDS("output/brain_development_ST_final_tutorial.rds")
DefaultAssay(seurat_object) <- "RNA"
seurat_object@assays <- seurat_object@assays[1]
seurat_object@meta.data <- seurat_object@meta.data[,1:5]
seurat_object@reductions <- list()
seurat_object@images <- list()

coembed_adult <- merge(x = seurat_object, y = temp_object)

embeddings <- c()
for(i in 4:5){
	embedding <- seurat_list[[i]]@reductions$ref.harmony.umap@cell.embeddings
	rownames(embedding) <- paste0(names(seurat_list)[i], "_", rownames(embedding))
	embeddings <- rbind(embeddings, embedding)
}

all(rownames(embeddings) == colnames(seurat_object))

embeddings <- rbind(embeddings, temp_object@reductions$harmony_umap@cell.embeddings)

all(rownames(embeddings) == colnames(coembed_adult))

coembed_adult[["ref.umap"]] <- CreateDimReducObject(embeddings = embeddings, key = "ref.umap_", assay = 'RNA')
coembed_adult$cell_identity <- factor(coembed_adult$cell_identity, levels = levels(coembed_adult$cell_identity))

meta_data <- c()
for(i in 4:5){
	tmp_meta <- seurat_list[[i]]@meta.data[,c("SCT_cluster", "cluster_cellcharter", "predicted.celltype_2023_NSMB.score", "predicted.celltype_2023_NSMB")]
	rownames(tmp_meta) <- paste0(names(seurat_list)[i], "_", rownames(tmp_meta))
	meta_data <- rbind(meta_data, tmp_meta)
}
all(rownames(meta_data) == colnames(seurat_object))
coembed_adult <- AddMetaData(coembed_adult, metadata = meta_data)
#!----------------------------------------------------------------------------------------------------------------------------------------
tmp_colors <- selected_colors
names(tmp_colors) <- paste0("C", 0:(length(tmp_colors) - 1))

coembed_adult$samples <- coembed_adult$orig.ident
coembed_adult$samples[which(!coembed_adult$samples %in% c("Adult7", "Adult15"))] <- "NSMB"

selected_cells <- list()
for(i in unique(coembed_adult$samples)){
	tmp_selected_cells <- list()
	for(j in levels(coembed_adult$SCT_cluster)){
		temp <- which(coembed_adult$samples == i & coembed_adult$SCT_cluster == j)
		if(length(temp) > 0){
			tmp_selected_cells[[j]] <- colnames(coembed_adult)[temp]
		}
	}
	selected_cells[[i]] <- tmp_selected_cells
}

tmp_coembed <- subset(coembed_adult, subset = samples %in% c("Adult7", "NSMB"))
DimPlot(tmp_coembed, reduction = "ref.umap", cells.highlight = selected_cells[["Adult7"]], cols.highlight = rev(tmp_colors[match(names(selected_cells[["Adult7"]])[order(names(selected_cells[["Adult7"]]))], names(tmp_colors))]), sizes.highlight = 0.2, raster = FALSE) + SetAxes()
ggsave("Adult7_ref_Adult_NSMB.pdf")

predictions <- table(seurat_list$Adult7$SCT_cluster, seurat_list$Adult7$predicted.celltype_2023_NSMB)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)
p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
    low = "#ffffc8", high = "#7d0025") + xlab("SCT clusters") + ylab("Predicted cell type label with 2023 NSMB") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
print(p1)

tmp_coembed <- subset(coembed_adult, subset = samples %in% c("Adult15", "NSMB"))
DimPlot(tmp_coembed, reduction = "ref.umap", cells.highlight = selected_cells[["Adult15"]], cols.highlight = rev(tmp_colors[match(names(selected_cells[["Adult15"]])[order(names(selected_cells[["Adult15"]]))], names(tmp_colors))]), sizes.highlight = 0.2, raster = FALSE) + SetAxes()
ggsave("Adult15_ref_Adult_NSMB.pdf")

predictions <- table(seurat_list$Adult15$SCT_cluster, seurat_list$Adult15$predicted.celltype_2023_NSMB)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)
p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
    low = "#ffffc8", high = "#7d0025") + xlab("SCT clusters") + ylab("Predicted cell type label with 2023 NSMB") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
print(p1)
#!-----------------------------------------------------------------
# development merge annotation
wrap_plots(map(seq(0.1, 1, 0.1), function(x) DimPlot(seurat_object, reduction = "umap_harmony", group.by = paste0("harmony_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)

Idents(seurat_object) <- "harmony_snn_res.1.5"
dir.create("SCT_harmony_merge_Res")
selected_cells <- CellsByIdentities(seurat_object)
for(i in 1:length(selected_cells)){
	p0 <- Seurat::SpatialPlot(seurat_object, cells.highlight = selected_cells[i], cols.highlight = c("#FFFF00", alpha("lightgrey", 0)), crop = FALSE, shape = 22, combine = T, ncol = 3, alpha = NULL) & NoLegend()
	print(p0)
	ggsave(paste0("SCT_harmony_merge_Res/", names(selected_cells)[i], ".png"))
}

Seurat::SpatialPlot(seurat_object, features = "Sftpc", crop = FALSE, shape = 22, alpha = c(0.3, 1)) + scale_fill_gradientn(colours = paletteContinuous(set = "sambaNight"))
Seurat::SpatialPlot(seurat_object, features = i, crop = FALSE, shape = 22, alpha = c(0.3, 1)) + scale_fill_gradientn(colours = paletteContinuous(set = "sambaNight")) + ggtitle("log(UMI)")

seurat_object$seurat_clusters <- seurat_object$harmony_snn_res.1
seurat_object <- FindSubCluster(seurat_object, graph.name = "harmony_snn", resolution = 0.2, cluster = 6, subcluster.name = "sub_clusters")
Idents(seurat_object) <- "seurat_clusters"

seurat_object$cell_identity <- as.character(seurat_object$sub_clusters)
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c("6_0"))] <- "Cortex_DP_AP" # cortex_DP_apical_progenitors radial_glia "Pax6", "Aldoc", "Fabp7", "Nes", "Vim", "Tnc", "Hmga2", "Lrp4", "Gli2", "Gli3", "Creb5"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c("6_1", "6_2"))] <- "Cortex_MP_AP" # cortex_MP_apical_progenitors "Dct", "Eya1", "Adamts19"

seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(24, 21, 2))] <- "Cortex_IPC" # cortex_intermediate_progenitor_cells "Eomes", "Nrn1", "Slc17a6", "Sstr2", "Sema3c", "Nrn1", "Slc17a6", "Ust", "Lhx2", "Plcb1", "Scrt2", "Sorl1", "Bcar1"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(4))] <- "Cortex_DP_MigN" # cortex_migrating_neurons_cells "Adamts2", "Sema3c", "Pou3f2", "Unc5d", "Nrp1" # 38871984
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(19))] <- "Cortex_MP_MigN" # cortex_migrating_neurons_cells "Fgfr1", "Calb2", "Zfp385b", "Cemip", "Nr3c2", "Camk4", "Wscd2", "Zbtb16", "Dpf3"

seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(15))] <- "Cortex_Layer_6b" # cortex_deep_layer_excitatory_neurons_cells "Tmem132d", "Pappa2", "Cdh18", "Sv2b", "Kcnab1", "Abca8a", "Nxph3", "Ptger3", "Gas7", "Col19a1", "Osbpl10", "Pcsk5", "Met", "Crtac1", "Pou6f2"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(1))] <- "Cortex_DP_CThPN" # cortex_DP_deep_layer_excitatory_neurons_cells "Hs3st4", "Nfe2l3", "Rimbp2", "Rgs6", "Necab1", "Lmo7", "Cnih3", "Spock1", "Wnt7b", "Adamts3", "Cdh6", "Klhl1", "Npas2", "Cdh9", "Mpp6", "Npy", "Pik3r1", "Kctd12", "Foxp2", "Tbr1"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(10))] <- "Cortex_MP_CThPN" # cortex_MP_deep_layer_excitatory_neurons_cells CThPN "Prickle1", "Nr4a3", "Pde1a", "Tspan18", "Sh3gl2", "Adcy1", "B3gat1", "Prickle1", "Zfpm2", "Tle4", "Neto2", "Slc24a2", "Rspo3", "Kitl", "Adamts1", "Fezf2", "Galnt9", "Kcnk1"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(13))] <- "Cortex_DP_SCPN" # Cortex_DP_migrating_excitatory_neurons_cells "Tafa1", "Nwd2", "Cntn6", "Tcerg1l", "Pex5l", "Hcn1", "Cntnap5a", "Lrp1b", "Lmo3", "Ncam2", "Adamts20", "Syt4", "Fat3", "Rorb", "Sorcs1"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(12))] <- "Cortex_MP_SCPN" # cortex_MP_migrating_excitatory_neurons_cells SCPN/CTPN "Mical2", "L3mbtl4", "Rasgrf2", "Trmt9b", "Ntf3", "Usp43", "Msra", "Ldb2", "Eps8", "Fbn2", "Sla", "Arhgap28", "Id2", "Chst8"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(7))] <- "Cortex_UL_CPN" # cortex_upper_layer_excitatory_neurons_cells CPN "Ccbe1", "Tafa2", "Bcl6", "Cdh12", "Tmtc1", "Jakmip1", "Slc16a7", "Inhba", "Mlip", "St6gal2", "Chrna7", "Kif26b", "Dok5", "Eml1", "Ptpro", "Zbtb18", "Limch1", "Itpr1", "Aff2", "Tgfbr1", "Trim67"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(11))] <- "Cortex_CR" # cortex_Cajal-Retzius "Mkx", "Dlgap2", "Reln", "Mafb", "Dgkg", "Nxph2", "Kcnc2", "Nxph1", "Elfn1", "Zbtb16", "Ptchd4", "Cbfa2t3", "Cacng3", "Iqsec1"

seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(23))] <- "MGE_InN" # MGE_Interneuron "Lhx6", "Lhx8", "Nkx2-1", "Adamts5", "Rai2"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(17))] <- "LGE_AP" # LGE_astrocytes "Pdzph1", "Daam2", "Aqp4", "Hes5", "Aldh1l1", "Ascl1"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(14, 26))] <- "LGE_IPC" # LGE_Interneuron "Dlx6", "Dlx1", "Dlx2", "Cdca7"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(8))] <- "LGE_MSN_prog" # LGE_MSN_progenitor "Sp9", "Sp8", "Drd2", "Slit3", "Tac1", "Six3"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(5))] <- "CPU_pre_MSN" # "Drd2", "Gpc6", "Sox8", "Nebl", "Map3k1", "Ikzf1", "Oprm1", "Rasgef1b", "Tmeff2", "Galntl6"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(0))] <- "CPU_MSN" # "Drd1", "Tac1", "Ppp1r1b", "Ngef", "Gnal", "Pde10a", "Foxp1" # 34727523
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(25))] <- "CPU_Endo"

seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(27))] <- "ChP_MCC" # multiciliated cells: "Gmnc", "Lmx1a", "Otx2", "Ins2", "Slc12a2", "Slc4a5", "Aqp1", "Vat1l"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(18))] <- "LCS_IMN" # LCS_immature_neuron "Tfap2d", "Lhfpl3", "Nr2f1"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(22))] <- "CLA_ExN" # CLA_excitatory_neurons Claustrum "Nr4a2", "Kcnip4", "Sv2b", "Il1rapl2", "Galnt14", "Prss12", "Ppp1r14c", "Lypd6b", "Unc13c", "Mtus2", "Grin3a", "Grm1", "Lypd6", "Samd5", "Chrm3", "Angpt1"

seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(20))] <- "LS_AP" # LS_astrocytes_progenitor "Slit2", "Kif21a", "Adgrv1"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(9))] <- "LS_IPC" # LS_astrocytes "Zic4", "Zic1", "Trhde", "Fgd5", "Zic2", "Galnt13"

seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(3))] <- "Endo" # Endothelial_cells "Col4a1", "Kdr", "Fli1", "Lama4", "Ets1", "Itga1", "Fli1", "Rgs5", "Ptprb", "Mecom", "Pecam1", "Cldn5", "Cdh5"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(16))] <- "Fibrob" # Fibroblast "Ahnak", "Col3a1", "Col5a1", "Col6a3", "Dcn", "Lum", "Col6a2", "Col6a1", "Lgals1", "Col1a2", "Svep1"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(28))] <- "MG" # microglia "Ptprc", "Ly86", "Csf1r", "Cd86", "Adgre1"

seurat_object$cell_identity <- factor(seurat_object$cell_identity, levels = c("Cortex_DP_AP", "Cortex_MP_AP", "Cortex_IPC", "Cortex_MpMN", "Cortex_RMSN", "Cortex_Layer_6b", "Cortex_DP_CThPN", "Cortex_MP_CThPN", "Cortex_DP_SCPN", "Cortex_MP_SCPN", "Cortex_UL_CPN", "Cortex_CR", "MGE_InN", "LGE_AP", "LGE_IPC", "LGE_MSN_prog", "CPU_pre_MSN", "CPU_MSN", "CPU_Endo", "ChP_MCC", "LCS_IMN", "CLA_ExN", "LS_AP", "LS_IPC", "Endo", "Fibrob", "MG"))

seurat_object$cell_type_allnames <- as.character(seurat_object$cell_type)
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "AP")] <- "apical_progenitors"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "IPC")] <- "intermediate_progenitor_cells"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "MpMN")] <- "multipolar_migrating_neurons_cells"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "RMSN")] <- "Rostral_Migratory_Stream_neurons_cells"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "CThPN")] <- "corticothalamic_projection_neurons"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "SCPN")] <- "subcerebral_projection_neurons"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "UL_CPU")] <- "upper_layer_callosal_projection_neurons"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "CR")] <- "Cajal-Retzius_cells"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "MGE_InN")] <- "MGE_Interneuron"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "LGE_Astro")] <- "LGE_astrocytes"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "LGE_InN")] <- "LGE_Interneuron"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "LGE_MSN_prog")] <- "LGE_MSN_progenitor"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "CPU_D2_MSN")] <- "CPU_D2_MSN"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "CPU_D1_MSN")] <- "CPU_D1_MSN"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "CPU_Endo")] <- "CPU_Endothelial_cells"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "ChP_MCC")] <- "ChP_multiciliated_cells"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "LCS_IMN")] <- "LCS_immature_neuron"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "CLA_ExN")] <- "CLA_Claustrum_excitatory_neurons"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "LSP")] <- "LS_astrocytes_progenitor"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "LSA")] <- "LS_astrocytes"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "Endo")] <- "Endothelial_cells"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "MG")] <- "microglia"
seurat_object$cell_type_allnames[which(seurat_object$cell_type == "Fibrob")] <- "Fibroblast"

Idents(seurat_object) <- "cell_identity"
DefaultAssay(seurat_object) <- "SCT"

key_markers <- c("Pax6", "Eomes", "Hmga2", "Lrp4", "Dct", "Adamts19", "Eya1", "Nrn1", "Slc17a6", "Sema3c", "Scrt2", "Lhx2", "Bcar1", "Adamts2", "Pou3f2", "Unc5d", "Fgfr1", "Cemip", "Tbr1", "Hs3st4", "Nfe2l3", "Tmem132d", "Pappa2", "Kcnab1", "Abca8a", "Npas2", "Wnt7b", "Klhl1", "Npy", "Kctd12", "Nr4a3", "Slc24a2", "Pde1a", "Fezf2", "Kcnk1", "Cntn6", "Tafa1", "Tcerg1l", "Rorb", "Mical2", "Trmt9b", "Ntf3", "Id2", "Tafa2", "Bcl6", "Inhba", "Reln", "Dgkg", "Elfn1", "Lhx6", "Rai2", "Lhx8", "Nkx2-1", "Pdzph1", "Daam2", "Hes5", "Ascl1", "Dlx6", "Dlx1", "Dlx2", "Sp9", "Slit3", "Six3", "Drd2", "Sox8", "Nebl", "Ikzf1", "Drd1", "Tac1", "Ngef", "Gmnc", "Lmx1a", "Aqp1", "Tfap2d", "Nr2f1", "Nr4a2", "Sv2b", "Slit2", "Zic4", "Trhde", "Galnt13", "Lama4", "Ets1", "Ptprb", "Col3a1", "Col5a1", "Dcn", "Ptprc", "Ly86", "Csf1r")

DotPlot(seurat_object, features = key_markers, cols = c("lightgrey", "red"), col.min = 0.1, dot.min = 0.1) + RotatedAxis() + labs(y = "cell identity", x = "Features") + NoLegend()
#!----------------------------------------------------------------------------------------------------------------------------------------------------------
library(CytoTRACE2)
seurat_object <- cytotrace2(seurat_object, is_seurat = TRUE, slot_type = "counts", species = 'mouse')

p1 <- plot_cytotrace2(seurat_object, reduction = "umap_harmony")

(p1$CytoTRACE2_UMAP + p1$CytoTRACE2_Potency_UMAP + p1$CytoTRACE2_Relative_UMAP) & SetAxes()
ggsave("ST_development_CytoTRACE2.pdf")
Seurat::SpatialPlot(seurat_object, shape = 22, features = "CytoTRACE2_Relative", crop = F) & scale_fill_gradientn(colours = viridis::plasma(15, direction = 1))

Idents(seurat_object) <- "Region"
VlnPlot(seurat_object, features = "CytoTRACE2_Relative", split.by = "orig.ident", cols = c("#34D916", "#00D4E6", "#1E90FF"), pt.size = 0, idents = c("Cortex", "LGE", "CPU", "LS")) & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.1, position = position_dodge(0.9), show.legend = F)
qsave(seurat_object, "output/brain_development_ST_final_tutorial.qs")
#!-----------------------------------------------------------------------------------------------------------------------------------------------------------
#seurat_object$cell_identity1 <- seurat_object$cell_identity
#seurat_object$cell_identity <- as.character(seurat_object$cell_identity)

#seurat_object$cell_identity[which(seurat_object$cell_identity == "LGE_MGE_AP")] <- "LGE_AP"
#seurat_object$cell_identity[which(seurat_object$cell_identity == "Cortex_SCPN" & seurat_object$orig.ident == "E14_5")] <- "Cortex_DP_SCPN"

#seurat_object$cell_identity <- factor(seurat_object$cell_identity, levels = intersect(levels(seurat_object$cell_identity1), unique(seurat_object$cell_identity)))
#!----------------------------------------------------------------------------------------------------------------------------------------------------------
library(ComplexHeatmap)
library(MetaNeighbor)
library(RColorBrewer)
cols = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0 , 1,length = 101)

seurat_object <- qread("output/brain_development_ST_final_tutorial.qs")

seurat_object$sample_Cluster <- paste0(seurat_object$orig.ident, "|", seurat_object$Cluster)

seurat_object$Region_Cluster <- NA
seurat_object$Region_Cluster[which(seurat_object$Cluster %in% paste0("C", c(1:5)))] <- "GZ"
seurat_object$Region_Cluster[which(seurat_object$Cluster %in% paste0("C", c(6:14, 19)))] <- "CP"
seurat_object$Region_Cluster[which(seurat_object$Cluster %in% paste0("C", c(17, 18)))] <- "MGE"
seurat_object$Region_Cluster[which(seurat_object$Cluster %in% paste0("C", c(15, 16, 20, 21)))] <- "LGE"
seurat_object$Region_Cluster[which(seurat_object$Cluster %in% paste0("C", c(22:24)))] <- "CPU"
seurat_object$Region_Cluster[which(seurat_object$Cluster %in% paste0("C", c(25)))] <- "ChP"
seurat_object$Region_Cluster[which(seurat_object$Cluster %in% paste0("C", c(26:29)))] <- "PIR"
seurat_object$Region_Cluster[which(seurat_object$Cluster %in% paste0("C", c(30:32)))] <- "LS"
seurat_object$Region_Cluster[which(seurat_object$Cluster %in% paste0("C", c(33:36)))] <- "Others"
seurat_object$Region_Cluster <- factor(seurat_object$Region_Cluster, levels = c("GZ", "CP", "MGE", "LGE", "CPU", "ChP", "PIR", "LS", "Others"))

DefaultAssay(seurat_object) <- "SCT"
seurat_sce <- as.SingleCellExperiment(seurat_object, assay = "SCT")
global_hvgs = VariableFeatures(seurat_object)

study_ID = "orig.ident"
cell_type = "cell_identity"
study_ID = seurat_sce[[study_ID]]
cell_type = seurat_sce[[cell_type]]
aurocs = MetaNeighborUS(var_genes = global_hvgs, dat = seurat_sce, study_id = study_ID, cell_type = cell_type, fast_version = T)

pdf(file = paste0("MetaNeighborUS_Spatial_ST_cell_identity_heatmap1.pdf"), width = 8, height = 8)
	MetaNeighbor::plotHeatmap(aurocs, cex = 0.3)
dev.off()

study_ID = "orig.ident"
cell_type = "Cluster"
study_ID = seurat_sce[[study_ID]]
cell_type = seurat_sce[[cell_type]]
aurocs = MetaNeighborUS(var_genes = global_hvgs, dat = seurat_sce, study_id = study_ID, cell_type = cell_type, fast_version = T)

pdf(file = paste0("MetaNeighborUS_Spatial_ST_Cluster_heatmap1.pdf"), width = 8, height = 8)
	MetaNeighbor::plotHeatmap(aurocs, cex = 0.3)
dev.off()

CytoTRACE2_score <- c()
for(i in rownames(aurocs)){
	CytoTRACE2_score <- c(CytoTRACE2_score, mean(seurat_object$CytoTRACE2_Relative[which(seurat_object$sample_Cluster == i)]))
}
CytoTRACE2_Score <- as.character(CytoTRACE2_score)
CytoTRACE2_Score[which(CytoTRACE2_score >= 0 & CytoTRACE2_score < 0.1)] <- "CytoT1"
CytoTRACE2_Score[which(CytoTRACE2_score >= 0.1 & CytoTRACE2_score < 0.2)] <- "CytoT2"
CytoTRACE2_Score[which(CytoTRACE2_score >= 0.2 & CytoTRACE2_score < 0.3)] <- "CytoT3"
CytoTRACE2_Score[which(CytoTRACE2_score >= 0.3 & CytoTRACE2_score < 0.4)] <- "CytoT4"
CytoTRACE2_Score[which(CytoTRACE2_score >= 0.4 & CytoTRACE2_score < 0.5)] <- "CytoT5"
CytoTRACE2_Score[which(CytoTRACE2_score >= 0.5 & CytoTRACE2_score < 0.6)] <- "CytoT6"
CytoTRACE2_Score[which(CytoTRACE2_score >= 0.6 & CytoTRACE2_score < 0.7)] <- "CytoT7"
CytoTRACE2_Score[which(CytoTRACE2_score >= 0.7 & CytoTRACE2_score < 0.8)] <- "CytoT8"
CytoTRACE2_Score[which(CytoTRACE2_score >= 0.8 & CytoTRACE2_score < 0.9)] <- "CytoT9"
CytoTRACE2_Score[which(CytoTRACE2_score >= 0.9 & CytoTRACE2_score <= 1)] <- "CytoT10"
CytoTRACE2_Score <- factor(CytoTRACE2_Score, levels = paste0("CytoT", 1:10))

CytoTRACE2_color <- paletteContinuous(set = "beach", n = 10)
names(CytoTRACE2_color) <- levels(CytoTRACE2_Score)

tmp <- as.data.frame(strsplit(colnames(aurocs), split = '\\|'))
sample <- tmp[1,]
names(sample) <- NULL
major_types <- tmp[2,]
names(major_types) <- NULL

sample <- factor(sample, levels(seurat_object$orig.ident))

major_types <- factor(major_types, levels(seurat_object$Cluster))
Region <- unique(seurat_object@meta.data[,c("Region_Cluster", "Cluster")])
Region <- Region$Region_Cluster[match(major_types, Region$Cluster)]

ann_row_1 <- rowAnnotation(a = major_types, border = TRUE, col = list(a = ClusterIDColors[match(levels(major_types), names(ClusterIDColors))]), annotation_label = "Cluster", annotation_legend_param = list(title = "Cluster"))
ann_row_2 <- rowAnnotation(a = sample, border = TRUE, col = list(a = cluster_colors[match(levels(sample), names(cluster_colors))]), annotation_label = "sample", annotation_legend_param = list(title = "sample"))
ann_row_3 <- rowAnnotation(a = Region, border = TRUE, col = list(a = cluster_colors[match(levels(Region), names(cluster_colors))]), annotation_label = "Region", annotation_legend_param = list(title = "Region"))
ann_row_4 <- rowAnnotation(a = CytoTRACE2_Score, border = TRUE, col = list(a = CytoTRACE2_color), annotation_label = "CytoT2_Score", annotation_legend_param = list(title = "CytoT2_Score"))

ann_row <- c(ann_row_2, ann_row_1, ann_row_4, ann_row_3)
#ann_row@gap <- rep(unit(1, "mm"), length(ann_row))

ann_col_1 <- HeatmapAnnotation(a = major_types, border = TRUE, col = list(a = ClusterIDColors[match(levels(major_types), names(ClusterIDColors))]), annotation_label = "Cluster", show_legend = F, annotation_name_side = "left")
ann_col_2 <- HeatmapAnnotation(a = sample, border = TRUE, col = list(a = cluster_colors[match(levels(sample), names(cluster_colors))]), annotation_label = "sample", show_legend = F, annotation_name_side = "left")
ann_col_3 <- HeatmapAnnotation(a = Region, border = TRUE, col = list(a = cluster_colors[match(levels(Region), names(cluster_colors))]), annotation_label = "Region", show_legend = F, annotation_name_side = list(title = "left"))
ann_col_4 <- HeatmapAnnotation(a = CytoTRACE2_Score, border = TRUE, col = list(a = CytoTRACE2_color), annotation_label = "CytoT2_Score", show_legend = F, annotation_name_side = list(title = "left"))

ann_col <- c(ann_col_2, ann_col_1, ann_col_4, ann_col_3)
#ann_col@gap <- rep(unit(1, "mm"), length(ann_col))

ordering <- stats::as.dendrogram(orderCellTypes(aurocs), hang = 1)

ident_labels <- rownames(aurocs)

pdf(paste0("MetaNeighborUS_Spatial_ST_Cluster_heatmap12.pdf"), width = 14, height = 13)

Heatmap(aurocs, col = cols, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8), row_dend_width = unit(70, 'mm'), column_dend_height = unit(70, 'mm'), cluster_rows = ordering, cluster_columns = ordering, 
		left_annotation = ann_row, top_annotation = ann_col, row_labels = ident_labels, column_labels = ident_labels,
		# show_row_names = F, show_column_names = F,
		heatmap_legend_param = list(title_position = "topcenter", legend_width = unit(5, "cm"), legend_height = unit(3, "cm"), title = "AUROC"))

dev.off()

study_ID = "orig.ident"
cell_type = "Region"
study_ID = seurat_sce[[study_ID]]
cell_type = seurat_sce[[cell_type]]
aurocs = MetaNeighborUS(var_genes = global_hvgs, dat = seurat_sce, study_id = study_ID, cell_type = cell_type, fast_version = T)

pdf(file = paste0("MetaNeighborUS_ST_Region_heatmap1.pdf"), width = 8, height = 8)
	MetaNeighbor::plotHeatmap(aurocs, cex = 0.6)
dev.off()

study_ID = "orig.ident"
cell_type = "cell_type"
study_ID = seurat_sce[[study_ID]]
cell_type = seurat_sce[[cell_type]]
aurocs = MetaNeighborUS(var_genes = global_hvgs, dat = seurat_sce, study_id = study_ID, cell_type = cell_type, fast_version = T)

pdf(file = paste0("MetaNeighborUS_ST_cell_type_heatmap1.pdf"), width = 8, height = 8)
	MetaNeighbor::plotHeatmap(aurocs, cex = 0.6)
dev.off()

Idents(seurat_object) <- "cell_identity"

gsea_file <- file("GO_BP_neuron_associated_terms", open="rt")
gsea_gene_set <- list()
count <- 1
while(TRUE){ 
        line <- readLines(gsea_file, n = 1)
        if(length(line) == 0) break
        gene_set <- strsplit(line, "\t")[[1]]
        gsea_gene_set[[count]] <- intersect(gene_set[3:length(gene_set)], rownames(seurat_object))
        names(gsea_gene_set)[count] <- gene_set[1]
        count <- count + 1
}
close(gsea_file)

seurat_object <- UCell::AddModuleScore_UCell(seurat_object, features = gsea_gene_set, ncores = 1, name = "")
pdf("GO_BP_neuron_associated_terms_all_cells.pdf", width = 24)
for(i in names(gsea_gene_set)){
	p0 <- VlnPlot(seurat_object, group.by = "cell_type", split.by = "orig.ident", pt.size = 0, cols = c("#34D916", "#00D4E6", "#1E90FF"), features = i) & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.1, position = position_dodge(0.9), show.legend = F)
	print(p0)
}
dev.off()
#!----------------------------------------------------------------------------------------------------------------
library(ggsankey)

meta_data <- seurat_object@meta.data %>% make_long(Region, cell_type)
meta_data$node <- factor(meta_data$node, levels = rev(c(levels(seurat_object$Region), levels(seurat_object$cell_type))))
meta_data$next_node <- factor(meta_data$next_node, levels = rev(c(levels(seurat_object$Region), levels(seurat_object$cell_type))))

ggplot(meta_data, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) + geom_sankey(color = "black", lwd = 0.2, flow.alpha = 0.5, show.legend = FALSE, na.rm = FALSE) + scale_x_discrete(expand = c(0, 0.8)) + 
		scale_fill_manual(values = cluster_colors, drop = FALSE) + geom_sankey_text(size = 3.2, color = "black") + theme_void() + theme(axis.text.x = element_text())
ggsave("Development_QC_cell_identity/development_cell_type_Region_sankey_plot.pdf")

meta_data <- seurat_object@meta.data %>% make_long(Region, Cluster)
meta_data$node <- factor(meta_data$node, levels = rev(c(levels(seurat_object$Region), levels(seurat_object$Cluster))))
meta_data$next_node <- factor(meta_data$next_node, levels = rev(c(levels(seurat_object$Region), levels(seurat_object$Cluster))))

ggplot(meta_data, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) + geom_sankey(color = "black", lwd = 0.2, flow.alpha = 0.5, show.legend = FALSE, na.rm = FALSE) + scale_x_discrete(expand = c(0, 0.8)) + 
		scale_fill_manual(values = c(cluster_colors, ClusterIDColors), drop = FALSE) + geom_sankey_text(size = 3.2, color = "black") + theme_void() + theme(axis.text.x = element_text())
ggsave("Development_QC_cell_identity/development_Cluster_Region_sankey_plot.pdf")
#!----------------------------------------------------------------------------------------------------------------
Idents(seurat_object) <- "Region"
seurat_object <- PrepSCTFindMarkers(seurat_object)
all_markers <- FindAllMarkers(seurat_object, only.pos = T, min.pct = 0.1)
de_markers <- all_markers %>% filter(avg_log2FC > 1, pct.1 > 0.3, p_val_adj < 0.05)

de_markers$index <- 1:nrow(de_markers)
de_markers <- de_markers[-grep("Rik$|^Gm|^mt-", de_markers$gene),]
de_markers <- de_markers[,-1]
de_markers <- de_markers[order(de_markers$gene, de_markers$avg_log2FC, decreasing = TRUE),]
de_markers <- de_markers[!duplicated(de_markers$gene),]
de_markers <- de_markers[order(de_markers$index),]

features <- de_markers$gene

RNA_Matrix <- AverageExpression(seurat_object, assays = "SCT", features = features, group.by = "Region", slot = "counts", return.seurat = T)
RNA_Matrix$CellType <- colnames(RNA_Matrix)
Idents(RNA_Matrix) <- "CellType"

library(ClusterGVis)
CGVs <- prepareDataFromscRNA(object = RNA_Matrix, diffData = de_markers, showAverage = FALSE, group.by = "CellType", assays = "SCT", slot = "scale.data", scale.data = FALSE)
names(CGVs$wide.res) <- c("gene", levels(RNA_Matrix), "cluster")
CGVs$geneMode <- "average"

key_markers <- c("Veph1", "Emid1", "Dct","Pax6", "Eomes", "Ltbp1", "Slc17a6", "Lhx2", "Neurod2", "Ptger3", "Abca8a", "Tmem132d", "Pappa2", "Nfe2l3", "Nr4a3", "Pde1a", "Fezf2", "Tafa1", "Nwd2", "Ntf3", "L3mbtl4", "Slc26a7", "Rasgrf2", "Tafa2", "Inhba", "Vit", "Notch1", "Lhx6", "Nxph1", "Dlx1", "Slc18a2", "Dlx6", "Sp9", "Six3", "Dlx2", "Drd2", "Isl1", "Ikzf1", "Rxrg", "Cldn5", "Folr1", "Gmnc", "Nr4a2", "Zic4", "Ano1", "Il17rd", "Adgrl4", "Cd93", "Colec12")

pdf("Region_ST_heatmap.pdf", width = 10, height = 12, onefile = F)
visCluster(object = CGVs, ht.col.list = list(col_range = seq(-2, 2, length = 200), col_color = paletteContinuous(set = "solarExtra", n = 200)), plot.type = "both", cluster_rows = TRUE, column_names_rot = 45, show_row_dend = F, markGenes = key_markers, markGenes.side = "left", genes.gp = c('italic',fontsize = 12,col = "black"), sample.col = cluster_colors[match(levels(seurat_object), names(cluster_colors))], ctAnno.col = cluster_colors[match(levels(seurat_object), names(cluster_colors))], use_raster = TRUE, cluster.order =  1:12)
dev.off()
pdf("Region_ST_heatmap1.pdf", width = 8, height = 10, onefile = F)
ht <- GroupHeatmap(seurat_object, features = features, slot = "data", exp_method = "zscore", group.by = "Region", group_palcolor = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))]), cluster_rows = FALSE, cluster_columns = FALSE, heatmap_palette = "viridis", nlabel = 0, show_row_names = FALSE, row_names_side = "left", features_label = key_markers, use_raster = TRUE)
ht$plot
dev.off()

seurat_object$time_Region <- paste0(seurat_object$orig.ident, "_", seurat_object$Region)
seurat_object$time_Region <- factor(seurat_object$time_Region, intersect(paste0(rep(unique(seurat_object$orig.ident), 3), "_", rep(levels(seurat_object$Region), each = 3)), unique(seurat_object$time_Region)))

Idents(seurat_object) <- "time_Region"
all_markers <- FindAllMarkers(seurat_object, only.pos = T, min.pct = 0.1)
de_markers <- all_markers %>% filter(avg_log2FC > 0.5, pct.1 > 0.3, p_val_adj < 0.05)

de_markers$index <- 1:nrow(de_markers)
de_markers <- de_markers[-grep("Rik$|^Gm|^mt-", de_markers$gene),]
de_markers <- de_markers[,-1]
de_markers <- de_markers[order(de_markers$gene, de_markers$avg_log2FC, decreasing = TRUE),]
de_markers <- de_markers[!duplicated(de_markers$gene),]
de_markers <- de_markers[order(de_markers$index),]
de_markers$final_cluster <- gsub(de_markers$cluster, pattern = "E1[468]_5_", replacement = "")
write.csv(de_markers, "Time_Region_ST.csv", quote = F, row.names = F)
pdf("Time_Region_ST_heatmap1.pdf", width = 12, height = 12, onefile = F)
ht <- SCP::GroupHeatmap(seurat_object, features = de_markers$gene, slot = "data", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, feature_split = de_markers$final_cluster, feature_split_palcolor = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))]), heatmap_palette = "viridis", nlabel = 0, show_row_names = FALSE, row_names_side = "left", features_label = key_markers, use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot
dev.off()

Idents(seurat_object) <- "Region"
# identity conserved markers cross time
conserved_markers <- list()
for(i in levels(seurat_object$Region)){
	conserved_markers[[i]] <- FindConservedMarkers(seurat_object, ident.1 = i, grouping.var = "orig.ident", logfc.threshold = 0.25, min.pct = 0.1, min.diff.pct = 0.1, only.pos = T, verbose = FALSE, assay = "SCT")
}
for(i in levels(seurat_object$Region)){
	if(length(grep("E14_5_", colnames(conserved_markers[[i]]))) == 0){
		conserved_markers[[i]]$E14_5_p_val <- NA
		conserved_markers[[i]]$E14_5_avg_log2FC <- NA
		conserved_markers[[i]]$E14_5_pct.1 <- NA
		conserved_markers[[i]]$E14_5_pct.2 <- NA
		conserved_markers[[i]]$E14_5_p_val_adj <- 0
	}
	if(length(grep("E16_5_", colnames(conserved_markers[[i]]))) == 0){
		conserved_markers[[i]]$E16_5_p_val <- NA
		conserved_markers[[i]]$E16_5_avg_log2FC <- NA
		conserved_markers[[i]]$E16_5_pct.1 <- NA
		conserved_markers[[i]]$E16_5_pct.2 <- NA
		conserved_markers[[i]]$E16_5_p_val_adj <- 0
	}
	if(length(grep("E18_5_", colnames(conserved_markers[[i]]))) == 0){
		conserved_markers[[i]]$E18_5_p_val <- NA
		conserved_markers[[i]]$E18_5_avg_log2FC <- NA
		conserved_markers[[i]]$E18_5_pct.1 <- NA
		conserved_markers[[i]]$E18_5_pct.2 <- NA
		conserved_markers[[i]]$E18_5_p_val_adj <- 0
	}
	conserved_markers[[i]] <- conserved_markers[[i]][,c(paste0(rep(c("E14_5_", "E16_5_", "E18_5_"), each = 5), rep(c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"), 3)))]
	conserved_markers[[i]]$gene <- rownames(conserved_markers[[i]])
}

all_conserved_markers <- c()
for(i in levels(seurat_object$Region)){
	conserved_markers[[i]]$cluster <- i
	all_conserved_markers <- rbind(all_conserved_markers, conserved_markers[[i]])
}
for(i in 1:nrow(all_conserved_markers)){
	all_conserved_markers$max_pval_adj[i] <- max(all_conserved_markers[i,c("E14_5_p_val_adj", "E16_5_p_val_adj", "E18_5_p_val_adj")])
}
all_conserved_markers <- all_conserved_markers[which(all_conserved_markers$max_pval_adj < 0.05),]
de_conserved_markers <- all_conserved_markers[which(all_conserved_markers$E14_5_avg_log2FC > 0.5 | all_conserved_markers$E16_5_avg_log2FC > 0.5 | all_conserved_markers$E18_5_avg_log2FC > 0.5),]

all_conserved_markers <- FindAllMarkers(seurat_object, features = unique(de_conserved_markers$gene))
all_conserved_markers$index <- 1:nrow(all_conserved_markers)
all_conserved_markers <- all_conserved_markers[order(all_conserved_markers$gene, all_conserved_markers$avg_log2FC, decreasing = TRUE),]
all_conserved_markers <- all_conserved_markers[!duplicated(all_conserved_markers$gene),]
all_conserved_markers <- all_conserved_markers[order(all_conserved_markers$index),]

de_conserved_markers$final_cluster <- all_conserved_markers$cluster[match(de_conserved_markers$gene, all_conserved_markers$gene)]
write.csv(de_conserved_markers, "Region_all_conserved_markers_cell_type.csv", quote = F, row.names = F)

ht <- SCP::GroupHeatmap(seurat_object, features = all_conserved_markers$gene, slot = "data", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, feature_split = all_conserved_markers$cluster, feature_split_palcolor = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))]), heatmap_palette = "viridis", nlabel = 0, show_row_names = FALSE, row_names_side = "left", features_label = key_markers, use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot
pdf("Region_cell_type_conserved_features.pdf", width = 12, height = 12)
ht$plot
dev.off()
#!----------------------------------------------------------------------------------------------------------------------------------------------------------
Idents(seurat_object) <- "cell_type"
DefaultAssay(seurat_object) <- "SCT"
AP_markers <- FindMarkers(seurat_object, ident.1 = c("AP", "LGE_AP", "LGE_MGE_AP", "LS_AP"), only.pos = T, min.pct = 0.1)
AP_markers <- AP_markers %>% filter(avg_log2FC > 1, pct.1 > 0.5, p_val_adj < 0.05)
AP_markers <- AP_markers[-grep("Rik$", rownames(AP_markers)),]

ht <- SCP::GroupHeatmap(seurat_object, features = rownames(AP_markers), slot = "counts", exp_method = "zscore", group.by = "cell_type", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$cell_type), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, feature_split_palcolor = list(cluster_colors[match(levels(seurat_object$cell_type), names(cluster_colors))]), heatmap_palette = "viridis", nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot

dir.create("AP_Specific_Markers")
DefaultAssay(seurat_object) <- "spARC_SCT"
for(i in rownames(AP_markers)){
	try({
		p0 <- Seurat::SpatialPlot(seurat_object, features = i, shape = 22, crop = F) & scale_fill_gradientn(colours = paletteContinuous(set = "sambaNight"))
		print(p0)
		ggsave(paste0("AP_Specific_Markers/", i, "_rna.png"), width = 24, height = 24)
		p1 <- Seurat::SpatialPlot(hic_object, shape = 22, crop = F, features = hic_object@misc$gene_infor_500k$gene_region[which(hic_object@misc$gene_infor_500k$gene_name == i)]) & scale_fill_gradientn(values = seq(0, 1, 0.01), colours = paletteContinuous(set = "blueYellow", n = 100))
		print(p1)
		ggsave(paste0("AP_Specific_Markers/", i, "_hic.png"), width = 24, height = 24)
	})
}
#!----------------------------------------------------------------------------------------------------------------------------------------------------------
all_markers <- FindAllMarkers(seurat_object, only.pos = T, min.pct = 0.1)
de_markers <- all_markers %>% filter(avg_log2FC > 1, pct.1 > 0.3, p_val_adj < 0.05)
de_markers$index <- 1:nrow(de_markers)
de_markers <- de_markers[order(de_markers$gene, de_markers$avg_log2FC, decreasing = TRUE),]
de_markers <- de_markers[!duplicated(de_markers$gene),]
de_markers <- de_markers[order(de_markers$index),]
de_markers$final_cluster <- de_markers$cluster[match(de_markers$gene, de_markers$gene)]
de_markers <- de_markers[-grep("Gm[0-9]*|Rik$|Rik[0-9]$", de_markers$gene),]
write.csv(de_markers, "all_cell_identity_DE_markers.csv", quote = F, row.names = F)
de_markers$final_cluster <- names(ClusterNames)[match(de_markers$final_cluster, ClusterNames)]
de_markers$final_cluster <- factor(de_markers$final_cluster, levels = levels(seurat_object$Cluster))

key_markers <- c("Lpar1", "Slc15a2", "Emid1", "Dct","Pax6", "Eomes", "Slc17a6", "Nrn1", "Adamts2", "Lypd6", "Cyp7b1", "Abca8a", "Nr4a2", "Hs3st3b1", "Wnt7b", "Klhl1", "Npy", "Nr4a3", "Pde1a", "Fezf2", "Dpf3", "Adamts17", "Kit", "Tafa1", "Ntf3", "Trmt9b", "Rasgrf2", "Inhba", "Tafa2", "Reln", "Dtl","Hmga2", "Pdzph1", "Aldh1l1", "Lhx6", "Nxph1", "Dlx6", "Sp9", "Dlx2", "Six3", "Drd2", "Isl1", "Ikzf1", "Rarb", "Rxrg", "Cldn5", "Folr1", "Gmnc", "Tfap2d", "Ppp1r14c","Unc13c", "Lypd6b", "Zic4", "Ano1", "Il17rd", "Adgrl4", "Cd93", "Dcn", "Col5a1", "Spp1", "C1qb", "Colec12")
DotPlot(seurat_object, features = key_markers, cols = cluster_colors[1:3], dot.scale = 8, split.by = "orig.ident", dot.min = 0.3, col.min = 0.2) + coord_flip() + RotatedAxis() + NoLegend()
ggsave("key_markers_for_cell_identity.pdf")
DotPlot(seurat_object, features = key_markers, cols = cluster_colors[1:3], group.by = "Cluster", dot.scale = 8, split.by = "orig.ident", dot.min = 0.3, col.min = 0.2) + coord_flip() + RotatedAxis()
ggsave("key_markers_for_Cluster.pdf")
DotPlot(seurat_object, features = key_markers, cols = c("lightgrey", "red"), group.by = "Cluster", dot.scale = 8, dot.min = 0.3, col.min = 0.2) + coord_flip() + RotatedAxis()
ggsave("key_markers_for_Cluster_all.pdf", width = 14, height = 17)
DotPlot(seurat_object, features = key_markers, cols = c("lightgrey", "red"), group.by = "cell_identity", dot.scale = 8, dot.min = 0.3, col.min = 0.2) + coord_flip() + RotatedAxis()
ggsave("key_markers_for_cell_identity_all.pdf", width = 14, height = 20)

key_markers1 <- c("Veph1", "Emid1", "Dct","Pax6", "Eomes", "Ltbp1", "Slc17a6", "Lhx2", "Neurod2", "Ptger3", "Abca8a", "Tmem132d", "Pappa2", "Nfe2l3", "Nr4a3", "Pde1a", "Fezf2", "Tafa1", "Nwd2", "Ntf3", "L3mbtl4", "Slc26a7", "Rasgrf2", "Tafa2", "Inhba", "Vit", "Notch1", "Lhx6", "Nxph1", "Dlx1", "Slc18a2", "Dlx6", "Sp9", "Six3", "Dlx2", "Drd2", "Isl1", "Ikzf1", "Rxrg", "Cldn5", "Folr1", "Gmnc", "Nr4a2", "Zic4", "Ano1", "Il17rd", "Adgrl4", "Cd93", "Colec12")
key_markers <- unique(c(key_markers, key_markers1, "Veph1", "Emid1", "Dct","Pax6", "Eomes", "Ltbp1", "Slc17a6", "Neurod2", "Ptger3", "Abca8a", "Trpc3", "Met", "Nr4a2", "Tmem132d", "Pappa2", "Hs3st3b1", "Wnt7b", "Ntsr1", "Npas2", "Klhl1", "Npy", "Nr4a3", "Pde1a", "Fezf2", "Tafa1", "Nwd2", "Tcerg1l", "Adamts20", "L3mbtl4", "Ntf3", "Trmt9b", "Slc26a7", "Rasgrf2", "Tafa2", "Inhba", "Dtl","Hmga2", "Pdzph1", "Mt3", "Aldh1l1", "Lhx6", "Nxph1", "Six3", "Dlx1", "Dlx6", "Sp9", "Dlx2", "Drd2", "Isl1", "Ikzf1", "Rxrg", "Cldn5", "Folr1", "Gmnc", "Zic4", "Ano1", "Il17rd", "Adgrl4", "Cd93", "Colec12", "Satb2", "Dlx2", "Gad2", "Foxp1", "Pde10a", "Sox2", "Hes5", "Pax6", "Notch1", "Notch2", "Rbfox1", "Ank3", "Celf2", "Zbtb20"))

ht <- SCP::GroupHeatmap(seurat_object, features = de_markers$gene, feature_split = de_markers$final_cluster, slot = "counts", assay = "SCT", exp_method = "zscore", group.by = "Cluster", group_palcolor = list(ClusterIDColors), feature_split_palcolor = list(ClusterIDColors), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palcolor = c("#0099CC", "white", "#CC0033"), nlabel = 0, show_row_names = F, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE), features_label = key_markers)
ht$plot
pdf("all_DE_Clusters.pdf", width = 16, height = 18)
ht$plot
dev.off()

DefaultAssay(seurat_object) <- "spARC_SCT"
DefaultAssay(hic_object) <- "spARC_scAB250kb_scale"
dir.create("key_features")
for(i in key_markers){
	try({
		p0 <- Seurat::SpatialPlot(seurat_object, features = i, shape = 22, crop = F) & scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")))
		p1 <- Seurat::SpatialPlot(hic_object, shape = 22, crop = F, features = unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name == i & gene_infor_250k$is_promoter == "promoter")])) & scale_fill_gradientn(values = seq(0, 1, 0.01), colours = c("#1CA8FF","#3AD6FE","#91E9FF","#E9F0C4","#FFDA3B","#FF3606","#E30407"))
		print(p0/p1)
		ggsave(paste0("key_features/", i, "_1.pdf"), width = 24, height = 24)
		p0 <- FeaturePlot(seurat_object, reduction = "spatial_umap", features = i) & scale_colour_gradientn(colours = paletteContinuous(set = "solarExtra"))
		p1 <- FeaturePlot(hic_object, features = unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name == i & gene_infor_250k$is_promoter == "promoter")]), reduction = "spatial_umap") & scale_colour_gradientn(colours = paletteContinuous(set = "solarExtra"))
		print(p0/p1)
		ggsave(paste0("key_features/", i, "_2.pdf"), width = 20, height = 10)
	})
}

DefaultAssay(seurat_object) <- "SCT"
DefaultAssay(hic_object) <- "scAB250kb_scale"
for(i in key_markers){
	try({
		p0 <- VlnPlot(seurat_object, features = i, pt.size = 0, group.by = "Cluster", cols = ClusterIDColors) + NoLegend()
		p1 <- VlnPlot(hic_object, pt.size = 0, features = unique(gene_infor_250k$gene_region[which(gene_infor_250k$gene_name == i & gene_infor_250k$is_promoter == "promoter")]), group.by = "banksy_cluster", cols = hic_object@misc$banksy_cluster_colors) + ylab("scAB values") & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.3, position = position_dodge(0.9), show.legend = F) & NoLegend()
		print(p0 + p1)
		ggsave(paste0("key_features/", i, "_vlnplot.pdf"), width = 20, height = 14)
	})
}
#!----------------------------------------------------------------------------------------------------------------------------------------------------------
# cellhint harmonisation with development data and identity the conserved cell type(cellhint group) markers

cellhint_group <- read.table("cellhint/cellhint_group", header = T, sep = "\t")
seurat_object$cellhint_group_name <- cellhint_group$cellhint_group_name[match(seurat_object$high_hierarchy, cellhint_group$high_hierarchy)]
seurat_object$cellhint_group_name <- factor(seurat_object$cellhint_group_name, levels = intersect(levels(seurat_object$cell_identity), unique(seurat_object$cellhint_group_name)))
seurat_object$cellhint_group <- names(ClusterNames)[match(seurat_object$cellhint_group_name, ClusterNames)]
seurat_object$cellhint_group <- factor(seurat_object$cellhint_group, levels = intersect(ClusterID, unique(seurat_object$cellhint_group)))
cellhint_group_ID <- paste0("group", 1:length(levels(seurat_object$cellhint_group)))
names(cellhint_group_ID) <- levels(seurat_object$cellhint_group)
cellhint_group_ID <- cellhint_group_ID[match(seurat_object$cellhint_group, names(cellhint_group_ID))]
names(cellhint_group_ID) <- NULL
seurat_object$cellhint_group_ID <- cellhint_group_ID
seurat_object$cellhint_group_ID <- factor(seurat_object$cellhint_group_ID, levels = paste0("group", 1:length(levels(seurat_object$cellhint_group))))

seurat_object$cellhint <- paste0(seurat_object$cellhint_group, "_", seurat_object$cellhint_group_name)
seurat_object$cellhint <- factor(seurat_object$cellhint, levels = paste0(levels(seurat_object$cellhint_group), "_", levels(seurat_object$cellhint_group_name)))

Idents(seurat_object) <- "cellhint"

conserved_markers <- list()
for(i in levels(seurat_object$cellhint)){
	conserved_markers[[i]] <- FindConservedMarkers(seurat_object, ident.1 = i, grouping.var = "orig.ident", logfc.threshold = 0.25, min.pct = 0.1, min.diff.pct = 0.1, only.pos = T, verbose = FALSE)
}

for(i in levels(seurat_object$cellhint)){
	if(length(grep("E14_5_", colnames(conserved_markers[[i]]))) == 0){
		conserved_markers[[i]]$E14_5_p_val <- NA
		conserved_markers[[i]]$E14_5_avg_log2FC <- NA
		conserved_markers[[i]]$E14_5_pct.1 <- NA
		conserved_markers[[i]]$E14_5_pct.2 <- NA
		conserved_markers[[i]]$E14_5_p_val_adj <- NA
	}
	if(length(grep("E16_5_", colnames(conserved_markers[[i]]))) == 0){
		conserved_markers[[i]]$E16_5_p_val <- NA
		conserved_markers[[i]]$E16_5_avg_log2FC <- NA
		conserved_markers[[i]]$E16_5_pct.1 <- NA
		conserved_markers[[i]]$E16_5_pct.2 <- NA
		conserved_markers[[i]]$E16_5_p_val_adj <- NA
	}
	if(length(grep("E18_5_", colnames(conserved_markers[[i]]))) == 0){
		conserved_markers[[i]]$E18_5_p_val <- NA
		conserved_markers[[i]]$E18_5_avg_log2FC <- NA
		conserved_markers[[i]]$E18_5_pct.1 <- NA
		conserved_markers[[i]]$E18_5_pct.2 <- NA
		conserved_markers[[i]]$E18_5_p_val_adj <- NA
	}
	conserved_markers[[i]] <- conserved_markers[[i]][,c(paste0(rep(c("E14_5_", "E16_5_", "E18_5_"), each = 5), rep(c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"), 3)), "max_pval", "minimump_p_val")]
	conserved_markers[[i]]$gene <- rownames(conserved_markers[[i]])
}
all_conserved_markers <- c()
for(i in levels(seurat_object$cellhint)){
	conserved_markers[[i]]$cluster <- i
	all_conserved_markers <- rbind(all_conserved_markers, conserved_markers[[i]])
}
all_conserved_markers <- all_conserved_markers[which(all_conserved_markers$max_pval < 0.05),]
write.csv(all_conserved_markers, "cellhint/all_conserved_markers.csv", quote = F, row.names = F)

key_markers <- c("Veph1", "Emid1", "Dct","Pax6", "Eomes", "Ltbp1", "Slc17a6", "Lhx2", "Neurod2", "Ptger3", "Abca8a", "Tmem132d", "Pappa2", "Nfe2l3", "Nr4a3", "Pde1a", "Fezf2", "Tafa1", "Nwd2", "Ntf3", "L3mbtl4", "Slc26a7", "Rasgrf2", "Tafa2", "Inhba", "Vit", "Notch1", "Lhx6", "Nxph1", "Dlx1", "Slc18a2", "Dlx6", "Sp9", "Six3", "Dlx2", "Drd2", "Isl1", "Ikzf1", "Rxrg", "Cldn5", "Folr1", "Gmnc", "Nr4a2", "Zic4", "Ano1", "Il17rd", "Adgrl4", "Cd93", "Colec12")
DotPlot(seurat_object, features = key_markers, cols = cluster_colors[1:3], dot.scale = 8, split.by = "orig.ident", dot.min = 0.2, col.min = 0.2) + coord_flip() + RotatedAxis() + NoLegend()
ggsave("Development_QC_cell_identity/conserved_marker_with_cellhint_group.pdf")

ht <- GroupHeatmap(seurat_object, features = key_markers, slot = "data", exp_method = "zscore", group.by = "cell_identity", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$cell_identity), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = FALSE, cluster_columns = FALSE, add_dot = TRUE, heatmap_palette = "viridis", nlabel = 0, show_row_names = TRUE, row_names_side = "left")
ht$plot

features <- unique(all_conserved_markers$gene)
seurat_selected <- ScaleData(seurat_object, features = features, assay = "SCT", scale.max = 10, split.by = "orig.ident")
colnames(seurat_selected@meta.data)[13] <- "CellCycle"
pdf("Development_QC_cell_identity/cellhint.pdf", width = 27, height = 20)
plot_doheatmap(dataset = seurat_selected, markers = features, sort_var = "cellhint_group_name", anno_var = c("cellhint_group_name", "cell_identity", "Cluster", "orig.ident", "CellCycle"), hm_limit = seq(-2, 2, length = 256), hm_colors = paletteContinuous(set = "solarExtra"), anno_colors = list(cluster_colors[match(levels(seurat_selected$cellhint_group_name), names(cluster_colors))], cluster_colors[match(levels(seurat_selected$cell_identity), names(cluster_colors))], ClusterIDColors[match(levels(seurat_selected$Cluster), names(ClusterIDColors))], cluster_colors[1:3], c("#F8766D","#00BA38","#619CFF")), column_split = "cellhint_group_name",
				row_font_size = 10, label_markers = key_markers, cor_row_label_line = "blue", lwd_row_label_line = 0.3)	
dev.off()
#!----------------------------------------------------------------------------------------------------------------------------------------------------------
seurat_object$time_cell_identity <- paste0(seurat_object$orig.ident, "_", seurat_object$cell_identity)
seurat_object$time_cell_identity <- factor(seurat_object$time_cell_identity, levels = intersect(paste0(rep(unique(seurat_object$orig.ident), length(levels(seurat_object$cell_identity))), "_", rep(levels(seurat_object$cell_identity), each = 3)), unique(seurat_object$time_cell_identity)))

selected_clusters <- levels(seurat_object$cell_identity)[grep("Cortex_", levels(seurat_object$cell_identity))]
selected_clusters <- setdiff(selected_clusters, c("Cortex_CR", "Cortex_MGE_derived_InN"))
ggplot(seurat_object@meta.data[which(seurat_object$cell_identity %in% selected_clusters),], aes(time_cell_identity, fill = Phase)) + geom_bar(position = "fill") + ylab("percent of cells") + xlab("") + 
    theme(plot.background = element_rect(fill = NA), legend.key = element_rect(size = 2), legend.background = element_rect(color = NA), axis.text.y = element_text(size = 10, face='bold'), 
    axis.text.x = element_text(size = 10, face='bold', angle = 45, hjust = 1, vjust = 1)) +  
    scale_fill_manual(values = c("#F8766D","#00BA38","#619CFF"))
ggsave("Development_QC_cell_identity/cell_cycle_distribution_with_cell_identity_during_cortex_development.pdf")

seurat_object$time_cell_identity <- factor(seurat_object$time_cell_identity, levels = intersect(paste0(rep(unique(seurat_object$orig.ident), each = length(levels(seurat_object$cell_identity))), "_", rep(levels(seurat_object$cell_identity), 3)), unique(seurat_object$time_cell_identity)))
ggplot(seurat_object@meta.data[which(seurat_object$cell_identity %in% selected_clusters),], aes(time_cell_identity, fill = Phase)) + geom_bar(position = "fill") + ylab("percent of cells") + xlab("") + 
    theme(plot.background = element_rect(fill = NA), legend.key = element_rect(size = 2), legend.background = element_rect(color = NA), axis.text.y = element_text(size = 10, face='bold'), 
    axis.text.x = element_text(size = 10, face='bold', angle = 45, hjust = 1, vjust = 1)) +  
    scale_fill_manual(values = c("#F8766D","#00BA38","#619CFF"))
ggsave("Development_QC_cell_identity/cell_cycle_distribution_with_cell_identity_during_cortex_development1.pdf")

tmp_data <- table(seurat_object$time_cell_identity, seurat_object$Phase)
tmp_data <- tmp_data/rowSums(tmp_data)
meta_data <- data.frame(celltype = rep(colnames(tmp_data), each = nrow(tmp_data)), values = as.vector(tmp_data), group = rep(rownames(tmp_data), ncol(tmp_data)))
meta_data$celltype <- factor(meta_data$celltype, levels = colnames(tmp_data))
meta_data$group <- factor(meta_data$group, levels = rownames(tmp_data))

ggplot(meta_data, aes(x = group, y = values, fill = celltype)) + geom_bar(stat = "identity", position = "dodge") + ylab("percent of cells") + xlab("") + geom_text(aes(label = paste0(round(meta_data$values * 100, 1), "%")), position = position_dodge(0.9), vjust = -0.25) +
    theme(plot.background = element_rect(fill = NA), legend.position = "right", legend.key = element_rect(size = 2), legend.background = element_rect(color = NA), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12)) +
    scale_fill_manual(values = c("#F8766D","#00BA38","#619CFF"))
#!----------------------------------------------------------------------------------------------------------------------------------------------------------
manual_colors <- list(cluster_colors = as.list(cluster_colors), ClusterIDColors = as.list(ClusterIDColors), domainIDColors = as.list(domainIDColors))
node_colors <- data.frame(dataset = "E14_5", cell_type = levels(seurat_list$E14_5), ClusterID = names(ClusterNames)[match(levels(seurat_list$E14_5), ClusterNames)], cluster_colors = cluster_colors[match(levels(seurat_list$E14_5), names(cluster_colors))])
node_colors <- rbind(node_colors, data.frame(dataset = "E16_5", cell_type = levels(seurat_list$E16_5), ClusterID = names(ClusterNames)[match(levels(seurat_list$E16_5), ClusterNames)], cluster_colors = cluster_colors[match(levels(seurat_list$E16_5), names(cluster_colors))]))
node_colors <- rbind(node_colors, data.frame(dataset = "E18_5", cell_type = levels(seurat_list$E18_5), ClusterID = names(ClusterNames)[match(levels(seurat_list$E18_5), ClusterNames)], cluster_colors = cluster_colors[match(levels(seurat_list$E18_5), names(cluster_colors))]))
seurat2scanpy(seurat_object, assay = "SCT", manual_color = manual_colors, major_umap = "umap_harmony", others = node_colors, savefile = paste0("Trajectory_analysis/anndata/brain_development_ST_final.h5ad"))
#!----------------------------------------------------------------------------------------------------------------------------------------------------------
tmp_object <- subset(seurat_object, cells = colnames(seurat_object)[-which(seurat_object$cell_identity %in% c("ChP_MCC", "CPU_Endo", "Endo", "Fibrob", "MG", "Stromal"))])
source("~/scripts/seurat2scanpy/shiny_st.R")
options(browser = "/usr/bin/firefox")
tmp_object <- shiny_st(seurat = tmp_object, isVisium = F, assay = "SCT", image = "E14_5")

tmp_object <- subset(tmp_object, subset = cell_identity != "removed_cells")
tmp_object$cell_identity <- factor(tmp_object$cell_identity, levels = intersect(levels(tmp_object$cell_identity), unique(tmp_object$cell_identity)))
Idents(tmp_object) <- "cell_identity"

DefaultAssay(tmp_object) <- "SCT"
write.csv(tmp_object@meta.data, "Trajectory_analysis/anndata/selected_cells_meta.csv", quote = F)
seurat2scanpy(tmp_object, assay = "SCT", manual_color = manual_colors, major_umap = "umap_harmony", savefile = paste0("Trajectory_analysis/anndata/brain_development_ST_final_selected_cells.h5ad"))

seurat_list$E14_5 <- RenameCells(seurat_list$E14_5, add.cell.id = "E14_5")
seurat_list$E14_5 <- subset(seurat_list$E14_5, cells = colnames(tmp_object)[tmp_object$orig.ident == "E14_5"])
seurat_list[[1]]$cluster <- as.character(seurat_list[[1]]$cell_identity)
seurat_list[[1]]$cluster[which(seurat_list[[1]]$cluster %in% c("Cortex_DP_AP", "Cortex_MP_AP", "LGE_MGE_AP"))] <- "start_clusters"
seurat_list$E16_5 <- RenameCells(seurat_list$E16_5, add.cell.id = "E16_5")
seurat_list$E16_5 <- subset(seurat_list$E16_5, cells = colnames(tmp_object)[tmp_object$orig.ident == "E16_5"])
seurat_list[[2]]$cluster <- as.character(seurat_list[[2]]$cell_identity)
seurat_list[[2]]$cluster[which(seurat_list[[2]]$cluster %in% c("Cortex_DP_AP", "Cortex_MP_AP", "LGE_AP"))] <- "start_clusters"
seurat_list$E18_5 <- RenameCells(seurat_list$E18_5, add.cell.id = "E18_5")
seurat_list$E18_5 <- subset(seurat_list$E18_5, cells = colnames(tmp_object)[tmp_object$orig.ident == "E18_5"])
seurat_list[[3]]$cluster <- as.character(seurat_list[[3]]$cell_identity)
seurat_list[[3]]$cluster[which(seurat_list[[3]]$cluster %in% c("Cortex_DP_AP", "Cortex_MP_AP", "LGE_AP"))] <- "start_clusters"
seurat2scanpy(seurat_list$E14_5, assay = "SCT", manual_color = manual_colors, major_umap = "SCT_umap", savefile = paste0("Trajectory_analysis/anndata/brain_development_ST_E14_5_final_selected_cells.h5ad"))
seurat2scanpy(seurat_list$E16_5, assay = "SCT", manual_color = manual_colors, major_umap = "SCT_umap", savefile = paste0("Trajectory_analysis/anndata/brain_development_ST_E16_5_final_selected_cells.h5ad"))
seurat2scanpy(seurat_list$E18_5, assay = "SCT", manual_color = manual_colors, major_umap = "SCT_umap", savefile = paste0("Trajectory_analysis/anndata/brain_development_ST_E18_5_final_selected_cells.h5ad"))
#!----------------------------------------------------------------------------------------------------------------------------------------------------------
seurat_object <- subset(tmp_object, cells = colnames(tmp_object)[grep("Cortex_", tmp_object$cell_identity)])
seurat_object <- subset(seurat_object, subset = cell_identity != "Cortex_CR")
seurat_object <- subset(seurat_object, subset = cell_identity != "Cortex_MGE_derived_InN")
seurat_object <- shiny_st(seurat = seurat_object, isVisium = F, assay = "SCT", image = "E14_5")
seurat_object <- subset(seurat_object, subset = banksy_cluster %in% c("D1", "D2", "D3", "D4", "D5", "D6"))
seurat_object <- subset(seurat_object, subset = orig.ident != "removed_cells")
seurat_object$cell_identity <- factor(seurat_object$cell_identity, levels = intersect(levels(seurat_object$cell_identity), unique(seurat_object$cell_identity)))
seurat_object$banksy_cluster <- factor(seurat_object$banksy_cluster, levels = intersect(levels(seurat_object$banksy_cluster), unique(seurat_object$banksy_cluster)))
seurat_object$banksy_domains <- factor(seurat_object$banksy_domains, levels = intersect(levels(seurat_object$banksy_domains), unique(seurat_object$banksy_domains)))

seurat_object$Region <- as.character(seurat_object$cell_identity)
seurat_object$Region[which(seurat_object$Region %in% c("Cortex_DP_AP", "Cortex_MP_AP"))] <- "VZ"
seurat_object$Region[which(seurat_object$Region %in% c("Cortex_IPC"))] <- "SVZ"
seurat_object$Region[which(seurat_object$Region %in% c("Cortex_MigN", "Cortex_DP_MigN", "Cortex_MP_MigN"))] <- "IZ"
seurat_object$Region[which(seurat_object$Region %in% c("Cortex_Layer_6b"))] <- "Layer_VIb"
seurat_object$Region[which(seurat_object$Region %in% c("Cortex_CThPN", "Cortex_DP_CThPN", "Cortex_MP_CThPN"))] <- "Layer_VI"
seurat_object$Region[which(seurat_object$Region %in% c("Cortex_SCPN", "Cortex_DP_SCPN", "Cortex_MP_SCPN"))] <- "Layer_V"
seurat_object$Region[which(seurat_object$Region %in% c("Cortex_UL_CPN"))] <- "Layer_II_III_IV"

seurat_object <- shiny_st(seurat = seurat_object, isVisium = F, assay = "SCT", image = "E14_5")
seurat_object$Region <- factor(seurat_object$Region, levels = c("VZ0", "VZ", "SVZ", "IZ", "Layer_VIb", "Layer_VI", "Layer_V", "Layer_II_III_IV"))

temp_object <- CreateSeuratObject(GetAssayData(seurat_object, assay = "RNA", layer = "counts"), min.cells = 1, min.features = 1)
temp_object@meta.data <- seurat_object@meta.data
temp_object[["spatial_umap"]] <- seurat_object[["spatial_umap"]]
temp_object@images <- seurat_object@images
seurat_object <- temp_object

seurat_object <- Split_Layers(seurat_object, split.by = "orig.ident")
seurat_object <- SCTransform(seurat_object, vst.flavor = "v2", method = "glmGamPoi", vars.to.regress = c("S.Score", "G2M.Score"))
seurat_object <- RunPCA(seurat_object, npcs = 30, verbose = FALSE)
# one-liner to run Integration
seurat_object <- IntegrateLayers(object = seurat_object, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = 'harmony', normalization.method = "SCT", verbose = FALSE)
seurat_object <- RunUMAP(seurat_object, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony", return.model = T)
seurat_object <- FindNeighbors(seurat_object, dims = 1:30, reduction = 'harmony', annoy.metric = "cosine", graph.name = c('harmony_nn', 'harmony_snn'))
seurat_object <- FindClusters(seurat_object, cluster.name = "SCT_harmony_cluster", graph.name = 'harmony_snn')
Idents(seurat_object) <- "cell_identity"

DefaultAssay(seurat_object) <- "SCT"
seurat_object <- JoinLayers(seurat_object, assay = "RNA")

library(CytoTRACE2)
seurat_object <- cytotrace2(seurat_object, is_seurat = TRUE, slot_type = "counts", species = 'mouse')

p1 <- plot_cytotrace2(seurat_object, reduction = "umap_harmony")

(p1$CytoTRACE2_UMAP + p1$CytoTRACE2_Potency_UMAP + p1$CytoTRACE2_Relative_UMAP) & SetAxes()
ggsave("Development_QC_cell_identity/Cortex_CytoTRACE2.pdf")

seurat_object <- spARC_Seurat(seurat_object, group_by = "orig.ident", assay = "SCT")
qsave(seurat_object, "output/brain_development_ST_final_tutorial_cortex.qs")
seurat2scanpy(seurat_object, assay = "SCT", manual_color = manual_colors, major_umap = "umap_harmony", savefile = paste0("Trajectory_analysis/anndata/brain_development_ST_final_cortex.h5ad"))

DefaultAssay(seurat_object) <- "SCT"
seurat_object$time_Region <- paste0(seurat_object$orig.ident, "_", seurat_object$Region)
seurat_object$time_Region <- factor(seurat_object$time_Region, levels = c("E14_5_GZ", "E14_5_CP", "E16_5_GZ", "E16_5_CP", "E18_5_GZ", "E18_5_CP"))
Idents(seurat_object) <- "time_Region"

seurat_object <- PrepSCTFindMarkers(seurat_object)

E14_5_GZ_CP <- FindMarkers(seurat_object, ident.1 = "E14_5_CP", ident.2 = "E14_5_GZ", min.pct = 0.1)
E16_5_GZ_CP <- FindMarkers(seurat_object, ident.1 = "E16_5_CP", ident.2 = "E16_5_GZ", min.pct = 0.1)
E18_5_GZ_CP <- FindMarkers(seurat_object, ident.1 = "E18_5_CP", ident.2 = "E18_5_GZ", min.pct = 0.1)

E14_5_E16_5_CP <- FindMarkers(seurat_object, ident.1 = "E16_5_CP", ident.2 = "E14_5_CP", min.pct = 0.1)
E14_5_E18_5_CP <- FindMarkers(seurat_object, ident.1 = "E18_5_CP", ident.2 = "E14_5_CP", min.pct = 0.1)
E16_5_E18_5_CP <- FindMarkers(seurat_object, ident.1 = "E18_5_CP", ident.2 = "E16_5_CP", min.pct = 0.1)

E14_5_E16_5_GZ <- FindMarkers(seurat_object, ident.1 = "E16_5_GZ", ident.2 = "E14_5_GZ", min.pct = 0.1)
E14_5_E18_5_GZ <- FindMarkers(seurat_object, ident.1 = "E18_5_GZ", ident.2 = "E14_5_GZ", min.pct = 0.1)
E16_5_E18_5_GZ <- FindMarkers(seurat_object, ident.1 = "E18_5_GZ", ident.2 = "E16_5_GZ", min.pct = 0.1)

E14_5_GZ_CP <- E14_5_GZ_CP %>% filter(abs(avg_log2FC) >= 1, p_val_adj < 0.05)
E16_5_GZ_CP <- E16_5_GZ_CP %>% filter(abs(avg_log2FC) >= 1, p_val_adj < 0.05)
E18_5_GZ_CP <- E18_5_GZ_CP %>% filter(abs(avg_log2FC) >= 1, p_val_adj < 0.05)

E14_5_E16_5_CP <- E14_5_E16_5_CP %>% filter(abs(avg_log2FC) >= 1, p_val_adj < 0.05)
E14_5_E18_5_CP <- E14_5_E18_5_CP %>% filter(abs(avg_log2FC) >= 1, p_val_adj < 0.05)
E16_5_E18_5_CP <- E16_5_E18_5_CP %>% filter(abs(avg_log2FC) >= 1, p_val_adj < 0.05)

E14_5_E16_5_GZ <- E14_5_E16_5_GZ %>% filter(abs(avg_log2FC) >= 1, p_val_adj < 0.05)
E14_5_E18_5_GZ <- E14_5_E18_5_GZ %>% filter(abs(avg_log2FC) >= 1, p_val_adj < 0.05)
E16_5_E18_5_GZ <- E16_5_E18_5_GZ %>% filter(abs(avg_log2FC) >= 1, p_val_adj < 0.05)

features <- unique(c(rownames(E14_5_GZ_CP), rownames(E16_5_GZ_CP), rownames(E18_5_GZ_CP), rownames(E14_5_E16_5_CP), rownames(E14_5_E18_5_CP), rownames(E16_5_E18_5_CP), rownames(E14_5_E16_5_GZ), rownames(E14_5_E18_5_GZ), rownames(E16_5_E18_5_GZ)))
features <- features[-grep("Rik$", features)]
features <- features[-grep("Gm[0-9]*", features)]

ht <- SCP::GroupHeatmap(seurat_object, features = features, slot = "counts", exp_method = "zscore", group.by = "Region", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$Region), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palette = "viridis", nlabel = 0, show_row_names = F, row_names_side = "left", n_split = 6, use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot
pdf("GZ_CP_DE_Clusters.pdf", width = 6, height = 8)
ht$plot
dev.off()
write.csv(ht$feature_metadata, "GZ_CP_DE_Clusters.csv", quote = F, row.names = F)
de_markers <- ht$feature_metadata
dir.create("Cortex_Region_features")
for(i in unique(de_markers$feature_split)){
        gene_list <- as.character(de_markers[which(de_markers$feature_split == i),'features'])
        enrichRes <- enrichAnalysis(genelist = gene_list, geneType = "SYMBOL", species = "mmu", database = c("go", "kegg", "MSigDb"), GO.model = "BP", MSigDb.signature = "H", sampleName = paste0("enrich_", i), minGSSize = 5, outpath = "Cortex_Region_features")
}

seurat_object$time_Region <- factor(seurat_object$time_Region, levels = c("E14_5_GZ", "E16_5_GZ", "E18_5_GZ", "E14_5_CP", "E16_5_CP", "E18_5_CP"))
RNA_Matrix <- AverageExpression(seurat_object, assays = "SCT", features = features, group.by = "time_Region", slot = "counts", return.seurat = T)

library(ClusterGVis)

data_matrix <- GetAssayData(RNA_Matrix, slot = "data")
colnames(data_matrix) <- levels(seurat_object$time_Region)

cm <- clusterData(exp = expm1(data_matrix), scaleData = TRUE, cluster.method = "kmeans", cluster.num = 8) # kmeans
visCluster(object = cm, ht.col.list = list(col_range = seq(-2, 2, length = 200), col_color = paletteContinuous(set = "solarExtra", n = 200)), plot.type = "both", cluster_rows = TRUE, column_names_rot = 45, use_raster = TRUE, show_row_dend = F)
pdf("RNA_heatmap.pdf", width = 8, height = 8, onefile = F)
visCluster(object = CGVs, ht.col.list = list(col_range = seq(-2, 2, length = 200), col_color = paletteContinuous(set = "solarExtra", n = 200)), plot.type = "both", cluster_rows = TRUE, column_names_rot = 45, use_raster = TRUE, cluster.order =  1:7)
dev.off()
#!------------------------------------------------------------------------------------------------------------------
seurat_object <- subset(tmp_object, cells = colnames(tmp_object)[which(tmp_object$cell_identity %in% c("LGE_MGE_AP", "LGE_AP", "LGE_IPC", "LGE_InN", "CPU_pre_MSN", "CPU_MSN"))])
seurat_object <- shiny_st(seurat = seurat_object, isVisium = F, assay = "SCT", image = "E14_5")
seurat_object <- subset(seurat_object, subset = orig.ident != "removed_cells")
seurat_object$cell_identity[which(seurat_object$cell_identity == "LGE_MGE_AP")] <- "LGE_AP"
seurat_object$cell_identity <- factor(seurat_object$cell_identity, levels = intersect(levels(seurat_object$cell_identity), unique(seurat_object$cell_identity)))

temp_object <- CreateSeuratObject(GetAssayData(seurat_object, assay = "RNA", layer = "counts"), min.cells = 1, min.features = 1)
temp_object@meta.data <- seurat_object@meta.data
temp_object[["spatial_umap"]] <- seurat_object[["spatial_umap"]]
temp_object@images <- seurat_object@images
seurat_object <- temp_object

seurat_object <- Split_Layers(seurat_object, split.by = "orig.ident")
seurat_object <- SCTransform(seurat_object, vst.flavor = "v2", method = "glmGamPoi", vars.to.regress = c("S.Score", "G2M.Score"))
seurat_object <- RunPCA(seurat_object, npcs = 30, verbose = FALSE)
# one-liner to run Integration
seurat_object <- IntegrateLayers(object = seurat_object, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = 'harmony', normalization.method = "SCT", verbose = FALSE)
seurat_object <- RunUMAP(seurat_object, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony", return.model = T)
seurat_object <- FindNeighbors(seurat_object, dims = 1:30, reduction = 'harmony', annoy.metric = "cosine", graph.name = c('harmony_nn', 'harmony_snn'))
seurat_object <- FindClusters(seurat_object, resolution = seq(0.1, 1, 0.1), graph.name = 'harmony_snn')
wrap_plots(purrr::map(seq(0.1, 1, 0.1), function(x) DimPlot(seurat_object, reduction = "umap_harmony", group.by = paste0("harmony_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)
seurat_object$seurat_clusters <- seurat_object$harmony_snn_res.0.4
Idents(seurat_object) <- "cell_identity"

DefaultAssay(seurat_object) <- "SCT"
seurat_object <- JoinLayers(seurat_object, assay = "RNA")

library(CytoTRACE2)
seurat_object <- cytotrace2(seurat_object, is_seurat = TRUE, slot_type = "counts", species = 'mouse')

p1 <- plot_cytotrace2(seurat_object, reduction = "umap_harmony")

(p1$CytoTRACE2_UMAP + p1$CytoTRACE2_Potency_UMAP + p1$CytoTRACE2_Relative_UMAP) & SetAxes()
ggsave("Development_QC_cell_identity/MSN_CytoTRACE2.pdf")

qsave(seurat_object, "output/brain_development_ST_final_tutorial_MSN.qs")
seurat2scanpy(seurat_object, assay = "SCT", manual_color = manual_colors, major_umap = "umap_harmony", savefile = paste0("Trajectory_analysis/anndata/brain_development_ST_final_MSN.h5ad"))
#!----------------------------------------------------------------------------------------------------------------------------------------------------------
# deal Adult15 spatial RNA
tmp_object <- seurat_list$Adult15

tmp_object <- SCTransform(tmp_object, vst.flavor = "v2", method = "glmGamPoi", vars.to.regress = c("nFeature_RNA", "S.Score", "G2M.Score"), verbose = TRUE)
tmp_object <- RunPCA(tmp_object, npcs = 30)
tmp_object <- RunUMAP(tmp_object, reduction = "pca", dims = 1:30, min.dist = 0.1, reduction.name = "SCT_umap", return.model = T)
tmp_object <- FindNeighbors(tmp_object, reduction = "pca", dims = 1:30, annoy.metric = "cosine", graph.name = c('SCT_nn', 'SCT_snn'))
tmp_object <- FindClusters(tmp_object, cluster.name = "SCT_cluster", resolution = 1.5, graph.name = "SCT_snn")
tmp_object <- FindClusters(tmp_object, resolution = seq(0.1, 2, 0.1), graph.name = "SCT_snn")

wrap_plots(map(seq(0.6, 2, 0.1), function(x) DimPlot(tmp_object, reduction = "SCT_umap", group.by = paste0("SCT_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)

tmp_object <- FindClusters(tmp_object, cluster.name = "SCT_cluster", resolution = 0.8, graph.name = "SCT_snn")
DimPlot(tmp_object, reduction = "SCT_umap", group.by = "SCT_cluster", cols = selected_colors, label = T) & NoLegend()
Seurat::SpatialPlot(tmp_object, shape = 22, cols = selected_colors, group.by = "SCT_cluster") & NoLegend()

Idents(tmp_object) <- "SCT_cluster"
selected_cells <- CellsByIdentities(tmp_object)
Seurat::SpatialPlot(tmp_object, cells.highlight = selected_cells[setdiff(names(selected_cells), "NA")], cols.highlight = c("#FFFF00", alpha("lightgrey", 0)), crop = FALSE, shape = 22, facet.highlight = T, combine = T, ncol = 7, alpha = NULL) & NoLegend()

library(Banksy)
locs <- Seurat::GetTissueCoordinates(tmp_object)[,seq_len(2)]
colnames(locs) <- c("array_row", "array_col")
tmp_object@meta.data <- cbind(tmp_object@meta.data, locs)
tmp_object <- RunBanksy(tmp_object, lambda = 0.2, assay = 'SCT', slot = 'scale.data', dimx = "array_row", dimy = "array_col", features = 'variable', k_geom = 24, assay_name = "banksy_SCT")
npcs = 30
tmp_object <- RunPCA(tmp_object, assay = "banksy_SCT", reduction.name = "pca_banksy", npcs = npcs, features = rownames(tmp_object))
tmp_object <- RunUMAP(tmp_object, reduction = "pca_banksy", dims = 1:npcs, min.dist = 0.3, reduction.name = "banksy_umap", return.model = T, n.neighbors = 24)
tmp_object <- FindNeighbors(tmp_object, reduction = "pca_banksy", dims = 1:npcs, annoy.metric = "cosine", graph.name = c('banksy_nn', 'banksy_snn'), k.param = 24)

tmp_object <- FindClusters(tmp_object, resolution = seq(0.1, 2, 0.1), graph.name = "banksy_snn")
wrap_plots(map(seq(0.6, 2, 0.1), function(x) DimPlot(tmp_object, reduction = "banksy_umap", group.by = paste0("banksy_snn_res.", x), cols = selected_colors, label = T) + NoLegend()), ncol = 5)

tmp_object <- FindClusters(tmp_object, cluster.name = "banksy_cluster", resolution = 1.9, graph.name = "banksy_snn")
tmp_object$banksy_cluster <- as.character(tmp_object$banksy_cluster)
tmp_object$banksy_cluster[which(tmp_object$banksy_snn_res.2 == 5)] <- 19
tmp_object$banksy_cluster <- factor(tmp_object$banksy_cluster, levels = 0:19)

DimPlot(tmp_object, label = T, cols = selected_colors, group.by = "banksy_cluster", reduction = "banksy_umap") + NoLegend()
Seurat::SpatialPlot(tmp_object, cols = selected_colors, group.by = "banksy_cluster", crop = FALSE, shape = 22)
Idents(tmp_object) <- "banksy_cluster"
selected_cells <- CellsByIdentities(tmp_object)
Seurat::SpatialPlot(tmp_object, cells.highlight = selected_cells[setdiff(names(selected_cells), "NA")], cols.highlight = c("#FFFF00", alpha("lightgrey", 0)), crop = FALSE, shape = 22, facet.highlight = T, combine = T, ncol = 7, alpha = NULL) & NoLegend()

predictions <- table(tmp_object$banksy_cluster, tmp_object$SCT_cluster)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)

ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "SCT_cluster",
    low = "#ffffc8", high = "#7d0025") + xlab("banksy_cluster") + ylab("") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

tmp_object$cell_identity <- as.character(tmp_object$SCT_cluster)
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(19))] <- "ITL1" # "Camk2a", "Slc1a2", "Glul"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(1))] <- "ITL23GL" # Excitatory neurons, Cortex Layer 2/3 "Rasgrf2", "Lrrtm4", "Cux2", "Pdzrn3", "Slc9a9", "Unc5d", | "Calb1", "Cux2", "Lamp5"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(2))] <- "ITL4GL" # Excitatory spiny neuron "Rorb", "Myo16", "Kcnh5", "Astn2", "Zmat4", "Unc5d", "Brinp3" | "Satb2", "Rorb", "Ovol2", "Rspo1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(7))] <- "ITL5GL" # "Slc24a3", "Plcxd2", "Fat3" | "Satb2", "Rorb"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(13))] <- "PTGL" # Excitatory neurons, Cortex PT, Cortex Layer 5b "Parm1", "Pcsk5" | "Fezf2", "Kcng1", "Npsr1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(0))] <- "CTGL_ITL6GL" # Excitatory neurons, Cortex CT "Zfpm2", "Pde1a", "Nos1ap", "Cdh18", "Neto1", "Ncald", "Thsd7b", | "Foxp2", "Satb2", "Sulf1", "Garnl3"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(18))] <- "OPC" # "Lhfpl3", "Vcan" | "Pdgfra", "Olig2", "C1ql1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(6))] <- "ODC1" # "Plp1", "Mag", "Prr5l", "Dock10", "Trf" | "Mal", "Mog", "Opalin", "Ppp1r14a"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(5))] <- "ODC2" # "Mbp", "Ccdc106", "Mobp", "Fht1", "Tlk1", "Prkch", "Plekhb1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(10))] <- "CPU_MSN" # "Gng7", "Rgs9", "Pde10a", "Adcy5", "Penk", "Drd2", "Rarb", "Gpr88", "Ppp1r1b", "Dach1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(11))] <- "Astrocytes" # "Mertk", "Gli2", "Plpp3", "Prex2", "Lsamp", "Slc39a12", "Nwd1", "Zbtb20", "Ptprz1", "Nhsl1" | "Fam107a", "S1pr1", "Plpp3", "Gpc5"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(14))] <- "SSTGA" # MGE-derived neurogliaform cells, Sst positive "Grin3a", "Synpr", "Cdh13", "Grip1", "Reln", "Sox6", "Sst", "Npas3" | "Elfn1", "Lypd6b", "Sst"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(12))] <- "PVGA" # "Erbb4", "Cntnap4", "Gad1", "Kcnc1", "Grip1", "Eya4", "Syt2", "Pvalb", "Zfp536" | "Pvalb", "Erbb4"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(15))] <- "VIPGA" # "Adarb2", "Grip1", "Erbb4", "Zfp536", "Adra1a", "Dner" | "Vip", "Rgs12", "Npas1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(17))] <- "VLMC" # Vascular and leptomeningeal cells "Ptgds", "Ahnak", "Apod" | "Slc6a13", "Tbx18"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(9))] <- "Endo" # "Slco1a4", "Atp10a", "Flt1", "Adgrl4", "Mecom", "Ptprb"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(16))] <- "MG" # "Inpp5d", "Apbb1ip", "Ctss", "Zfhx3", "Tmcc3" | "C1qb", "P2ry12", "Siglech"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(3))] <- "unknown" # "Calm1"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(4))] <- "unknown" # "Vwa3a"
tmp_object$cell_identity[which(tmp_object$cell_identity %in% c(8))] <- "unknown"

tmp_object$cell_identity <- factor(tmp_object$cell_identity, levels = c("ITL1", "ITL23GL", "ITL4GL", "ITL5GL", "PTGL", "CTGL_ITL6GL", "OPC", "ODC1", "ODC2", "CPU_MSN", "Astrocytes", "SSTGA", "PVGA", "VIPGA", "VLMC", "Endo", "MG", "unknown"))
Seurat::SpatialPlot(tmp_object, group.by = "cell_identity", shape = 22, cols = cluster_colors)
DimPlot(tmp_object, label = T, group.by = "cell_identity", cols = cluster_colors) + NoLegend() + SetAxes()

tmp_object$Cluster <- as.character(tmp_object$SCT_cluster)
tmp_object$Cluster[which(tmp_object$Cluster %in% c(19))] <- "C1"
tmp_object$Cluster[which(tmp_object$Cluster %in% c(1))] <- "C2"
tmp_object$Cluster[which(tmp_object$Cluster %in% c(2))] <- "C3"
tmp_object$Cluster[which(tmp_object$Cluster %in% c(7))] <- "C4"
tmp_object$Cluster[which(tmp_object$Cluster %in% c(13))] <- "C5"
tmp_object$Cluster[which(tmp_object$Cluster %in% c(0))] <- "C6"
tmp_object$Cluster[which(tmp_object$Cluster %in% c(18))] <- "C7"
tmp_object$Cluster[which(tmp_object$Cluster %in% c(6))] <- "C8"
tmp_object$Cluster[which(tmp_object$Cluster %in% c(5))] <- "C9"
tmp_object$Cluster[which(tmp_object$Cluster %in% c(10))] <- "C10"
tmp_object$Cluster[which(tmp_object$Cluster %in% c(11))] <- "C11"
tmp_object$Cluster[which(tmp_object$Cluster %in% c(14))] <- "C12"
tmp_object$Cluster[which(tmp_object$Cluster %in% c(12))] <- "C13"
tmp_object$Cluster[which(tmp_object$Cluster %in% c(15))] <- "C14"
tmp_object$Cluster[which(tmp_object$Cluster %in% c(17))] <- "C15"
tmp_object$Cluster[which(tmp_object$Cluster %in% c(9))] <- "C16"
tmp_object$Cluster[which(tmp_object$Cluster %in% c(16))] <- "C17"
tmp_object$Cluster[which(tmp_object$Cluster %in% c(3, 4, 8))] <- "C18"

tmp_object$Cluster <- factor(tmp_object$Cluster, levels = paste0("C", 1:18))
Cluster_colors <- cluster_colors[match(levels(tmp_object$cell_identity), names(cluster_colors))]
names(Cluster_colors) <- levels(tmp_object$Cluster)
Seurat::SpatialPlot(tmp_object, group.by = "Cluster", shape = 22, cols = Cluster_colors)

tmp_object$banksy_domains <- as.character(tmp_object$banksy_cluster)
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c(13))] <- "ITL1"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c(10))] <- "ITL23GL"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c(4, 8))] <- "ITL23GL"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c(3))] <- "ITL4GL"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c(2, 19))] <- "ITL5GL"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c(0))] <- "PTGL"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c(1))] <- "CTGL_ITL6GL"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c(5))] <- "CTGL_ITL6GL"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c(9))] <- "ODC1"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c(7))] <- "ODC2"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c(6))] <- "CPU_MSN"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c(11))] <- "Astrocytes"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c(16))] <- "SSTGA"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c(14))] <- "PVGA"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c(17))] <- "VIPGA"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c(12))] <- "Endo"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c(15))] <- "unknown"
tmp_object$banksy_domains[which(tmp_object$banksy_domains %in% c(18))] <- "unknown"
tmp_object$banksy_domains <- factor(tmp_object$banksy_domains, levels = c("ITL1", "ITL23GL", "ITL4GL", "ITL5GL", "PTGL", "CTGL_ITL6GL", "ODC1", "ODC2", "CPU_MSN", "Astrocytes", "SSTGA", "PVGA", "VIPGA", "Endo", "unknown"))

tmp_object$banksy_Cluster <- as.character(tmp_object$banksy_cluster)
tmp_object$banksy_Cluster[which(tmp_object$banksy_Cluster %in% c(13))] <- "D1"
tmp_object$banksy_Cluster[which(tmp_object$banksy_Cluster %in% c(4, 8, 10))] <- "D2"
tmp_object$banksy_Cluster[which(tmp_object$banksy_Cluster %in% c(3))] <- "D3"
tmp_object$banksy_Cluster[which(tmp_object$banksy_Cluster %in% c(2, 19))] <- "D4"
tmp_object$banksy_Cluster[which(tmp_object$banksy_Cluster %in% c(0))] <- "D5"
tmp_object$banksy_Cluster[which(tmp_object$banksy_Cluster %in% c(1, 5))] <- "D6"
tmp_object$banksy_Cluster[which(tmp_object$banksy_Cluster %in% c(9))] <- "D7"
tmp_object$banksy_Cluster[which(tmp_object$banksy_Cluster %in% c(7))] <- "D8"
tmp_object$banksy_Cluster[which(tmp_object$banksy_Cluster %in% c(6))] <- "D9"
tmp_object$banksy_Cluster[which(tmp_object$banksy_Cluster %in% c(11))] <- "D10"
tmp_object$banksy_Cluster[which(tmp_object$banksy_Cluster %in% c(16))] <- "D11"
tmp_object$banksy_Cluster[which(tmp_object$banksy_Cluster %in% c(14))] <- "D12"
tmp_object$banksy_Cluster[which(tmp_object$banksy_Cluster %in% c(17))] <- "D13"
tmp_object$banksy_Cluster[which(tmp_object$banksy_Cluster %in% c(12))] <- "D14"
tmp_object$banksy_Cluster[which(tmp_object$banksy_Cluster %in% c(15, 18))] <- "D15"

tmp_object$banksy_Cluster <- factor(tmp_object$banksy_Cluster, levels = paste0("D", 1:15))

banksy_Cluster_colors <- cluster_colors[match(levels(tmp_object$banksy_domains), names(cluster_colors))]
names(banksy_Cluster_colors) <- levels(tmp_object$banksy_Cluster)
Seurat::SpatialPlot(tmp_object, group.by = "banksy_Cluster", shape = 22, cols = banksy_Cluster_colors, crop = FALSE)
DimPlot(tmp_object, label = T, reduction = "banksy_umap", group.by = "banksy_Cluster", cols = banksy_Cluster_colors) + NoLegend() + SetAxes()
DimPlot(tmp_object, label = T, reduction = "banksy_umap", group.by = "banksy_domains", cols = cluster_colors) + NoLegend() + SetAxes()
#!-----------------------------------------------------------------
temp_object <- readRDS("output/FC_mouse_snRNA_2023_NSMB.rds")
temp_object$cell_identity <- factor(temp_object$cell_identity, levels = c("OBGL", "ITL23GL", "ITL45GL", "ITL5GL", "NPGL", "PTGL", "CTGL", "ITL6GL", "OPC", "OGC", "CLAGL", "D12MSN", "ASC", "OBGA", "SSTGA", "PVGA", "VIPGA", "STRGA", "VLMC", "MGL"))

anchors <- FindTransferAnchors(reference = temp_object, query = tmp_object, dims = 1:30, reference.reduction = "harmony", reference.assay = "RNA", query.assay = "RNA", reduction = "pcaproject")

seurat_object <- MapQuery(anchorset = anchors, reference = temp_object, query = seurat_object, refdata = list(celltype = "cell_identity"), reference.reduction = "harmony", new.reduction.name = "ref.harmony", reduction.model = "harmony_umap", projectumap.args = list(reduction.name = "ref.harmony.umap"))
colnames(seurat_object@meta.data)[204] <- "predicted.celltype_2023_NSMB.score"
colnames(seurat_object@meta.data)[205] <- "predicted.celltype_2023_NSMB"

tmp_object$predicted.celltype_2023_NSMB <- factor(tmp_object$predicted.celltype_2023_NSMB, levels = intersect(levels(temp_object$cell_identity), unique(tmp_object$predicted.celltype_2023_NSMB)))

DefaultAssay(tmp_object) <- "SCT"
Idents(tmp_object) <- "cell_identity"

Seurat::SpatialPlot(tmp_object, shape = 22, cols = ClusterIDColors, group.by = "Cluster")
Seurat::SpatialPlot(tmp_object, shape = 22, cols = domainIDColors, group.by = "banksy_cluster")
Seurat::SpatialPlot(tmp_object, shape = 22, cols = cluster_colors)
DimPlot(tmp_object, cols = ClusterIDColors, group.by = "Cluster", label = T, repel = T) + SetAxes()
DimPlot(tmp_object, cols = cluster_colors, label = T, repel = T) + SetAxes()
DimPlot(tmp_object, reduction = "banksy_umap", cols = domainIDColors, group.by = "banksy_cluster", label = T, repel = T) + SetAxes()

library(ggsankey)

meta_data <- tmp_object@meta.data %>% make_long(Cluster, banksy_Cluster)
meta_data$node <- factor(meta_data$node, levels = c(levels(tmp_object$Cluster), levels(tmp_object$banksy_Cluster)))
meta_data$next_node <- factor(meta_data$next_node, levels = c(levels(tmp_object$Cluster), levels(tmp_object$banksy_Cluster)))

ggplot(meta_data, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) + geom_sankey(color = "black", lwd = 0.2, flow.alpha = 0.5, show.legend = FALSE, na.rm = FALSE) + scale_x_discrete(expand = c(0, 0.8)) + 
		scale_fill_manual(values = c(Cluster_colors, banksy_Cluster_colors), drop = FALSE) + geom_sankey_text(size = 3.2, color = "black") + theme_void() + theme(axis.text.x = element_text())
ggsave("cell2location_Res/Adult15_SCT_spot_cluster_banksy_domian_sankey_plot.pdf")
#!--------------------------------------------------------------------------------
library(UCell)
source("/home/yiyelinfeng/scripts/Rscripts/lung_project/IPF/spatial-lung-fibrosis/scripts/custom_colors.R")

selected_clusters <- colnames(tmp_object@meta.data)[63:82]

tmp_meta <- as.matrix(tmp_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(tmp_object$banksy_Cluster, "_", tmp_object$banksy_domains)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)
tmp_meta0 <- tmp_meta0[paste0(levels(tmp_object$banksy_Cluster), "_", levels(tmp_object$banksy_domains)), ]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_2023_NSMB_", names(seurat_list)[5], "_banksy_zscore.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_2023_NSMB_", names(seurat_list)[5], "_banksy_zscore1.pdf"), width = 10, height = 10)

tmp_meta <- as.matrix(tmp_object@meta.data[,selected_clusters[2:8]])
rownames(tmp_meta) <- paste0(tmp_object$banksy_Cluster, "_", tmp_object$banksy_domains)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta <- tmp_meta[grep("GL$", rownames(tmp_meta)),]
tmp_meta0 <- apply(tmp_meta, 2, function(x){scale(x)})
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)
tmp_meta0 <- tmp_meta0[intersect(paste0(levels(tmp_object$banksy_Cluster), "_", levels(tmp_object$banksy_domains)), rownames(tmp_meta0)),]

pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_2023_NSMB_", names(seurat_list)[5], "_banksy_zscore_selected0.pdf"), width = 10, height = 10)
pheatmap::pheatmap(t(tmp_meta0), cellwidth = 20, cellheight = 20, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = F, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_2023_NSMB_", names(seurat_list)[5], "_banksy_zscore1_selected0.pdf"), width = 10, height = 10)

selected_clusters <- colnames(tmp_object@meta.data)[85:205]
tmp_meta <- as.matrix(tmp_object@meta.data[,selected_clusters])
rownames(tmp_meta) <- paste0(tmp_object$banksy_Cluster, "_", tmp_object$banksy_domains)
tmp_meta <- limma::avereps(tmp_meta)
tmp_meta <- apply(tmp_meta, 2, function(x){scale(x)})
rownames(tmp_meta) <- unique(paste0(tmp_object$banksy_Cluster, "_", tmp_object$banksy_domains))
tmp_meta0 <- normalize_to_01(as.vector(tmp_meta))
tmp_meta0 <- matrix(tmp_meta0, nrow = nrow(tmp_meta))
rownames(tmp_meta0) <- rownames(tmp_meta)
colnames(tmp_meta0) <- colnames(tmp_meta)

pheatmap::pheatmap(tmp_meta0, cellwidth = 15, cellheight = 15, show_colnames = T, show_rownames = T, cluster_cols = T, cluster_rows = T, color = col_scale_mako_custom, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_2023_Nature_", names(seurat_list)[5], "_banksy_zscore.pdf"), width = 30, height = 10)
pheatmap::pheatmap(tmp_meta0, cellwidth = 15, cellheight = 15, show_colnames = T, show_rownames = T, cluster_cols = T, cluster_rows = T, color = col_scale_bupu, filename = paste0("Development_QC_cell_identity/cell_identity_c2l_density_2023_Nature_", names(seurat_list)[5], "_banksy_zscore1.pdf"), width = 30, height = 10)

Idents(tmp_object) <- "banksy_domains"

key_markers <- c("Calm1","Camk2a", "Slc1a2", "Lrrtm4", "Rasgrf2", "Slc9a9", "Cux2", "Unc5d", "Astn2", "Rorb", "Myo16", "Brinp3", "Satb2", "Slc24a3", "Fat3", "Zfpm2", "Pde1a", "Nos1ap", "Cdh18", "Neto1", "Ncald", "Garnl3", "Plp1", "Mag", "Prr5l", "Dock10", "Mbp", "Mobp", "Fht1", "Tlk1", "Pde10a", "Gng7", "Adcy5", "Penk", "Drd2", "Lsamp", "Mertk", "Plpp3", "Prex2", "Grip1", "Cdh13", "Grin3a", "Reln", "Sst", "Erbb4", "Gad1", "Pvalb", "Syt2", "Adarb2", "Adra1a", "Vip", "Slco1a4", "Atp10a", "Flt1")

DotPlot(tmp_object, group.by = "banksy_domains", features = key_markers, cols = c("lightgrey", "red"), col.min = 0.1, dot.min = 0.1) + coord_flip() + RotatedAxis() + labs(y = "banksy domains", x = "Features") + NoLegend()
ggsave("Development_QC_cell_identity/Adult_cell_identity_key_markers.pdf")
#!-------------------------------------------------------------------------------------------------------------------------
library(semla)
tmp_object <- qread("output/brain_development_ST_final_tutorial_cortex.qs")
seurat_list <- qread("output/seurat_list_SCT.qs")
seurat_list$E14_5 <- RenameCells(seurat_list$E14_5, add.cell.id = "E14_5")
seurat_list$E16_5 <- RenameCells(seurat_list$E16_5, add.cell.id = "E16_5")
seurat_list$E18_5 <- RenameCells(seurat_list$E18_5, add.cell.id = "E18_5")

seurat_object <- UpdateSeuratForSemla(seurat_list$E14_5)
seurat_object <- LoadImages(seurat_object, image_height = 1080)
temp <- seurat_object@meta.data[,c("row", "col")]
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
MapMultipleFeatures(seurat_object, image_use = "raw", shape = "tile", pt_size = 6, max_cutoff = 0.95, min_cutoff = 0.05, override_plot_dims = TRUE,  colors = c("#017351", "#332288", "#44AA99", "#C9DB74", "#DDCC77", "#CC6677","#AA4499", "#16F2F2"), features = c("ITL23GL", "ITL45GL", "ITL5GL", "PTGL", "CTGL", "ITL6GL", "D12MSN", "ASC"))
VlnPlot(seurat_object, cols = c("#017351", "#332288", "#44AA99", "#C9DB74", "#DDCC77", "#CC6677","#AA4499", "#16F2F2"), features = c("ITL23GL", "ITL45GL", "ITL5GL", "PTGL", "CTGL", "ITL6GL", "D12MSN", "ASC"), pt.size = 0, stack = T, flip = T) + ylab("cell2location score")

MapMultipleFeatures(seurat_object, image_use = "raw", shape = "tile", pt_size = 6, max_cutoff = 0.95, min_cutoff = 0.05, override_plot_dims = TRUE,  colors = c("#017351", "#332288", "#44AA99", "#117733", "#DDCC77"), features = c("L2/3_IT_ENT_Glut", "L5_ET_CTX_Glut", "L4/5_IT_CTX_Glut", "L6_CT_CTX_Glut", "L6b_CTX_Glut"))
VlnPlot(seurat_object, cols = c("#017351", "#332288", "#44AA99", "#117733", "#DDCC77"), features = c("L2/3_IT_ENT_Glut", "L5_ET_CTX_Glut", "L4/5_IT_CTX_Glut", "L6_CT_CTX_Glut", "L6b_CTX_Glut"), pt.size = 0, stack = T, flip = T) + ylab("cell2location score")
#!--------------------------------------------------------------------------------------------------------------------------
source("~/scripts/seurat2scanpy/shiny_st.R")
options(browser = "/usr/bin/firefox")
seurat_object$Region <- seurat_object$banksy_domains
seurat_object <- shiny_st(seurat = seurat_object, isVisium = F, assay = "SCT", image = "Adult15")
seurat_object$Region <- factor(seurat_object$Region, levels = c("Layer_I", "Layer_II_III", "Layer_IV", "Layer_V", "Layer_VI", "CC", "CPU"))
seurat_object$banksy_type <- "excitory_neuron"
seurat_object$banksy_type[which(seurat_object$banksy_domains %in% c("SSTGA", "SSTGA", "PVGA", "VIPGA"))] <- "inhibitory_neuron"
seurat_object$banksy_type[which(seurat_object$banksy_domains %in% c("ODC1", "ODC2"))] <- "ODC"
seurat_object$banksy_type[which(seurat_object$banksy_domains %in% c("Astrocytes"))] <- "Astrocytes"
seurat_object$banksy_type[which(seurat_object$banksy_domains %in% c("CPU_MSN"))] <- "MSN"
seurat_object$banksy_type[which(seurat_object$banksy_domains %in% c("Endo", "unknown"))] <- "others"
seurat_object$banksy_type <- factor(seurat_object$banksy_type, levels = rev(c("excitory_neuron", "inhibitory_neuron", "ODC", "Astrocytes", "MSN", "others")))

ggplot(seurat_object@meta.data[which(seurat_object$banksy_domains != "unknown" & seurat_object$Region %in% c("Layer_I", "Layer_II_III", "Layer_IV", "Layer_V", "Layer_VI")),], aes(Region, fill = banksy_type)) + geom_bar(position = "fill") + 
	ylab("percent of cells") + xlab("") + 
    theme(plot.background = element_rect(fill = NA), legend.key = element_rect(size = 2), legend.background = element_rect(color = NA), axis.text.y = element_text(size = 10, face='bold'), 
    axis.text.x = element_text(size = 10, face='bold', angle = 45, hjust = 1, vjust = 1)) +  
    scale_fill_manual(values = rev(c("#009E73", "#E7298A", "#0072B2", "#16F2F2", "#D38B5C", "#886C00")))
ggsave("Development_QC_cell_identity/Adult15_Region_banksy_domains_barplot.pdf")

tmp_data <- table(seurat_object$Region, seurat_object$banksy_domains)
tmp_data <- tmp_data/rowSums(tmp_data)
meta_data <- data.frame(celltype = rep(colnames(tmp_data), each = nrow(tmp_data)), values = as.vector(tmp_data), group = rep(rownames(tmp_data), ncol(tmp_data)))
meta_data$celltype <- factor(meta_data$celltype, levels = colnames(tmp_data))
meta_data$group <- factor(meta_data$group, levels = rownames(tmp_data))
meta_data <- meta_data[which(meta_data$celltype %in% c("ITL1", "ITL23GL", "ITL4GL", "ITL5GL", "PTGL", "CTGL_ITL6GL", "ODC1", "Astrocytes", "Endo") & meta_data$group %in% c("Layer_I", "Layer_II_III", "Layer_IV", "Layer_V", "Layer_VI")),]

ggplot(meta_data, aes(x = group, y = values, fill = celltype)) + geom_bar(stat = "identity", position = "dodge") + ylab("percent of cells") + xlab("") + 
	coord_cartesian(ylim = c(0, 0.1)) +
	geom_text(aes(label = paste0(round(meta_data$values * 100, 1), "%")), position = position_dodge(0.9), vjust = -0.25) +
    theme(plot.background = element_rect(fill = NA), legend.position = "right", legend.key = element_rect(size = 2), legend.background = element_rect(color = NA), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12)) +
    scale_fill_manual(values = cluster_colors)
ggplot(meta_data, aes(x = group, y = values, fill = celltype)) + geom_bar(stat = "identity", position = "dodge") + ylab("percent of cells") + xlab("") + 
	coord_cartesian(ylim = c(0.15, 0.8)) +
	geom_text(aes(label = paste0(round(meta_data$values * 100, 1), "%")), position = position_dodge(0.9), vjust = -0.25) +
    theme(plot.background = element_rect(fill = NA), legend.position = "right", legend.key = element_rect(size = 2), legend.background = element_rect(color = NA), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12)) +
    scale_fill_manual(values = cluster_colors)
ggsave("Development_QC_cell_identity/Adult15_cortex_Region_banksy_domains_barplot.pdf")

seurat_object <- subset(seurat_object, subset = banksy_domains %in% c("ITL23GL", "ITL4GL", "ITL5GL", "PTGL", "CTGL_ITL6GL", "ODC1", "Astrocytes", "SSTGA", "PVGA", "VIPGA"))
seurat_object <- shiny_st(seurat = seurat_object, isVisium = F, assay = "SCT", image = "Adult15")
seurat_object <- subset(seurat_object, subset = banksy_domains != "removed_cells")
DefaultAssay(seurat_object) <- "RNA"
seurat_object@assays <- seurat_object@assays[1]

tmp_object <- qread("output/brain_development_ST_final_tutorial_cortex.qs")
DefaultAssay(tmp_object) <- "RNA"
tmp_object@assays <- tmp_object@assays[1]

seurat_object <- RenameCells(seurat_object, add.cell.id = "Adult15")
seurat_object@meta.data <- seurat_object@meta.data[,c(1:18, 59:62)]
seurat_object$cell_identity <- seurat_object$banksy_domains
seurat_object$Cluster <- seurat_object$banksy_Cluster
temp_object <- CreateSeuratObject(GetAssayData(seurat_object, assay = "RNA", layer = "counts"), min.cells = 1, min.features = 1)

temp_object <- merge(tmp_object, seurat_object)
temp_object <- JoinLayers(temp_object)

tmp <- seurat_object@meta.data[,c("array_col", "array_row")]
tmp[,1] <- tmp[,1] + 5400 + 1092
tmp[,2] <- 1080 - tmp[,2]
tmp <- rbind(tmp_object@reductions$spatial_umap@cell.embeddings, as.matrix(tmp))
all(rownames(tmp) == colnames(temp_object))
temp_object[["spatial_umap"]] <- CreateDimReducObject(embeddings = tmp)
seurat_object <- temp_object
seurat_object <- Split_Layers(seurat_object, split.by = "orig.ident")
seurat_object <- SCTransform(seurat_object, vst.flavor = "v2", method = "glmGamPoi", vars.to.regress = c("nFeature_RNA", "S.Score", "G2M.Score"))
seurat_object <- RunPCA(seurat_object, npcs = 30, verbose = FALSE)
# one-liner to run Integration
seurat_object <- IntegrateLayers(object = seurat_object, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = 'harmony', normalization.method = "SCT", verbose = FALSE)
seurat_object <- RunUMAP(seurat_object, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony", return.model = T)

library(Banksy)
library(SeuratWrappers)
library(harmony)
locs <- c()
for(i in names(seurat_object@images)){
	locs <- rbind(locs, Seurat::GetTissueCoordinates(seurat_object, image = i)[,seq_len(2)])
}
all(rownames(locs) == colnames(seurat_object))
colnames(locs) <- c("array_row", "array_col")
seurat_object$array_row <- locs$array_row
seurat_object$array_col <- locs$array_col
seurat_object$orig.ident <- factor(seurat_object$orig.ident, levels = unique(seurat_object$orig.ident))
plot(FetchData(seurat_object, c('array_row', 'array_col')), col = factor(seurat_object$orig.ident))
seurat_object <- RunBanksy(seurat_object, lambda = 0.2, assay = 'SCT', dimx = "array_row", dimy = "array_col", features = 'variable', k_geom = 24, assay_name = "banksy_SCT")
npcs = 30
seurat_object <- RunPCA(seurat_object, assay = "banksy_SCT", reduction.name = "pca_banksy", npcs = npcs, features = rownames(seurat_object))
seurat_object <- RunHarmony(seurat_object, group.by.vars = "orig.ident", reduction.name = 'pca_banksy', reduction.save = 'banksy_harmony')
seurat_object <- RunUMAP(seurat_object, reduction = "banksy_harmony", dims = 1:npcs, min.dist = 0.3, reduction.name = "banksy_umap", return.model = T, n.neighbors = 24)
#!--------------------------------------------------------------------------------------
seurat_object$cell_identity <- factor(seurat_object$cell_identity, levels = c(levels(tmp_object$cell_identity), "ITL23GL", "ITL4GL", "ITL5GL", "PTGL", "CTGL_ITL6GL"))
seurat_object$Cluster <- factor(seurat_object$Cluster, levels = c(levels(tmp_object$Cluster), "D2", "D3", "D4", "D5", "D6"))

Idents(seurat_object) <- "cell_identity"
DefaultAssay(seurat_object) <- "SCT"
seurat_object <- JoinLayers(seurat_object, assay = "RNA")

qsave(seurat_object, "output/brain_ST_final_tutorial_cortex.qs")
seurat2scanpy(seurat_object, assay = "SCT", manual_color = manual_colors, major_umap = "banksy_umap", savefile = paste0("Trajectory_analysis/anndata/brain_ST_final_cortex.h5ad"))
#!--------------------------------------------------------------------------------------------------------------------------
seurat_object <- qread("output/brain_ST_final_tutorial_cortex.qs")
seurat_object <- subset(seurat_object, subset = orig.ident != "Adult15")
DefaultAssay(seurat_object) <- "SCT"
seurat_object$orig.ident <- factor(seurat_object$orig.ident, levels = unique(seurat_object$orig.ident))
seurat_object$cell_type <- factor(seurat_object$cell_type, levels = levels(seurat_object$cell_type)[1:7])
seurat_object <- PrepSCTFindMarkers(seurat_object)

Idents(seurat_object) <- "cell_type"
all_markers <- FindAllMarkers(seurat_object, min.pct = 0.1, only.pos = T)
de_markers <- all_markers %>% filter(avg_log2FC > 0.5, pct.1 > 0.3, p_val_adj < 0.05)
write.csv(de_markers, "Cortex_cell_type_DE_markers.csv", quote = F, row.names = F)
de_markers$index <- 1:nrow(de_markers)
de_markers <- de_markers[order(de_markers$gene, de_markers$avg_log2FC, decreasing = TRUE),]
de_markers <- de_markers[!duplicated(de_markers$gene),]
de_markers <- de_markers[order(de_markers$index),]
de_markers$final_cluster <- de_markers$cluster[match(de_markers$gene, de_markers$gene)]


key_markers <- c("Veph1", "Dct","Pax6", "Eomes", "Ltbp1", "Slc17a6", "Lhx2", "Neurod2", "Ptger3", "Abca8a", "Tmem132d", "Pappa2", "Nr4a3", "Fezf2", "Tafa1", "Slc26a7", "Rasgrf2", "Tafa2", "Inhba", "Sema3c", "Unc5d", "Adamts2", "Hs3st4", "Pdzd2", "Tle4", "Dpp10", "Tmem200a", "Pou6f2", "Nr4a2", "Cux1", "Cux2", "Frem2", "Smad3", "Tnc", "Megf10", "Lrp4")

ht <- SCP::GroupHeatmap(seurat_object, features = de_markers$gene, slot = "counts", exp_method = "zscore", group.by = "cell_type", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$cell_type), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, feature_split = de_markers$final_cluster, feature_split_palcolor = list(cluster_colors[match(levels(seurat_object$cell_type), names(cluster_colors))]), heatmap_palette = "viridis", nlabel = 0, show_row_names = FALSE, row_names_side = "left", features_label = key_markers, use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot
pdf("Cortex_cell_type_DE_features.pdf", width = 12, height = 12)
ht$plot
dev.off()
#!-------------------------------------------------------
# identity conserved markers cross time
conserved_markers <- list()
for(i in levels(seurat_object$cell_type)){
	conserved_markers[[i]] <- FindConservedMarkers(seurat_object, ident.1 = i, grouping.var = "orig.ident", logfc.threshold = 0.25, min.pct = 0.1, min.diff.pct = 0.1, only.pos = T, verbose = FALSE, assay = "SCT")
}

for(i in levels(seurat_object$cell_type)){
	if(length(grep("E14_5_", colnames(conserved_markers[[i]]))) == 0){
		conserved_markers[[i]]$E14_5_p_val <- NA
		conserved_markers[[i]]$E14_5_avg_log2FC <- NA
		conserved_markers[[i]]$E14_5_pct.1 <- NA
		conserved_markers[[i]]$E14_5_pct.2 <- NA
		conserved_markers[[i]]$E14_5_p_val_adj <- 0
	}
	if(length(grep("E16_5_", colnames(conserved_markers[[i]]))) == 0){
		conserved_markers[[i]]$E16_5_p_val <- NA
		conserved_markers[[i]]$E16_5_avg_log2FC <- NA
		conserved_markers[[i]]$E16_5_pct.1 <- NA
		conserved_markers[[i]]$E16_5_pct.2 <- NA
		conserved_markers[[i]]$E16_5_p_val_adj <- 0
	}
	if(length(grep("E18_5_", colnames(conserved_markers[[i]]))) == 0){
		conserved_markers[[i]]$E18_5_p_val <- NA
		conserved_markers[[i]]$E18_5_avg_log2FC <- NA
		conserved_markers[[i]]$E18_5_pct.1 <- NA
		conserved_markers[[i]]$E18_5_pct.2 <- NA
		conserved_markers[[i]]$E18_5_p_val_adj <- 0
	}
	conserved_markers[[i]] <- conserved_markers[[i]][,c(paste0(rep(c("E14_5_", "E16_5_", "E18_5_"), each = 5), rep(c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"), 3)))]
	conserved_markers[[i]]$gene <- rownames(conserved_markers[[i]])
}

all_conserved_markers <- c()
for(i in levels(seurat_object$cell_type)){
	conserved_markers[[i]]$cluster <- i
	all_conserved_markers <- rbind(all_conserved_markers, conserved_markers[[i]])
}
for(i in 1:nrow(all_conserved_markers)){
	all_conserved_markers$max_pval_adj[i] <- max(all_conserved_markers[i,c("E14_5_p_val_adj", "E16_5_p_val_adj", "E18_5_p_val_adj")])
}
all_conserved_markers <- all_conserved_markers[which(all_conserved_markers$max_pval_adj < 0.05),]
de_conserved_markers <- all_conserved_markers[which(all_conserved_markers$E14_5_avg_log2FC > 0.5 | all_conserved_markers$E16_5_avg_log2FC > 0.5 | all_conserved_markers$E18_5_avg_log2FC > 0.5),]

all_conserved_markers <- FindAllMarkers(seurat_object, features = unique(de_conserved_markers$gene))
all_conserved_markers$index <- 1:nrow(all_conserved_markers)
all_conserved_markers <- all_conserved_markers[order(all_conserved_markers$gene, all_conserved_markers$avg_log2FC, decreasing = TRUE),]
all_conserved_markers <- all_conserved_markers[!duplicated(all_conserved_markers$gene),]
all_conserved_markers <- all_conserved_markers[order(all_conserved_markers$index),]
de_conserved_markers$final_cluster <- all_conserved_markers$cluster[match(de_conserved_markers$gene, all_conserved_markers$gene)]
write.csv(de_conserved_markers, "Cortex_de_conserved_markers_cell_type.csv", quote = F, row.names = F)

key_markers <- c("Veph1", "Dct","Pax6", "Eomes", "Ltbp1", "Slc17a6", "Lhx2", "Neurod2", "Ptger3", "Abca8a", "Tmem132d", "Pappa2", "Nr4a3", "Fezf2", "Tafa1", "Slc26a7", "Rasgrf2", "Tafa2", "Inhba", "Sema3c", "Unc5d", "Adamts2", "Hs3st4", "Pdzd2", "Tle4", "Dpp10", "Tmem200a", "Pou6f2", "Nr4a2", "Cux1", "Cux2", "Frem2", "Smad3", "Tnc", "Megf10", "Lrp4")

ht <- SCP::GroupHeatmap(seurat_object, features = de_conserved_markers$gene, slot = "counts", exp_method = "zscore", group.by = "cell_type", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$cell_type), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, feature_split = de_conserved_markers$cluster, feature_split_palcolor = list(cluster_colors[match(levels(seurat_object$cell_type), names(cluster_colors))]), heatmap_palette = "viridis", nlabel = 0, show_row_names = FALSE, row_names_side = "left", features_label = key_markers, use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot
pdf("Cortex_cell_type_conserved_features.pdf", width = 12, height = 12)
ht$plot
dev.off()
#!AP DE analysis-----------------------------------------------------------------------------------
seurat_object$time_cell_type <- paste0(seurat_object$orig.ident, "_", seurat_object$cell_type)
seurat_object$time_cell_type <- factor(seurat_object$time_cell_type, levels = intersect(paste0(rep(levels(seurat_object$orig.ident), length(levels(seurat_object$cell_type))), "_", rep(levels(seurat_object$cell_type), each = 3)), unique(seurat_object$time_cell_type)))
Idents(seurat_object) <- "time_cell_type"

E14_AP_markers <- FindMarkers(seurat_object, ident.1 = "E14_5_AP", ident.2 = c("E16_5_AP", "E18_5_AP"), min.pct = 0.1, only.pos = T)
E16_AP_markers <- FindMarkers(seurat_object, ident.1 = "E16_5_AP", ident.2 = c("E14_5_AP", "E18_5_AP"), min.pct = 0.1, only.pos = T)
E18_AP_markers <- FindMarkers(seurat_object, ident.1 = "E18_5_AP", ident.2 = c("E14_5_AP", "E16_5_AP"), min.pct = 0.1, only.pos = T)
E14_AP_markers$gene <- rownames(E14_AP_markers)
E14_AP_markers$cluster <- "E14_5_AP"
E16_AP_markers$gene <- rownames(E16_AP_markers)
E16_AP_markers$cluster <- "E16_5_AP"
E18_AP_markers$gene <- rownames(E18_AP_markers)
E18_AP_markers$cluster <- "E18_5_AP"

all_markers <- rbind(E14_AP_markers, E16_AP_markers, E18_AP_markers)
de_markers <- all_markers %>% filter(avg_log2FC > 0.5, pct.1 > 0.3, p_val_adj < 0.05)
write.csv(de_markers, "AP_cell_type_DE_markers.csv", quote = F, row.names = F)

cell_type_markers <- read.csv("Cortex_cell_type_DE_markers.csv")
de_markers <- de_markers[which(de_markers$gene %in% c(de_conserved_markers$gene[which(de_conserved_markers$final_cluster == "AP")], cell_type_markers$gene[which(cell_type_markers$cluster == "AP")])),]

de_markers <- de_markers[-grep("mt-", de_markers$gene),]
de_markers <- de_markers[-grep("Rik$", de_markers$gene),]
write.csv(de_markers, "AP_cell_type_DE_markers_filtered.csv", quote = F, row.names = F)

ht <- SCP::GroupHeatmap(seurat_object, features = unique(de_markers$gene), slot = "counts", exp_method = "zscore", group.by = "cell_type", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$cell_type), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, n_split = 5, heatmap_palette = "viridis", nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot
pdf("Cortex_AP_features.pdf", width = 8, height = 10)
ht$plot
dev.off()
write.csv(ht$feature_metadata, "Cortex_AP_features.csv", quote = F)
de_markers <- ht$feature_metadata
dir.create("Cortex_AP_features")
for(i in unique(de_markers$feature_split)){
        gene_list <- as.character(de_markers[which(de_markers$feature_split == i),'features'])
        enrichRes <- enrichAnalysis(genelist = gene_list, geneType = "SYMBOL", species = "mmu", database = c("go", "kegg", "MSigDb"), GO.model = "BP", MSigDb.signature = "H", sampleName = paste0("enrich_", i), minGSSize = 5, outpath = "Cortex_AP_features")
}
#! moscot analysis------------------------------------------------------------------------------------------------------------------------------------------
seurat_object <- qread("output/brain_ST_final_tutorial_cortex.qs")

transition_E14_E16 <- as.matrix(read.csv("Trajectory_analysis/moscot/tutorials/figures_all_cortex/E14_5_to_E16_5_transition_matrix_moscot_TemporalProblem.csv", row.names = 1))
transition_E14_E16[which(transition_E14_E16 <= 0.15)] <- 0
transition_E16_E18 <- as.matrix(read.csv("Trajectory_analysis/moscot/tutorials/figures_all_cortex/E16_5_to_E18_5_transition_matrix_moscot_TemporalProblem.csv", row.names = 1))
transition_E16_E18[which(transition_E16_E18 <= 0.15)] <- 0
transition_E16_E18["Layer_6b", "CThPN"] <- 0
transition_E16_E18["SCPN", "CThPN"] <- 0
transition_E18_P56 <- as.matrix(read.csv("Trajectory_analysis/moscot/tutorials/figures_all_cortex/E18_5_to_P56_transition_matrix_moscot_TemporalProblem.csv", row.names = 1))
transition_E18_P56[which(transition_E18_P56 <= 0.15)] <- 0

meta_data <- table(seurat_object$cell_type, seurat_object$orig.ident)
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
meta_data$node <- factor(meta_data$node, levels = levels(seurat_object$cell_type))
meta_data$next_node <- factor(meta_data$next_node, levels = levels(seurat_object$cell_type))

ggplot(meta_data, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) + geom_sankey(color = "black", lwd = 0.1, flow.alpha = 0.5, show.legend = FALSE, na.rm = FALSE) + scale_x_discrete(expand = c(0, 0.8)) + 
		scale_fill_manual(values = cluster_colors, drop = FALSE) + 
		geom_sankey_text(size = 3, color = "black") + theme_void() + theme(axis.text.x = element_text())
ggsave("Trajectory_analysis/moscot/tutorials/figures_all_cortex/cortex_transition_matrix_moscot_TemporalProblem_sankeyplot.pdf")

meta_data <- read.csv("Trajectory_analysis/moscot/tutorials/figures_all_cortex/brain_development_ST_final_all_cortex_moscot_cellrank.csv", row.names = 1)
seurat_object@meta.data <- cbind(seurat_object@meta.data, meta_data[,32:119])
Seurat::SpatialPlot(seurat_object, shape = 22, crop = F, features = "AP_push_E14_5_to_E16_5", images = c("E14_5", "E16_5")) & scale_fill_gradientn(colours = viridis::viridis(15, direction = -1))

TF_file_list <- list.files("Trajectory_analysis/moscot/tutorials/figures_all_cortex_new/")
TF_file_list <- TF_file_list[grep("*push.csv", TF_file_list)]
TF_file_list <- gsub(TF_file_list, pattern = ".csv", replacement = "")
core_list <- gsub(TF_file_list, pattern = "_push", replacement = "")
core_list <- gsub(core_list, pattern = "_drivers_tf", replacement = "")
core_genes <- list()
for(i in 1:length(TF_file_list)){
	tmp_data <- read.csv(paste0("Trajectory_analysis/moscot/tutorials/figures_all_cortex_new/", TF_file_list[i], ".csv"))
	colnames(tmp_data) <- c("genes", "corr", "pval", "qval", "ci_low", "ci_high")
	tmp_data$FDR <- log10(tmp_data$qval)
	tmp_data <- tmp_data[-grep("[0-9]Rik$", tmp_data$genes),]
	#core_genes[[i]] <- tmp_data$genes[1:30]
	core_genes[[core_list[i]]] <- tmp_data$genes[which(tmp_data$corr >= 0.2)]
	#core_gene <- unique(c(core_gene, core_genes[[i]]))

	#p0 <- ggplot(tmp_data[which(tmp_data$corr > 0 & tmp_data$qval < 0.05),], aes(FDR, corr)) + 
	#	geom_point(aes(size = corr, color = qval)) + scale_color_continuous(low = "#e06663", high = "#327eba", name = 'p.adjust', guide = guide_colorbar(reverse = TRUE)) +
	#	scale_size_continuous(range = c(1, 3)) + theme_classic(base_size = 10) + xlab("log10(qval)") + ggtitle(i) +
	#	theme(panel.grid = element_blank(), axis.text = element_text(size = 10, color = "black", face = "bold"), axis.title = element_text(size = 12, color = "black", face = "bold")) + 
	#	geom_text_repel(data = filter(tmp_data, genes %in% core_gene), max.overlaps = getOption("ggrepel.max.overlaps", default = 15), aes(label = genes), size = 4, color = 'black')
	#print(p0)
	#ggsave(paste0("Trajectory_analysis/moscot/tutorials/figures_all_cortex/", i, ".pdf"))
}
for(i in 1:length(core_genes)){
	label_size = 10
	if(length(core_genes[[i]]) < 100){
		height = 10
		n = length(core_genes[[i]])
		label_size = 3
	}else if(length(core_genes[[i]]) < 200){
		height = 12
		n = 20
	}else if(length(core_genes[[i]]) < 400){
		height = 14
		n = 20
	}else{
		height = 16
		n = 20
	}
	ht <- SCP::GroupHeatmap(seurat_object, features = core_genes[[i]], slot = "counts", exp_method = "zscore", group.by = "cell_type", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$cell_type), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palette = "viridis", nlabel = 0, show_row_names = FALSE, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE), features_label = head(core_genes[[i]], n = n), label_size = label_size)

	pdf(paste0("Cortex_moscot_", names(core_genes)[i], ".pdf"), width = 8, height = 12)
		print(ht$plot)
	dev.off()
}

core_genes <- core_genes[grep("_E14_5_E16_5_E18_5|Layer_6b|UL_CPN", names(core_genes))]
for(i in names(core_genes)){
	core_genes[[i]] <- core_genes[[i]][1:50]
}
names(core_genes) <- gsub(names(core_genes), pattern = "_E14_5_E16_5_E18_5|_E16_5_E18_5", replacement = "")

seurat_object <- UCell::AddModuleScore_UCell(seurat_object, features = core_genes, ncores = 1, name = "")
seurat_object$AP_normalized <- normalize_to_01(seurat_object$AP)
seurat_object$IPC_normalized <- normalize_to_01(seurat_object$IPC)
seurat_object$MigN_normalized <- normalize_to_01(seurat_object$MigN)
seurat_object$Layer_6b_normalized <- normalize_to_01(seurat_object$Layer_6b)
seurat_object$CThPN_normalized <- normalize_to_01(seurat_object$CThPN)
seurat_object$SCPN_normalized <- normalize_to_01(seurat_object$SCPN)
seurat_object$UL_CPN_normalized <- normalize_to_01(seurat_object$UL_CPN)

VlnPlot(seurat_object, features = paste0(names(core_genes), "_normalized"), group.by = "cell_type", pt.size = 0, split.by = "orig.ident", cols = c("#34D916", "#00D4E6", "#1E90FF")) & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.3, position = position_dodge(0.9), show.legend = F)
VlnPlot(seurat_object, features = "MigN_normalized", group.by = "cell_type", pt.size = 0, split.by = "orig.ident", cols = c("#34D916", "#00D4E6", "#1E90FF")) & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.3, position = position_dodge(0.9), show.legend = F)

marker_genes <- core_genes
hic_object <- qread("output/brain_development_Spatial_HiC_final_tutorial_cortex_new.qs")
DefaultAssay(hic_object) <- "scAB250k"
for(i in names(marker_genes)){
	marker_genes[[i]] <- intersect(rownames(hic_object), unique(hic_object@misc$gene_infor_250k$gene_region[which(hic_object@misc$gene_infor_250k$gene_name %in% marker_genes[[i]])]))
	hic_object@meta.data[,i] <- colMeans(as.matrix(GetAssayData(hic_object, assay = "scAB250k", slot = "data")[marker_genes[[i]],]))
}

hic_object$AP_normalized <- normalize_to_01(hic_object$AP)
hic_object$IPC_normalized <- normalize_to_01(hic_object$IPC)
hic_object$MigN_normalized <- normalize_to_01(hic_object$MigN)
hic_object$Layer_6b_normalized <- normalize_to_01(hic_object$Layer_6b)
hic_object$CThPN_normalized <- normalize_to_01(hic_object$CThPN)
hic_object$SCPN_normalized <- normalize_to_01(hic_object$SCPN)
hic_object$UL_CPN_normalized <- normalize_to_01(hic_object$UL_CPN)

VlnPlot(hic_object, features = paste0(names(core_genes), "_normalized"), pt.size = 0, group.by = "orig.ident", cols = c("#34D916", "#00D4E6", "#1E90FF"), ncol = 7) & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.3) & NoLegend()
VlnPlot(hic_object, features = paste0(names(core_genes), "_normalized"), pt.size = 0, group.by = "banksy_enrich", split.by = "orig.ident", cols = c("#34D916", "#00D4E6", "#1E90FF")) & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.3, position = position_dodge(0.9)) & NoLegend()
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

Idents(seurat_object) <- "cell_type"
dir.create("Trajectory_analysis/moscot/tutorials/figures_all_cortex/tmp/selected_TFs")
pdf("Trajectory_analysis/moscot/tutorials/figures_all_cortex/selected_TFs/selected_TFs_VlnPlot.pdf", width = 20)
for(i in core_genes){
	p0 <- 	VlnPlot(seurat_object, features = i, pt.size = 0, split.by = "orig.ident")
	print(p0)
}
dev.off()

DefaultAssay(seurat_object) <- "spARC_SCT"
for(i in core_genes$core_genes){
	try({
		p0 <- Seurat::SpatialPlot(seurat_object, features = i, shape = 22, crop = F) & scale_fill_gradientn(colours = paletteContinuous(set = "sambaNight"))
		p1 <- Seurat::SpatialPlot(hic_object, shape = 22, crop = F, features = hic_object@misc$gene_infor_500k$gene_region[which(hic_object@misc$gene_infor_500k$gene_name == i)]) & scale_fill_gradientn(values = seq(0, 1, 0.01), colours = paletteContinuous(set = "blueYellow", n = 100))
		print(p0/p1)
		ggsave(paste0("Trajectory_analysis/moscot/tutorials/figures_all_cortex/tmp/selected_TFs/selected_TFs_", i, ".png"), width = 24, height = 24)
	})
}
#!-------------------------------------------------------
meta_data <- read.csv("Trajectory_analysis/moscot/tutorials/figures_all_cortex/brain_development_ST_final_all_cortex_moscot_cellrank.csv", row.names = 1)
meta_data <- meta_data[colnames(seurat_object),]
meta_data <- meta_data[,c(32:36, 44:50)]
seurat_object@meta.data <- cbind(seurat_object@meta.data, meta_data)
#! E14_5 AP push--------------------------------------------------------------------------
Seurat::SpatialPlot(seurat_object, shape = 22, crop = F, features = "AP_mapping_14_5_to_16_5", images = c("E16_5")) & scale_fill_gradientn(colours = viridis::viridis(15, direction = -1))

seurat_object$E14_AP_push <- NA
seurat_object$E14_AP_push[which(seurat_object$cell_type == "AP" & seurat_object$orig.ident == "E14_5")] <- "source"
seurat_object$E14_AP_push[which(seurat_object$AP_mapping_14_5_to_16_5 > 0.0001 & seurat_object$orig.ident == "E16_5"  & seurat_object$cell_type == "AP")] <- "AP"
seurat_object$E14_AP_push[which(seurat_object$AP_mapping_14_5_to_16_5 > 0.0001 & seurat_object$orig.ident == "E16_5"  & seurat_object$cell_type == "IPC")] <- "IPC"
seurat_object$E14_AP_push[which(seurat_object$AP_mapping_14_5_to_16_5 > 0.0001 & seurat_object$orig.ident == "E16_5"  & seurat_object$cell_type == "CThPN")] <- "CThPN"
seurat_object$E14_AP_push[which(seurat_object$AP_mapping_14_5_to_16_5 > 0.0001 & seurat_object$orig.ident == "E16_5"  & seurat_object$cell_type == "SCPN")] <- "SCPN"
seurat_object$E14_AP_push <- factor(seurat_object$E14_AP_push, levels = c("AP", "IPC", "CThPN", "SCPN", "source"))
Seurat::SpatialPlot(seurat_object, shape = 22, crop = F, group.by = "E14_AP_push", images = c("E16_5"), cols = cluster_colors)

E14_AP_push <- list()
for(i in c("AP", "IPC", "CThPN", "SCPN")){
	seurat_object$tmp <- as.character(seurat_object$E14_AP_push)
	seurat_object$tmp[which(seurat_object$tmp == "source")] <- i
	Idents(seurat_object) <- "tmp"
	E14_AP_push[[i]] <- FindConservedMarkers(seurat_object, ident.1 = i, grouping.var = "orig.ident", logfc.threshold = 0.25, min.pct = 0.1, min.diff.pct = 0.1, only.pos = T, verbose = FALSE, assay = "SCT")
}
E14_AP_push_conserved_markers <- c()
for(i in c("AP", "IPC", "CThPN", "SCPN")){
	E14_AP_push[[i]]$gene <- rownames(E14_AP_push[[i]])
	E14_AP_push[[i]]$cluster <- i
	E14_AP_push_conserved_markers <- rbind(E14_AP_push_conserved_markers, E14_AP_push[[i]])
}
E14_AP_push_conserved_markers$max_pval_adj <- NA
for(i in 1:nrow(E14_AP_push_conserved_markers)){
	E14_AP_push_conserved_markers$max_pval_adj[i] <- max(E14_AP_push_conserved_markers[i,c("E14_5_p_val_adj", "E16_5_p_val_adj")])
}

E14_AP_push_conserved_markers <- E14_AP_push_conserved_markers[which(E14_AP_push_conserved_markers$max_pval_adj < 0.05),]
E14_AP_push_conserved_markers <- E14_AP_push_conserved_markers[which(E14_AP_push_conserved_markers$E14_5_avg_log2FC > 0.5 | E14_AP_push_conserved_markers$E16_5_avg_log2FC > 0.5),]
write.csv(E14_AP_push_conserved_markers, "E14_AP_push_conserved_markers.csv", quote = F, row.names = F)

ht <- SCP::GroupHeatmap(seurat_object, features = E14_AP_push_conserved_markers$gene[which(E14_AP_push_conserved_markers$cluster == "SCPN")], slot = "data", exp_method = "zscore", group.by = "cell_type", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$cell_type), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palette = "viridis", nlabel = 0, show_row_names = T, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot
Idents(seurat_object) <- "cell_type"
VlnPlot(seurat_object, features = E14_AP_push_conserved_markers$gene[which(E14_AP_push_conserved_markers$cluster == "CThPN")], pt.size = 0, split.by = "orig.ident", cols = c("#34D916", "#00D4E6", "#1E90FF"))
#! E14_5 IPC push--------------------------------------------------------------------------
Seurat::SpatialPlot(seurat_object, shape = 22, crop = F, features = "IPC_mapping_14_5_to_16_5", images = c("E16_5")) & scale_fill_gradientn(colours = viridis::viridis(15, direction = -1))
seurat_object$E14_IPC_push <- NA
seurat_object$E14_IPC_push[which(seurat_object$cell_type == "IPC" & seurat_object$orig.ident == "E14_5")] <- "source"
seurat_object$E14_IPC_push[which(seurat_object$IPC_mapping_14_5_to_16_5 > 0.0001 & seurat_object$orig.ident == "E16_5"  & seurat_object$cell_type == "IPC")] <- "IPC"
seurat_object$E14_IPC_push[which(seurat_object$IPC_mapping_14_5_to_16_5 > 0.0001 & seurat_object$orig.ident == "E16_5"  & seurat_object$cell_type == "MigN")] <- "MigN"
seurat_object$E14_IPC_push <- factor(seurat_object$E14_IPC_push, levels = c("IPC", "MigN", "source"))
Seurat::SpatialPlot(seurat_object, shape = 22, crop = F, group.by = "E14_IPC_push", images = c("E16_5"), cols = cluster_colors)

E14_IPC_push <- list()
for(i in c("IPC", "MigN")){
	seurat_object$tmp <- as.character(seurat_object$E14_IPC_push)
	seurat_object$tmp[which(seurat_object$tmp == "source")] <- i
	Idents(seurat_object) <- "tmp"
	E14_IPC_push[[i]] <- FindConservedMarkers(seurat_object, ident.1 = i, grouping.var = "orig.ident", logfc.threshold = 0.25, min.pct = 0.1, min.diff.pct = 0.1, only.pos = T, verbose = FALSE, assay = "SCT")
}
E14_IPC_push_conserved_markers <- c()
for(i in c("IPC", "MigN")){
	E14_IPC_push[[i]]$gene <- rownames(E14_IPC_push[[i]])
	E14_IPC_push[[i]]$cluster <- i
	E14_IPC_push_conserved_markers <- rbind(E14_IPC_push_conserved_markers, E14_IPC_push[[i]])
}
E14_IPC_push_conserved_markers$max_pval_adj <- NA
for(i in 1:nrow(E14_IPC_push_conserved_markers)){
	E14_IPC_push_conserved_markers$max_pval_adj[i] <- max(E14_IPC_push_conserved_markers[i,c("E14_5_p_val_adj", "E16_5_p_val_adj")])
}

E14_IPC_push_conserved_markers <- E14_IPC_push_conserved_markers[which(E14_IPC_push_conserved_markers$max_pval_adj < 0.05),]
E14_IPC_push_conserved_markers <- E14_IPC_push_conserved_markers[which(E14_IPC_push_conserved_markers$E14_5_avg_log2FC > 0.5 | E14_IPC_push_conserved_markers$E16_5_avg_log2FC > 0.5),]
write.csv(E14_IPC_push_conserved_markers, "E14_IPC_push_conserved_markers.csv", quote = F, row.names = F)

ht <- SCP::GroupHeatmap(seurat_object, features = E14_IPC_push_conserved_markers$gene[which(E14_IPC_push_conserved_markers$cluster == "MigN")], slot = "data", exp_method = "zscore", group.by = "cell_type", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$cell_type), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palette = "viridis", nlabel = 0, show_row_names = T, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot
Idents(seurat_object) <- "cell_type"
VlnPlot(seurat_object, features = E14_IPC_push_conserved_markers$gene[which(E14_IPC_push_conserved_markers$cluster == "MigN")], pt.size = 0, split.by = "orig.ident", cols = c("#34D916", "#00D4E6", "#1E90FF"))
#! E14_5 MigN push--------------------------------------------------------------------------
Seurat::SpatialPlot(seurat_object, shape = 22, crop = F, features = "MigN_mapping_14_5_to_16_5", images = c("E16_5")) & scale_fill_gradientn(colours = viridis::viridis(15, direction = -1))
seurat_object$E14_MigN_push <- NA
seurat_object$E14_MigN_push[which(seurat_object$cell_type == "MigN" & seurat_object$orig.ident == "E14_5")] <- "source"
seurat_object$E14_MigN_push[which(seurat_object$MigN_mapping_14_5_to_16_5 > 0.0001 & seurat_object$orig.ident == "E16_5"  & seurat_object$cell_type == "MigN")] <- "MigN"
seurat_object$E14_MigN_push[which(seurat_object$MigN_mapping_14_5_to_16_5 > 0.0001 & seurat_object$orig.ident == "E16_5"  & seurat_object$cell_type == "SCPN")] <- "SCPN"
seurat_object$E14_MigN_push[which(seurat_object$MigN_mapping_14_5_to_16_5 > 0.0001 & seurat_object$orig.ident == "E16_5"  & seurat_object$cell_type == "UL_CPN")] <- "UL_CPN"
seurat_object$E14_MigN_push <- factor(seurat_object$E14_MigN_push, levels = c("MigN", "SCPN", "UL_CPN", "source"))
Seurat::SpatialPlot(seurat_object, shape = 22, crop = F, group.by = "E14_MigN_push", images = c("E16_5"), cols = cluster_colors)

E14_MigN_push <- list()
for(i in c("MigN", "SCPN", "UL_CPN")){
	seurat_object$tmp <- as.character(seurat_object$E14_MigN_push)
	seurat_object$tmp[which(seurat_object$tmp == "source")] <- i
	Idents(seurat_object) <- "tmp"
	E14_MigN_push[[i]] <- FindConservedMarkers(seurat_object, ident.1 = i, grouping.var = "orig.ident", logfc.threshold = 0.25, min.pct = 0.1, min.diff.pct = 0.1, only.pos = T, verbose = FALSE, assay = "SCT")
}
E14_MigN_push_conserved_markers <- c()
for(i in c("MigN", "SCPN", "UL_CPN")){
	E14_MigN_push[[i]]$gene <- rownames(E14_MigN_push[[i]])
	E14_MigN_push[[i]]$cluster <- i
	E14_MigN_push_conserved_markers <- rbind(E14_MigN_push_conserved_markers, E14_MigN_push[[i]])
}
E14_MigN_push_conserved_markers$max_pval_adj <- NA
for(i in 1:nrow(E14_MigN_push_conserved_markers)){
	E14_MigN_push_conserved_markers$max_pval_adj[i] <- max(E14_MigN_push_conserved_markers[i,c("E14_5_p_val_adj", "E16_5_p_val_adj")])
}

E14_MigN_push_conserved_markers <- E14_MigN_push_conserved_markers[which(E14_MigN_push_conserved_markers$max_pval_adj < 0.05),]
E14_MigN_push_conserved_markers <- E14_MigN_push_conserved_markers[which(E14_MigN_push_conserved_markers$E14_5_avg_log2FC > 0.5 | E14_MigN_push_conserved_markers$E16_5_avg_log2FC > 0.5),]
write.csv(E14_MigN_push_conserved_markers, "E14_MigN_push_conserved_markers.csv", quote = F, row.names = F)

ht <- SCP::GroupHeatmap(seurat_object, features = E14_MigN_push_conserved_markers$gene[which(E14_MigN_push_conserved_markers$cluster == "UL_CPN")], slot = "data", exp_method = "zscore", group.by = "cell_type", split.by = "orig.ident", group_palcolor = list(cluster_colors[match(levels(seurat_object$cell_type), names(cluster_colors))]), cell_split_palcolor = list(c("#34D916", "#00D4E6", "#1E90FF")), cluster_rows = TRUE, cluster_columns = FALSE, heatmap_palette = "viridis", nlabel = 0, show_row_names = T, row_names_side = "left", use_raster = TRUE, raster_by_magick = TRUE, ht_params = list(show_row_dend = FALSE))
ht$plot
Idents(seurat_object) <- "cell_type"
VlnPlot(seurat_object, features = E14_MigN_push_conserved_markers$gene[which(E14_MigN_push_conserved_markers$cluster == "UL_CPN")], pt.size = 0, split.by = "orig.ident", cols = c("#34D916", "#00D4E6", "#1E90FF"))
#!--------------------------------------------------------------------------------------------------------------------------
seurat_object$cell_identity <- as.character(seurat_object$seurat_clusters)
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(16))] <- "VZ_NPC"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(2))] <- "Cortex_Layer_VI_ExN"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(13))] <- "Cortex_Layer_V_ExN"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(5))] <- "Cortex_Layer_IV_Spiny_ExN" # Excitatory spiny neurons
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(7))] <- "Cortex_Layer_II_III_ExN"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(11))] <- "Cortex_CR"

seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(9))] <- "Cortex_SST_InN" # somatostatin
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(12))] <- "Cortex_Pvalb_InN"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(14))] <- "Cortex_PV_InN" # parvalbumin
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(15))] <- "Cortex_Vip_InN" # vasoactive intestinal peptide
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(21))] <- "Cortex_ODC"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(1))] <- "CC_ODC" # Corpus Callosum oligodendrocytes
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(6))] <- "ODC"

seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(3))] <- "CPU_MSN"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(4))] <- "CPU_MSN"

seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(19))] <- "ChP_cilia"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(20))] <- "Ependymal"

seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(8))] <- "Endo"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(10))] <- "MG"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(18))] <- "Ery"

seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(0))] <- "unknown"
seurat_object$cell_identity[which(seurat_object$cell_identity %in% c(17))] <- "unknown"

seurat_object$cell_identity <- factor(seurat_object$cell_identity, levels = c("VZ_NPC", "Cortex_Layer_VI_ExN", "Cortex_Layer_V_ExN", "Cortex_Layer_IV_Spiny_ExN", "Cortex_Layer_II_III_ExN", "Cortex_CR", "Cortex_SST_InN", "Cortex_Pvalb_InN", "Cortex_PV_InN", "Cortex_Vip_InN", "Cortex_ODC", "CC_ODC", "ODC", "CPU_MSN", "ChP_cilia", "Ependymal", "Endo", "MG", "Ery", "unknown"))

seurat_object$Region <- as.character(seurat_object$Cluster)
seurat_object$Region[which(seurat_object$Region %in% paste0("C", c(1)))] <- "VZ"
seurat_object$Region[which(seurat_object$Region %in% paste0("C", c(2:11)))] <- "Cortex"

seurat_object$Region[which(seurat_object$Region %in% paste0("C", c(12)))] <- "CC"

seurat_object$Region[which(seurat_object$Region %in% paste0("C", c(14)))] <- "CPU"
seurat_object$Region[which(seurat_object$Region %in% paste0("C", c(15)))] <- "ChP"
seurat_object$Region[which(seurat_object$Region %in% paste0("C", c(16)))] <- "Ependymal"

seurat_object$Region[which(seurat_object$Region %in% paste0("C", c(13, 17:20)))] <- "others"

seurat_object$Region <- factor(seurat_object$Region, levels = c("VZ", "Cortex", "CC", "CPU", "ChP", "Ependymal", "others"))

seurat_object$cell_identity_allnames <- as.character(seurat_object$cell_identity)
seurat_object$cell_identity_allnames[which(seurat_object$cell_identity == "VZ_NPC")] <- "VZ_Neural_Stem_Cell"
seurat_object$cell_identity_allnames[which(seurat_object$cell_identity == "Cortex_Layer_VI_ExN")] <- "Cortex_Layer_VI_Excitatory_Neuron"
seurat_object$cell_identity_allnames[which(seurat_object$cell_identity == "Cortex_Layer_V_ExN")] <- "Cortex_Layer_V_Excitatory_Neuron"
seurat_object$cell_identity_allnames[which(seurat_object$cell_identity == "Cortex_Layer_IV_Spiny_ExN")] <- "Cortex_Layer_IV_Excitatory_Spiny_Neuron"
seurat_object$cell_identity_allnames[which(seurat_object$cell_identity == "Cortex_Layer_II_III_ExN")] <- "Cortex_Layer_II_III_Excitatory_Neuron"
seurat_object$cell_identity_allnames[which(seurat_object$cell_identity == "Cortex_CR")] <- "Cortex_Cajal-Retzius_cells"
seurat_object$cell_identity_allnames[which(seurat_object$cell_identity == "Cortex_SST_InN")] <- "Cortex_Somatostatin_Inhibitory_Neuron"
seurat_object$cell_identity_allnames[which(seurat_object$cell_identity == "Cortex_Pvalb_InN")] <- "Cortex_Pvalb_Inhibitory_Neuron"
seurat_object$cell_identity_allnames[which(seurat_object$cell_identity == "Cortex_PV_InN")] <- "Cortex_Parvalbumin_Inhibitory_Neuron"
seurat_object$cell_identity_allnames[which(seurat_object$cell_identity == "Cortex_Vip_InN")] <- "Cortex_Vasoactive_intestinal_peptide_Inhibitory_Neuron"
seurat_object$cell_identity_allnames[which(seurat_object$cell_identity == "Cortex_ODC")] <- "Cortex_Oligodendrocytes"
seurat_object$cell_identity_allnames[which(seurat_object$cell_identity == "CC_ODC")] <- "Corpus_Callosum_Oligodendrocytes"
seurat_object$cell_identity_allnames[which(seurat_object$cell_identity == "ODC")] <- "Oligodendrocytes"
seurat_object$cell_identity_allnames[which(seurat_object$cell_identity == "CPU_MSN")] <- "CPU_MSN"
seurat_object$cell_identity_allnames[which(seurat_object$cell_identity == "ChP_cilia")] <- "Choroid_plexus_cilia_cells"
seurat_object$cell_identity_allnames[which(seurat_object$cell_identity == "Ependymal")] <- "Ependymal_cells"
seurat_object$cell_identity_allnames[which(seurat_object$cell_identity == "Endo")] <- "Endothelial_cells"
seurat_object$cell_identity_allnames[which(seurat_object$cell_identity == "MG")] <- "microglia"
seurat_object$cell_identity_allnames[which(seurat_object$cell_identity == "Ery")] <- "Erythrocytes"

Idents(seurat_object) <- "cell_identity"
DefaultAssay(seurat_object) <- "SCT"

key_markers <- c("Tox3", "Pde1a", "Etl4", "Tle4", "Tafa1", "Ldb2", "Rorb", "Mef2c", "Cux2", "Rasgrf2", "Reln", "Dgkg", "Gad2", "Erbb4", "Pvalb", "Cdh13", "Grin3a", "Synpr", "Adarb2", "Mbp", "Mobp", "Tns3", "Plp1", "Penk", "Drd2", "Ttr", "Enpp2", "Tmem212", "Flt1", "Ptprc", "Ly86", "Hba-a1", "Hbb-bs", "Il31ra")

DotPlot(seurat_object, features = key_markers, cols = c("lightgrey", "red"), col.min = 0.1, dot.min = 0.1) + coord_flip() + RotatedAxis() + labs(y = "cell identity", x = "Features") + NoLegend()

qsave(seurat_object, "output/brain_Adult_ST_final_tutorial.qs")
#!----------------------------------------------------------------------------------------------------------------------------------------------------------
# StereoSeq data analysis
library(reticulate)
use_python("/home/yiyelinfeng/softwares/miniconda3/envs/Seurat5/bin")

sc <- import("scanpy")
sq <- import("squidpy")
np <- import("numpy")
pd <- import("pandas")
ad <- import("anndata")
setwd("StereoSeq/")

library(Seurat)
library(qs)
library(Matrix)
#!------------------------------------------------
adata <- sc$read_h5ad("Mouse_brain.h5ad")

count_matrix <- py_to_r(adata$layers[["count"]])
count_matrix <- as(count_matrix, "CsparseMatrix")
data_matrix <- py_to_r(adata$X)
data_matrix <- as(data_matrix, "CsparseMatrix")

rownames(count_matrix) <- rownames(adata$obs)
colnames(count_matrix) <- rownames(adata$var)

rownames(data_matrix) <- rownames(adata$obs)
colnames(data_matrix) <- rownames(adata$var)

annotation_colors <- py_to_r(adata$uns["annotation_colors"])[[1]]
spatial_umap <- adata$obsm[["spatial"]]
rownames(spatial_umap) <- rownames(adata$obs)
colnames(spatial_umap) <- c("spatial_1", "spatial_2")
spatial_umap[,2] <- -spatial_umap[,2]

seurat_brain <- CreateSeuratObject(counts = t(count_matrix), data = t(data_matrix), assay = "spatial", meta.data = py_to_r(adata$obs))
names(annotation_colors) <- unique(seurat_brain$annotation)
seurat_brain@misc$annotation_colors <- annotation_colors

seurat_brain[["spatial_umap"]] <- CreateDimReducObject(embeddings = spatial_umap, assay = "spatial")

seurat_brain$orig.ident <- "brain"
Idents(seurat_brain) <- "annotation"

DimPlot(seurat_brain, group.by = "annotation", cols = seurat_brain@misc$annotation_colors, pt.size = 0.5)
qsave(seurat_brain, "seurat_mouse_adult_brain_MOSTA.qs")
#!------------------------------------------------
adata <- sc$read_h5ad("Mouse_embryo_all_stage.h5ad")

count_matrix <- py_to_r(adata$layers[["count"]])
count_matrix <- as(count_matrix, "CsparseMatrix")
data_matrix <- py_to_r(adata$X)
data_matrix <- as(data_matrix, "CsparseMatrix")

rownames(count_matrix) <- rownames(adata$obs)
colnames(count_matrix) <- rownames(adata$var)

rownames(data_matrix) <- rownames(adata$obs)
colnames(data_matrix) <- rownames(adata$var)

annotation_colors <- py_to_r(adata$uns["annotation_colors"])[[1]]
spatial_umap <- adata$obsm[["spatial"]]
rownames(spatial_umap) <- rownames(adata$obs)
colnames(spatial_umap) <- c("spatial_1", "spatial_2")
spatial_umap[,2] <- -spatial_umap[,2]

seurat_object <- CreateSeuratObject(counts = t(count_matrix), data = t(data_matrix), assay = "spatial", meta.data = py_to_r(adata$obs))
names(annotation_colors) <- unique(seurat_object$annotation)
seurat_object@misc$annotation_colors <- annotation_colors

seurat_object[["spatial_umap"]] <- CreateDimReducObject(embeddings = spatial_umap, assay = "spatial")

seurat_object$orig.ident <- seurat_object$timepoint
Idents(seurat_object) <- "annotation"

DimPlot(seurat_object, group.by = "annotation", cols = seurat_object@misc$annotation_colors, pt.size = 0.5)
qsave(seurat_object, "seurat_mouse_embryo_MOSTA.qs")

library(UCell)

synapse_assembly <- read.csv("synapse_assembly_gene_list", header = F)[,1]

marker_genes <- list(synapse_assembly_full = synapse_assembly, synapse_assembly = key_genes)

seurat_E9 <- UCell::AddModuleScore_UCell(seurat_E9, features = marker_genes, ncores = 1, name = "")
seurat_E10 <- UCell::AddModuleScore_UCell(seurat_E10, features = marker_genes, ncores = 1, name = "")
seurat_E11 <- UCell::AddModuleScore_UCell(seurat_E11, features = marker_genes, ncores = 1, name = "")
seurat_E12 <- UCell::AddModuleScore_UCell(seurat_E12, features = marker_genes, ncores = 1, name = "")
seurat_E13 <- UCell::AddModuleScore_UCell(seurat_E13, features = marker_genes, ncores = 1, name = "")
seurat_E14 <- UCell::AddModuleScore_UCell(seurat_E14, features = marker_genes, ncores = 1, name = "")

synapse_assembly <- rbind(seurat_E9@meta.data[,c("synapse_assembly_full", "synapse_assembly")], seurat_E10@meta.data[,c("synapse_assembly_full", "synapse_assembly")], seurat_E11@meta.data[,c("synapse_assembly_full", "synapse_assembly")], seurat_E12@meta.data[,c("synapse_assembly_full", "synapse_assembly")], seurat_E13@meta.data[,c("synapse_assembly_full", "synapse_assembly")], seurat_E14@meta.data[,c("synapse_assembly_full", "synapse_assembly")])

FeaturePlot(seurat_object, reduction = "Spatial", features = "synapse_assembly_full", raster = FALSE) + SeuratAxes() + DarkTheme() + scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))[c(1:6, 9:11)])
ggsave("synapse_assembly_full.pdf")

FeaturePlot(seurat_object, reduction = "Spatial", features = "synapse_assembly", raster = FALSE) + SeuratAxes() + DarkTheme() + scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))[c(1:6, 9:11)])
ggsave("synapse_assembly.pdf")
#!----------------------------------------------------------------------------------------------------------------------------------------------------
# 10X HD analysis
seurat_object <- Load10X_Spatial(data.dir = "../../data/Visium_HD_Mouse_Brain_Fresh_Frozen/", bin.size = c(8), slice = "HD")

count.matrix <- Read10X_h5("../../data/Visium_HD_Mouse_Brain_Fresh_Frozen/binned_outputs/square_008um/filtered_feature_bc_matrix.h5", unique.features = F)
tmp_names <- as.data.frame(table(rownames(count.matrix)))
unique_genes <- as.character(tmp_names$Var1[which(tmp_names$Freq == 1)])

tmp_matrix <- count.matrix[-match(unique_genes, rownames(count.matrix)),]
count.matrix <- count.matrix[match(unique_genes, rownames(count.matrix)),]

ID <- as.character(rownames(tmp_matrix))
ID <- factor(ID, levels = unique(ID))
tmp_matrix <- as(rowsum(as.matrix(tmp_matrix), ID, reorder = FALSE, na.rm = TRUE), "dgCMatrix")
count.matrix <- rbind(count.matrix, tmp_matrix)

temp_object <- CreateSeuratObject(counts = count.matrix, project = "brain_HD", min.cells = 30, min.features = 100)
seurat_object <- subset(seurat_object, cells = colnames(temp_object))
temp_object@images$brain_HD <- seurat_object@images$HD.008um
temp_object@images$brain_HD@image <- png::readPNG("../../data/Visium_HD_Mouse_Brain_Fresh_Frozen/binned_outputs/square_008um/spatial/tissue_hires_image.png")
seurat_object <- temp_object
Seurat::SpatialPlot(seurat_object, features = "nCount_RNA", shape = 22, image.scale = "hires", crop = FALSE) & scale_fill_gradientn(colours = paletteContinuous(set = "sambaNight"))

seurat_object[["percent.mt"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^mt-")
seurat_object[["percent.ribo"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^Rp[sl]")

HB.genes <- rownames(seurat_object)[grep("^Hb[ab]-",rownames(seurat_object))]
seurat_object[["percent.hb"]] <- PercentageFeatureSet(object = seurat_object, features = HB.genes)
seurat_object <- subset(seurat_object, subset = percent.hb < 5 & nFeature_RNA < 3000 & nCount_RNA < 5000 & percent.mt < 15)
VlnPlot(seurat_object, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.hb"), pt.size = 0, ncol = 2) & geom_boxplot(color = "lightgrey", outlier.size = 0, width = 0.3) & NoLegend()
temp_object <- CreateSeuratObject(counts = GetAssayData(seurat_object, assay = "RNA", slot = "counts"), project = "brain_HD", min.cells = 1, min.features = 1)
#!--------------------------------------------------------------------------------------------
seurat_object <- qread("output/brain_Adult_HD_final_tutorial.qs")

seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)
seurat_object <- SketchData(object = seurat_object, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")

DefaultAssay(seurat_object) <- "sketch"
# perform clustering workflow
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object, assay = "sketch", reduction.name = "sketch_pca")
seurat_object <- FindNeighbors(seurat_object, assay = "sketch", reduction = "sketch_pca", annoy.metric = "cosine", dims = 1:50)
seurat_object <- FindClusters(seurat_object, cluster.name = "sketched_cluster", resolution = 3)
seurat_object <- RunUMAP(seurat_object, reduction = "sketch_pca", reduction.name = "sketch_umap", return.model = T, dims = 1:50)

seurat_object <- ProjectData(object = seurat_object, assay = "RNA", full.reduction = "full_pca", sketched.assay = "sketch", sketched.reduction = "sketch_pca", umap.model = "sketch_umap", dims = 1:50, refdata = list(seurat_cluster.projected = "sketched_cluster"))
DefaultAssay(seurat_object) <- "RNA"
Idents(seurat_object) <- "sketched_full_cluster"
DimPlot(seurat_object, reduction = "full_sketch_umap", label = T, raster = FALSE)
Seurat::SpatialPlot(seurat_object, image.scale = "hires", crop = FALSE)

library(SeuratWrappers)
library(Banksy)
seurat_object <- RunBanksy(seurat_object, lambda = 0.8, verbose = TRUE, assay = "RNA", slot = "data", features = "variable", k_geom = 24)
DefaultAssay(seurat_object) <- "BANKSY"
seurat_object <- RunPCA(seurat_object, assay = "BANKSY", reduction.name = "pca_banksy", npcs = 50, features = rownames(seurat_object))
seurat_object <- RunUMAP(seurat_object, reduction = "pca_banksy", dims = 1:50, min.dist = 0.3, reduction.name = "banksy_umap", return.model = T, n.neighbors = 24)
seurat_object <- FindNeighbors(seurat_object, reduction = "pca_banksy", dims = 1:50, annoy.metric = "cosine", graph.name = c('banksy_nn', 'banksy_snn'), k.param = 24)
seurat_object <- FindClusters(seurat_object, cluster.name = "banksy_cluster", graph.name = "banksy_snn", resolution = 0.5)

banksy_cells <- CellsByIdentities(seurat_object)
Seurat::SpatialPlot(seurat_object, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], shape = 22, image.scale = "hires", crop = FALSE, cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T, ncol = 8) + NoLegend()

library(Rmagic)
library(reticulate)
use_python("/home/yiyelinfeng/softwares/miniconda3/envs/R4/bin")
seurat_object <- spARC_Seurat(seurat_object, spatial_X = seurat_object@meta.data[,c("sdimy", "sdimx")])

cc.genes <- cc.genes.updated.2019
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
gene.mouse <- read.csv("/media/yiyelinfeng/data/Repository/baseFiles/ortholog_mouse.csv", header = T)
s.genes <- intersect(rownames(seurat_object), gene.mouse$Mouse.gene.name[match(s.genes, gene.mouse$Gene.name)])
g2m.genes <- intersect(rownames(seurat_object), gene.mouse$Mouse.gene.name[match(g2m.genes, gene.mouse$Gene.name)])
seurat_object <- CellCycleScoring(seurat_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
seurat_object$CC.Difference <- seurat_object$S.Score - seurat_object$G2M.Score

seurat_object <- RunPCA(seurat_object, npcs = 50, verbose = FALSE)
# one-liner to run Integration
#seurat_object <- IntegrateLayers(object = seurat_object, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = 'harmony', assay = "SCT", verbose = FALSE)
seurat_object <- FindNeighbors(seurat_object, reduction = "pca", dims = 1:50, annoy.metric = "cosine")
seurat_object <- FindClusters(seurat_object, resolution = seq(0.1, 0.8, 0.1))
seurat_object <- RunUMAP(seurat_object, reduction = "pca", dims = 1:50, reduction.name = "umap", return.model = T)

locs <- c()
for(i in names(seurat_object@images)){
	locs <- rbind(locs, Seurat::GetTissueCoordinates(seurat_object, image = i)[,seq_len(2)])
}
all(rownames(locs) == colnames(seurat_object))
colnames(locs) <- c("sdimy", "sdimx")
seurat_object@meta.data <- cbind(seurat_object@meta.data, locs)
#!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

