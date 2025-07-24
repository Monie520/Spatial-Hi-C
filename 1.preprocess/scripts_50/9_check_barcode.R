library(getopt)
command=matrix(c("combine_barcode","b",1,"character",
                 "work_path","p",1,"character",
                 "file_name","n",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)

args=getopt(command)
if (!is.null(args$help) || is.null(args$combine_barcode) || is.null(args$work_path) || is.null(args$file_name)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

work_path = args$work_path # "/asnas/liujiang_group/yuhao/SingleCell_Space/27_Adult_unstain_20um_RNA_1110/scRNA/fastq"
inpu_name = args$file_name # "RNA_Adult_1110_2.extract.barcode"
refe_path = args$combine_barcode # "/asnas/liujiang_group/yuhao/SingleCell_Space/27_Adult_unstain_20um_RNA_1110/out/combine_barcode.round2round1_index1_index2.reorder.txt"

# 
# work_path = "/media/maoni/data/CZP/spatial_hic/E14-brain-hic-15/0_check_linker"
# inpu_name = "/media/maoni/data/CZP/spatial_hic/E14-brain-hic-15/2_debarcode/barcodeB_A.read1.txt"
# refe_path =  "/media/maoni/data/CZP/spatial_hic/E14-brain-hic-15/barcode_ref/2_combine_barcode.round2round1_index1_index2.txt"


setwd(work_path)
library(ggplot2)
library(ComplexHeatmap)
# 1
ref_barcode = read.table(refe_path,stringsAsFactors = F)
our_barcode = read.table(inpu_name,stringsAsFactors = F)$V1
print("total barcode num")
length(our_barcode) # 66653845
total_barcode_num <- cbind("total barcode num", length(our_barcode))
write.table(total_barcode_num, paste(work_path,"/check_linker.log",sep=""), sep="\t",col.names=F,row.names=F, quote=F, append = T)

# 2
our_barcode=gsub('[_]', '', our_barcode)
our_barcode = our_barcode[our_barcode %in% ref_barcode$V1]
print("our barcode num")
length(our_barcode) # 65125203
barcode_in_ref_list <- cbind("our barcode num", length(our_barcode))
write.table(barcode_in_ref_list, paste(work_path,"/check_linker.log",sep=""), sep="\t",col.names=F,row.names=F, quote=F, append = T)

our_barcode = factor(our_barcode,levels = ref_barcode$V1)
ref_barcode$count = table(our_barcode)

nchannel = max(ref_barcode$V2)
all_combi.m = matrix(0,nchannel,nchannel) # barcodeA*barcodeB
row.names(all_combi.m) = paste0("A",1:nchannel)
colnames(all_combi.m) = paste0("B",1:nchannel)
tmp <- lapply(1:dim(ref_barcode)[1], function(i) all_combi.m[ref_barcode$V3[i],ref_barcode$V2[i]] <<- ref_barcode$count[i])


pdf("0_barcode_combination.pdf",width = 15,height = 15)
Heatmap(all_combi.m,cluster_rows=F,cluster_columns=F)
Heatmap(log2(all_combi.m+1),cluster_rows=F,cluster_columns=F)
dev.off()



library(circlize)
svg("1_barcode_combination.svg",width = 6,height = 6)
Heatmap(all_combi.m,cluster_rows=F,cluster_columns=F,rect_gp = gpar(col = "white", lwd = 5),
        col = colorRamp2(c(0, max(all_combi.m)/1000, max(all_combi.m)), c("blue", "green", "red")),
        show_heatmap_legend = F,show_row_names=F,show_column_names=F)
dev.off()


barcodeA_percent =rowSums(all_combi.m)/sum(all_combi.m)*100
barcodeB_percent =colSums(all_combi.m)/sum(all_combi.m)*100
result.df=data.frame(barcodeA_percent,barcodeB_percent)
write.csv(result.df, "4_barcode_proportion.csv")

# https://zhuanlan.zhihu.com/p/96444730


