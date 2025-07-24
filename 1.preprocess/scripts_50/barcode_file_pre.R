library(getopt)
command=matrix(c("work_path","p",1,"character",
                 "file_name_1","a",1,"character",
                 "file_name_2","b",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)

args=getopt(command)
if (!is.null(args$help) || is.null(args$work_path) || is.null(args$file_name_1) || is.null(args$file_name_2)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

work_path = args$work_path
barcodes_A = args$file_name_1
barcodes_B = args$file_name_2
# 
# work_path = "/media/maoni/data/CZP/spatial_hic/E14-brain-hic-15/barcode_ref"
# barcodes_A = "/media/maoni/data/CZP/spatial_hic/E14-brain-hic-15/barcode_ref/barcodes_A.txt"
# barcodes_B = "/media/maoni/data/CZP/spatial_hic/E14-brain-hic-15/barcode_ref/barcodes_B.txt"

setwd(work_path)
# 0_barcode_total.HiC.V1.txt (barcode file)
# 1_barcode_add_order.n1_n2.v1.txt
file0=read.table(barcodes_B,sep="\t",header=F,stringsAsFactors = F)
file1=read.table(barcodes_A,sep="\t",header=F,stringsAsFactors = F)
n=seq(1,dim(file0)[1],1)
file2=data.frame(n,n)
write.table(file2, "1_barcode_add_order.n1_n2.txt",sep="\t",col.names=F,row.names=F,quote=F)

# 2_combine_barcode.round2round1_index1_index2.v1.txt
file3=c()
for(i in 1:(dim(file0)[1])){
  col1=paste(file0[i,1], file1[,1],sep="")
  col2=i
  col3=file2[,1]
  col123=data.frame(col1,col2,col3)
  file3=rbind(file3,col123)
}
write.table(file3, "2_combine_barcode.round2round1_index1_index2.txt",sep="\t",col.names=F,row.names=F,quote=F)

