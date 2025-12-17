library(getopt)
command=matrix(c("sample", "s", 1, "character",
                 "QC_table_file", "q", 1, "character",
                 "check_linker_file", "l", 1, "character",
                 "stat_clean_file","f",1,"character",
                 "work_path","p",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)

args=getopt(command)
if (!is.null(args$help) || is.null(args$sample) || is.null(args$QC_table_file) || is.null(args$check_linker_file) || is.null(args$stat_clean_file) || is.null(args$work_path)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

sample = args$sample  # "E13-8-HiC"
work_path = args$work_path # "/date/guomaoni/CZP/spatial_hic_rna/HiC/E13-8-HiC" 
QC_table_file = args$QC_table_file # "/date/guomaoni/CZP/spatial_hic_rna/HiC/E13-8-HiC/3_pairtools/E13-8-HiC.dedup.stats"
check_linker_file = args$check_linker_file # "/date/guomaoni/CZP/spatial_hic_rna/HiC/E13-8-HiC/0_check_linker/check_linker.log"
stat_clean_file = args$stat_clean_file # "/date/guomaoni/CZP/spatial_hic_rna/HiC/E13-8-HiC/3_pairtools/E13-8-HiC.dedup.stats"

# sample = "E185-brain-hic-15"
# work_path = "/date/guomaoni/CZP/spatial_hic/E185-brain-hic-15"
# QC_table_file = "/date/guomaoni/CZP/spatial_hic/E185-brain-hic-15/3_hicMatrix/E185-brain-hic-15_10kb_QC/QC_table.txt"
# check_linker_file = "/date/guomaoni/CZP/spatial_hic/E185-brain-hic-15/0_check_linker/check_linker.log"
# stat_clean_file = "/date/guomaoni/CZP/spatial_hic/E185-brain-hic-15/5_bedpe/stat_clean.csv"

# linker
setwd(work_path)
spot_size = 50
stat_linker = read.table(check_linker_file, sep = "\t", header = FALSE)
stat_QC_table = t(read.table(QC_table_file, sep = "\t", header = FALSE))
stat_all_1 = as.matrix(c(sample,
                         spot_size,
                         as.character(stat_linker[3,]),
                         paste0("=", stat_linker[11,], "/", stat_linker[3,]),
                         paste0("=", stat_linker[13,], "/", stat_linker[11,]),
                         paste0("=", stat_QC_table[2,2], "/", stat_linker[11,]),
                         paste0("=", stat_QC_table[5,2], "/", stat_QC_table[2,2]),
                         paste0("=", stat_QC_table[15,2], "/", stat_QC_table[5,2]),
                         paste0("=", stat_QC_table[13,2], "/", stat_QC_table[5,2]),
                         paste0("=", stat_QC_table[6,2], "/", stat_QC_table[5,2]),
                         stat_QC_table[6,2],
                         paste0("=", stat_QC_table[16,2], "/", stat_QC_table[6,2]), 
                         paste0("=", stat_QC_table[17,2], "/", stat_QC_table[6,2]),
                         paste0("=", stat_QC_table[18,2], "/", stat_QC_table[6,2]),
                         stat_QC_table[18,2],
                         paste0("=", stat_QC_table[18,2], "/", stat_linker[3,])))

colnames(stat_all_1)=sample
rownames(stat_all_1)=c("sample", "um", "raw_reads", "ligation", "barcode", "bwa", "pairs", "duplication", "same fragment", "valid reads", "valid reads", "trans", "cis short <20kb", "cis long >=20kb", "cis long", "cis long/raw")

spot_NO = 50*50
stat = read.csv(stat_clean_file)
stat_all_2 = as.matrix(round(c(spot_NO, median(stat$raw_count), mean(stat$raw_count),median(stat$valid_contact), mean(stat$valid_contact),median(stat$cis_contact), mean(stat$cis_contact),median(stat$cis_more_10kb), mean(stat$cis_more_10kb),
                               sum(stat$raw_count), sum(stat$valid_contact), sum(stat$trans_contact), sum(stat$cis_less_10kb), sum(stat$cis_more_10kb))))
colnames(stat_all_2)=sample
rownames(stat_all_2)=c("pixel num", "median_raw_count", "mean_raw_count", "median_valid", "mean_valid", "median_cis", "mean_cis", "median_cis_long", "mean_cis_long", 
                       "total_raw_count", "total_valid", "total_trans", "total_cis_short", "total_cis_long")

stat_all = rbind(stat_all_1, stat_all_2)
stat_all = data.frame(rownames(stat_all), stat_all[,1])
colnames(stat_all) <- c("stat", sample)
write.table(stat_all, paste0(work_path, "/", sample, "_stat.txt"), sep="\t", col.names=F, row.names=F, quote=F)



