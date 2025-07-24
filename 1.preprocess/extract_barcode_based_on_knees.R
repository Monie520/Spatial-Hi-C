### BiocManager::install("GNET2")
library(GNET2)
library(ggplot2)
library(getopt)

command=matrix(c("sample", "s", 1, "character",
                 "stat_cropped_file","f1",1,"character",
                 "stat_cropped_filtered_file","f2",1,"character",
                 "knee_point_fig","k",1,"character",
                 "work_path","p",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)

args=getopt(command)
if (!is.null(args$help) || is.null(args$sample) || is.null(args$stat_cropped_file) || is.null(args$stat_cropped_filtered_file) || is.null(args$knee_point_fig) || is.null(args$work_path)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

sample = args$sample  # "E185-brain-hic-4"
work_path = args$work_path # "/date/guomaoni/CZP/spatial_hic/E185-brain-hic-4" 
stat_cropped_file = args$stat_cropped_file # "/date/guomaoni/CZP/spatial_hic/E185-brain-hic-4/6_crop/stat_cropped.csv"
stat_cropped_filtered_file = args$stat_cropped_filtered_file # "/media/maoni/data/CZP/spatial_hic/E185-brain-hic-4/6_crop/stat_cropped_filtered.csv"
knee_point_fig = args$knee_point_fig # "/media/maoni/data/CZP/spatial_hic/E185-brain-hic-4/6_crop/stat_cropped_filtered.csv"


### filter by raw valid_contact
setwd(work_path)
stat_cropped <- read.csv(stat_cropped_file)
stat_cropped <- stat_cropped[stat_cropped$valid_contact > 0, ]
stat_cropped$valid_contact_log <- log10(stat_cropped$valid_contact)

stat_cropped_sort <- stat_cropped[order(stat_cropped$valid_contact, decreasing = TRUE), ]
plot(stat_cropped_sort$valid_contact, type="l")
knee_point_up <- kneepointDetection(stat_cropped_sort$valid_contact)
knee_point_up

stat_cropped_sort <- stat_cropped[order(stat_cropped$valid_contact_log, decreasing = TRUE), ]
plot(stat_cropped_sort$valid_contact_log, type="l")
knee_point_down <- kneepointDetection(stat_cropped_sort$valid_contact_log)
knee_point_down

stat_cropped_sort_2 <- stat_cropped_sort[knee_point_up:(dim(stat_cropped_sort)[1]), ]
dim(stat_cropped_sort_2)

pdf(knee_point_fig, width = 7, height = 7)
plot(stat_cropped_sort_2$valid_contact,
     type="l", 
     lwd=2,
     ylab="valid_contact", 
     main=paste0(sample, "\n", "(total=", length(stat_cropped_sort_2$valid_contact), ", filtered=", knee_point_down-knee_point_up+1, ")"), 
     log="y")  # 设置y轴为对数刻度
abline(v = knee_point_down-knee_point_up, col = "red", lwd = 2, lty = 2)
dev.off()

pdf(knee_point_fig_2, width = 7, height = 7)
plot(stat_cropped_sort$valid_contact,
     type="l", 
     lwd=2, 
     ylab="valid_contact", 
     main=paste0(sample, "\n", "(total=", length(stat_cropped_sort$valid_contact), ", filtered=", knee_point_down-knee_point_up+1, ")"), 
     log="y")  # 设置y轴为对数刻度
abline(v = knee_point_up, col = "red", lwd = 2, lty = 2)
abline(v = knee_point_down, col = "red", lwd = 2, lty = 2)
dev.off()

stat_cropped_filter <- stat_cropped_sort[knee_point_up:knee_point_down, ]
write.csv(stat_cropped_filter, row.names = F, stat_cropped_filtered_file)

stat_cropped_filter_higashi <- stat_cropped[which(stat_cropped$pixel %in% stat_cropped_filter$pixel), 1:12]
write.csv(stat_cropped_filter_higashi, row.names = F, stat_cropped_filtered_file_higashi)

