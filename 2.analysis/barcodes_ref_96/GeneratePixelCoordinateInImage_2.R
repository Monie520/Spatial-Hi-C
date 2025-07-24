# Order 2 (real)
# 1x1 1x2 ... 1x96
# 2x1 2x2 ... 2x96
name = "big_chip"
hei = 1080
wid = 1080
path = "/media/maoni/data/CZP/spatial_hic/seurat_yuhao/core/96*96_blank_img"

comb.9696 = data.frame("i1"=rep(1:96,each=96),
                       "i2"=rep(1:96,times=96),stringsAsFactors = F)

comb.9696$comb_i = paste0(comb.9696$i1,"x",comb.9696$i2)
comb.9696$tissue = name
comb.9696$imagerow = (97 - comb.9696$i1 - 0.5)*hei/96
comb.9696$imagecol = (comb.9696$i2 - 0.5)*wid/96

write.table(comb.9696[,c("comb_i","tissue","i1","i2","imagerow","imagecol")],
            file = paste0(path,"/combine_barcode.round2round1_index1_index2.Seurat_2.txt"),
            row.names = F,col.names = F,sep = "\t",quote = F)
