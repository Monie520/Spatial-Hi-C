# Order 2 (real)
# 1x50 1x49 ... 1x1
# 2x50 2x49 ... 2x1
name = "small_chip"
hei = 1080
wid = 1080
path = "/media/maoni/data/CZP/spatial_hic/seurat_yuhao/core/50*50_blank_img"

comb.5050 = data.frame("i1"=rep(1:50,each=50),
                       "i2"=rep(1:50,times=50),stringsAsFactors = F)
comb.5050$comb_i = paste0(comb.5050$i1,"x",comb.5050$i2)
comb.5050$tissue = name
comb.5050$imagerow = (comb.5050$i1 - 0.5)*hei/50
comb.5050$imagecol = (51 - comb.5050$i2 - 0.5)*wid/50

write.table(comb.5050[,c("comb_i","tissue","i1","i2","imagerow","imagecol")],
            file = paste0(path,"/combine_barcode.round2round1_index1_index2.Seurat.txt"),
            row.names = F,col.names = F,sep = "\t",quote = F)

