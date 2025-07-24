
conda activate hicexplorer3.7

work_path=$1 # work_path=/media/maoni/data/CZP/spatial_hic/hicexplorer/E13
EXP=$2 # E13_all
method=$3  # pearson,spearman
figformat=$4  # png,pdf
range=500000:7000000 # set according to the size of interacting domains(Chenlab)

sample1=$5
sample2=$6

############################################ 200k raw corr  ######################################
            
### 200k (all samples)
cd $work_path/$EXP/h5/spatial_norm
hicCorrelate \
            -m `ls *_200kb.h5` \
            --method=$method \
            --log1p \
            --labels `ls *_200kb.h5` \
            --range $range \
            --zMin 0 \
            --zMax 1.0 \
            --plotNumbers \
            --outFileNameHeatmap $work_path/$EXP/plots/hicCorrelatePlot/Spatial_E13_all_heatmap_200k_$method.$figformat --colorMap YlGnBu \
            --outFileNameScatter $work_path/$EXP/plots/hicCorrelatePlot/Spatial_E13_all_scatterplot_200k_$method.$figformat
            # --plotFileFormat $figformat

### 200k (two samples)        
cd $work_path/$EXP/h5/spatial_norm
hicCorrelate \
            -m $sample1\_200kb.h5 $sample2\_200kb.h5 \
            --method=$method \
            --log1p \
            --labels $sample1 $sample2 \
            --range $range \
            --zMin 0 \
            --zMax 1.0 \
            --plotNumbers \
            --outFileNameHeatmap $work_path/$EXP/plots/hicCorrelatePlot/$sample1\_$sample2\_corr_heatmap_200k_$method.$figformat --colorMap YlGnBu \
            --outFileNameScatter $work_path/$EXP/plots/hicCorrelatePlot/$sample1\_$sample2\_corr_scatterplot_200k_$method.$figformat
            # --plotFileFormat $figformat
 
 
 
 
 
 
 
conda activate hicexplorer3.7

work_path=/media/maoni/data/CZP/spatial_hic/hicexplorer/E13
EXP=E13_all
method=pearson
figformat=pdf
range=500000:7000000 # set according to the size of interacting domains(Chenlab)


sample1=E1305-hic-24
sample2=E13-hic-13
sample3=E13-18-HiC
   
     
############################################ 200k normalize corr (Final)  ######################################            
            
### 200k (all samples)
cd $work_path/$EXP/h5/spatial_norm
hicCorrelate \
            -m `ls *_200kb_normalized.h5` \
            --method=$method \
            --log1p \
            --labels `ls *_200kb_normalized.h5` \
            --range $range \
            --zMin 0 \
            --zMax 1.0 \
            --plotNumbers \
            --outFileNameHeatmap $work_path/$EXP/plots/hicCorrelatePlot/norm2/Spatial_E13_all_normalized_heatmap_200k_$method.$figformat --colorMap YlGnBu \
            --outFileNameScatter $work_path/$EXP/plots/hicCorrelatePlot/norm2/Spatial_E13_all_normalized_scatterplot_200k_$method.$figformat
            # --plotFileFormat $figformat
            
### 200k (two samples)        
cd $work_path/$EXP/h5/spatial_norm
hicCorrelate \
            -m $sample2\_200kb_normalized.h5 $sample3\_200kb_normalized.h5 \
            --method=$method \
            --log1p \
            --labels $sample2 $sample3 \
            --range $range \
            --zMin 0 \
            --zMax 1.0 \
            --plotNumbers \
            --outFileNameHeatmap $work_path/$EXP/plots/hicCorrelatePlot/norm2/$sample2\_$sample3\_corr_normalized_heatmap_200k_$method.$figformat --colorMap YlGnBu \
            --outFileNameScatter $work_path/$EXP/plots/hicCorrelatePlot/norm2/$sample2\_$sample3\_corr_normalized_scatterplot_200k_$method.$figformat
            # --plotFileFormat $figformat    
            
            
            
