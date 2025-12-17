source /share/home/zhanglin/.bashrc

work_path=$1 # work_path=/media/maoni/data/CZP/spatial_hic/E13
EXP=$2 # E13_all
sample=$3  # E1305-hic-24
res=$4 # res=1Mb
chr=$5  # chr10

if [ -d $work_path ]; then

       ############################################### single pixel
       conda activate hicexplorer
       cd $work_path
       mkdir $work_path/$EXP/h5
       mkdir $work_path/$EXP/hic
       mkdir $work_path/$EXP/matrix
       mkdir $work_path/$EXP/RTI
       mkdir $work_path/$EXP/ABcom
       mkdir $work_path/$EXP/TADs
       mkdir $work_path/$EXP/Loops
       mkdir $work_path/$EXP/plots
       
### pre matrix
       hicTransform \
               --matrix $work_path/$EXP/h5/$sample\_$res\_Corrected.h5 \
               --outFileName $work_path/$EXP/ABcom/$sample\_$res\_Corrected_obs_exp.h5 \
               --method obs_exp
               
       hicTransform \
               --matrix $work_path/$EXP/h5/$sample\_$res\_Corrected_obs_exp.h5 \
               --outFileName $work_path/$EXP/ABcom/$sample\_$res\_Corrected_pearson_all.h5 \
               --method pearson
               
       hicPCA \
               --matrix $work_path/$EXP/h5/$sample\_$res\_Corrected_pearson_all.h5 \
               --outputFileName $work_path/$EXP/ABcom/$sample\_$res\_Corrected_pca1.bw $work_path/$EXP/ABcom/$sample\_$res\_Corrected_pca2.bw \
               --format bigwig
       hicPCA \
               --matrix $work_path/$EXP/h5/$sample\_$res\_Corrected_pearson_all.h5 \
               --outputFileName $work_path/$EXP/ABcom/$sample\_$res\_Corrected_pca1.bg $work_path/$EXP/ABcom/$sample\_$res\_Corrected_pca2.bg \
               --format bedgraph
               
### plot matrix and pc1 for corrected , pearson and O/E matrix
       hicPlotMatrix \
               --matrix $work_path/$EXP/h5/$sample\_$res\_Corrected.h5 \
               --outFileName $work_path/$EXP/plots/hicMatrixPlot/$sample\_$res\_Corrected_pca1.png \
               --perChr \
               --vMin 0 \
               --vMax 500 \
               --colorMap RdBu \
               --bigwig $work_path/$EXP/ABcom/$sample\_$res\_Corrected_pca1.bw
       
       hicPlotMatrix \
               --matrix $work_path/$EXP/ABcom/$sample\_$res\_Corrected_obs_exp.h5 \
               --outFileName $work_path/$EXP/plots/hicMatrixPlot/$sample\_$res\_Corrected_obs_exp_pca1.png \
               --perChr \
               --vMin 0 \
               --vMax 1 \
               --colorMap RdBu \
               --bigwig $work_path/$EXP/ABcom/$sample\_$res\_Corrected_pca1.bw
               
       hicPlotMatrix \
               --matrix $work_path/$EXP/ABcom/$sample\_$res\_Corrected_pearson_all.h5 \
               --outFileName $work_path/$EXP/plots/hicMatrixPlot/$sample\_$res\_Corrected_pearson_all_pca1.png \
               --perChr \
               --vMin -1 \
               --vMax 1 \
               --colorMap Reds \
               --bigwig $work_path/$EXP/ABcom/$sample\_$res\_Corrected_pca1.bw
               
               
### region plot matrix and pc1 for corrected , pearson and O/E matrix                      
       hicPlotMatrix \
               --matrix $work_path/$EXP/h5/$sample\_$res\_Corrected.h5 \
               --outFileName $work_path/$EXP/plots/hicMatrixPlot/$sample\_$res\_Corrected_pca1_$chr.png \
               --region $chr \
               --vMin 0 \
               --vMax 500 \
               --colorMap RdBu \
               --bigwig $work_path/$EXP/ABcom/$sample\_$res\_Corrected_pca1.bw
                            
       hicPlotMatrix \
               --matrix $work_path/$EXP/ABcom/$sample\_$res\_Corrected_obs_exp.h5 \
               --outFileName $work_path/$EXP/plots/hicMatrixPlot/$sample\_$res\_Corrected_obs_exp_pca1_$chr.png \
               --region $chr \
               --vMin 0 \
               --vMax 1 \
               --colorMap RdBu \
               --bigwig $work_path/$EXP/ABcom/$sample\_$res\_Corrected_pca1.bw
                              
       hicPlotMatrix \
               --matrix $work_path/$EXP/ABcom/$sample\_$res\_Corrected_pearson_all.h5 \
               --outFileName $work_path/$EXP/plots/hicMatrixPlot/$sample\_$res\_Corrected_pearson_all_pca1_$chr.png \
               --region $chr \
               --vMin -1 \
               --vMax 1 \
               --colorMap Reds \
               --bigwig $work_path/$EXP/ABcom/$sample\_$res\_Corrected_pca1.bw  
else
	echo "Please give the mother_path(of the project name) and the project name in order!
	Just like:bash xxx.sh [/media/xxx] [Spatial-E13]
fi

