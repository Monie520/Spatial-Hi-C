
work_path=$1
EXP=$2
sample=$3 # Spatial-E13-20-20220127 E1305-hic-24 E13-hic-com


if [ -d $work_path ]; then

       ############################################### single pixel
       conda activate hicexplorer3.7
       
       mkdir $work_path/$EXP/TADs
      
       cd $work_path/$EXP/h5
       
       # for sample in E13-hic-com_C1 # E13-hic-13_C7 E13-hic-13_C8 E13-hic-com_C7 E13-hic-com_C8
       # do
       hicFindTADs -m $work_path/$EXP/h5/${sample}_40kb_Corrected.h5 \
               --outPrefix $work_path/$EXP/TADs/${sample}_40kb_Corrected_TAD \
               --minDepth 300000 \
               --maxDepth 3000000 \
               --step 300000 \
               --minBoundaryDistance 400000 \
               --thresholdComparisons 0.01 \
               --correctForMultipleTesting fdr \
               --delta 0.01      
       # done
else
	echo "Please give the mother_path(of the project name) and the project name in order!
	Just like:bash xxx.sh [/media/xxx] [Spatial-E13]"
fi




