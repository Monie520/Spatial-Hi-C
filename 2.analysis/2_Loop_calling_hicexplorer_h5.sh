. /date/scripts/Pipelines/temp_bashrc

work_path=$1 # /date/guomaoni/CZP/spatial_hic/hicexplorer/Cere
EXP=$2  # Cere_all
sample=$3  # sample=Cere-hic-13 E13-hic-com
res=$4 # res=10kb


if [ -d $work_path ]; then

       ############################################### single pixel
       conda activate hicexplorer
       cd $work_path
       mkdir $work_path/$EXP/h5
       mkdir $work_path/$EXP/cool
       mkdir $work_path/$EXP/hic
       mkdir $work_path/$EXP/matrix
       mkdir $work_path/$EXP/RTI
       mkdir $work_path/$EXP/ABcom
       mkdir $work_path/$EXP/TADs
       mkdir $work_path/$EXP/Loops
       mkdir $work_path/$EXP/plots

       cd $work_path/$EXP/h5
       hicCorrectMatrix correct \
               --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
               --matrix $work_path/$EXP/h5/$sample\_$res.h5 \
               --correctionMethod KR \
              --perchr \
              --outFileName $work_path/$EXP/h5/$sample\_$res\_Corrected.h5

       hicConvertFormat -m $work_path/$EXP/h5/$sample\_$res\_Corrected.h5 \
               --inputFormat h5 \
               --outputFormat cool \
               -o $work_path/$EXP/cool/$sample\_$res\_Corrected.cool
       
       # call loop
       peakWidth=("2" "2" "3" "3" "4" "4" "6" "6" "6" "6" "6")
       windowSize=("5" "9" "10" "7" "10" "15" "10" "10" "10" "10" "15")
       maxLoopDistance=("2000000" "2000000" "2000000" "5000000" "2000000" "2000000" "2000000" "2000000" "2000000" "5000000" "2000000")
       obsExpThreshold=("1.5" "1.5" "1.5" "1.5" "1.5" "1.5" "1" "1" "1" "1" "1")
       pValuePreselection=("0.1" "0.1" "0.1" "0.1" "0.1" "0.1" "0.05" "0.05" "0.05" "0.05" "0.1")
       expected=("mean" "mean_nonzero" "mean" "mean" "mean" "mean" "mean_nonzero" "mean_nonzero_ligation" "mean_nonzero_ligation" "mean_nonzero" "mean")
       peakInteractionsThreshold=("10" "3" "3" "3" "3" "3" "5" "3" "5" "3" "3")
       
       for ((i = 0; i < ${#peakWidth[@]}; i++))
       do
       hicDetectLoops -m $work_path/$EXP/cool/$sample\_$res\_Corrected.cool \
               -o $work_path/$EXP/Loops/$sample\_$res\_Corrected.p${peakWidth[i]}win${windowSize[i]}.maxD_${maxLoopDistance[i]}.minOE${obsExpThreshold[i]}.p${pValuePreselection[i]}.expected_${expected[i]}.minF${peakInteractionsThreshold[i]}.v36.loops.bedgraph \
               --peakWidth ${peakWidth[i]} \
               --windowSize ${windowSize[i]} \
               --maxLoopDistance ${maxLoopDistance[i]} \
               --obsExpThreshold ${obsExpThreshold[i]} \
               --pValuePreselection ${pValuePreselection[i]} \
               --expected ${expected[i]} \
               --peakInteractionsThreshold ${peakInteractionsThreshold[i]} \
               --threads 6 \
               --threadsPerChromosome 4
       done
       
       # merge loop
       hicMergeLoops -i $work_path/$EXP/Loops/* -o $sample\_merged_loops.bedgraph -r 10000
       
       # loop strength
        
else
        echo "Please give the mother_path(of the project name) and the project name in order!
        Just like:bash xxx.sh [/media/xxx] [Spatial-E13]
fi



