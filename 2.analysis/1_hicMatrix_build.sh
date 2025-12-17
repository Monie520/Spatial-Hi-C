work_path=$1 # work_path=/media/maoni/data/CZP/spatial_hic/hicexplorer/E13
EXP=$2 # E13_all
sample=$3  # E1305-hic-24
res=$4 # res=40kb
resolution=$5 # resolution=40000
fold=`expr $resolution / 10000`

juicer_path=/media/maoni/data/Software/juicer/scripts/common
Reference=/media/maoni/data/CZP/HiC-CXP-pipeline/Reference
species=mouse
genome=mm10
enzyme=Mbol
q=30

if [ -d $work_path ]; then

       ############################################### single pixel
       conda activate hicexplorer3.7
       cd $work_path
       mkdir $work_path/$EXP/h5
       mkdir $work_path/$EXP/hic
       mkdir $work_path/$EXP/matrix
       mkdir $work_path/$EXP/RTI
       mkdir $work_path/$EXP/ABcom
       mkdir $work_path/$EXP/TADs
       mkdir $work_path/$EXP/Loops
       mkdir $work_path/$EXP/plots

       cd $work_path/$EXP/h5
             
       hicMergeMatrixBins \
               --matrix $work_path/$EXP/h5/$sample\_10kb.h5 --numBins $fold \
               --outFileName $work_path/$EXP/h5/$sample\_$res.h5
               
       hicCorrectMatrix correct \
               --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
               --matrix $work_path/$EXP/h5/$sample\_$res.h5 \
               --correctionMethod KR \
               --perchr \
               --outFileName $work_path/$EXP/h5/$sample\_$res\_Corrected.h5
      
       hicConvertFormat -m $work_path/$EXP/h5/$sample\_$res\_Corrected.h5 \
               --inputFormat h5 \
               --outputFormat ginteractions \
               -o $work_path/$EXP/hic/$sample\_$res\_Corrected.ginteractions
       
       hicConvertFormat -m $work_path/$EXP/h5/$sample\_$res\_Corrected.h5 \
               --inputFormat h5 \
               --outputFormat cool \
               -o $work_path/$EXP/cool/$sample\_$res\_Corrected.cool
                            
       awk -F "\t" '{print 0, $1, $2, 0, 0, $4, $5, 1, $7}' $work_path/$EXP/hic/$sample\_$res\_Corrected.ginteractions.tsv > $work_path/$EXP/hic/$sample\_$res\_Corrected.ginteractions.tsv.short
       
       sort -k2,2d -k6,6d $work_path/$EXP/hic/$sample\_$res\_Corrected.ginteractions.tsv.short > $work_path/$EXP/hic/$sample\_$res\_Corrected.ginteractions.tsv.short.sorted
       
       java -jar $juicer_path/juicer_tools.jar pre -r $resolution -q $q -f $Reference/$species/$genome/juicer/${genome}_${enzyme}.txt $work_path/$EXP/hic/$sample\_$res\_Corrected.ginteractions.tsv.short.sorted $work_path/$EXP/hic/$sample\_$res\_Corrected.ginteractions.tsv.short.sorted.q$q.hic $genome

       # transfer to matrix
       chr=chr6
       start=94000000
       end=100000000
       
       cd $work_path/$EXP/matrix
       h1d basic dump $work_path/$EXP/hic/$sample\_$res\_Corrected.ginteractions.tsv.short.sorted.q$q.hic $resolution all \
               --gt $Reference/$species/$genome/mm10_genome_table.txt \
               --normalize NONE \
               -o $sample\_$res \
               --datatype rawhic \
               --maxchr 19 \
               -n 8
               
       $work_path/$EXP/plots/heatmap      
       h1d basic plot $work_path/$EXP/matrix/$sample\_$res/$resolution/observed.NONE.$chr.matrix.gz \
               $resolution $chr --datatype matrix -o $work_path/$EXP/plots/heatmap/$sample\_$res\_$chr\_$start\_$end --plottype square \
               -s $start -e $end
	
else
	echo "Please give the mother_path(of the project name) and the project name in order!
	Just like:bash xxx.sh [/media/xxx] [Spatial-E13]
	Or contact with kim@wechat:13936277609"
fi

