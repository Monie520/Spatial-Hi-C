#!/usr/bin/bash

###########################################################################################################
####--------------------------------------Script Input Parameters--------------------------------------####
###########################################################################################################

# Check input script
if [ -e $1 ];
then
        source $1
else
        echo "Please check the existence of script files!"
        exit
fi


. /date/scripts/Pipelines/temp_bashrc
###########################################################################################################
####-------------------------------------Confirm Reference Species-------------------------------------####
###########################################################################################################

if [ $species == "human" ];
then
	ref_gtf=$gtf_human
	STAR_idx=$STARindex_human
elif [ $species == "mouse" ];
then
	ref_gtf=$gtf_mouse
	STAR_idx=$STARindex_mouse
elif [ $species == "mix" ];
then
	ref_gtf=$gtf_mix
	STAR_idx=$STARindex_mix
else
	echo "Your reference species does not exist!"
	exit
fi

###########################################################################################################
####--------------------------------------Check Directories/Files--------------------------------------####
###########################################################################################################

# Check output directory
if [ ! -d $outdir ];
then
	mkdir -p $outdir # -p, ensure creating sub-directories automatically
fi

cd $outdir

# Recreate log file
if [ -f $outdir/${sample}.checked.txt ];
then
        rm -rf $outdir/${sample}.checked.txt
        echo "" > $outdir/${sample}.checked.txt
fi

# Welcome
if [ -e $welcome ];
then
        cat $welcome
        cat "${welcome}" >> $outdir/${sample}.checked.txt
        echo "** pipeline start ** @ "`date`  >> $outdir/${sample}.checked.txt
fi

# Check input fq files
if [ -e $fq1 ] && [ -e $fq2 ];
then
        ln -sf $fq1 ${sample}_R1.fastq.gz
        ln -sf $fq2 ${sample}_R2.fastq.gz
        fq1=${outdir}/${sample}_R1.fastq.gz
        fq2=${outdir}/${sample}_R2.fastq.gz
else
        echo "Please check the existence of fastq files!"
        exit
fi

# Check barcode files
if [ -e $barcodes_A ] && [ -e $barcodes_B ];
then
        ln -sf $barcodes_A barcodesA.txt
        ln -sf $barcodes_B barcodesB.txt
        barcodes_A=${outdir}/barcodesA.txt
        barcodes_B=${outdir}/barcodesB.txt
        
        sed -i 's/\r//g' $barcodes_A # s/ A/B /g, global, each line
        sed -i 's/\r//g' $barcodes_B
        
        nA=$(cat $barcodes_A | wc -l)
        nB=$(cat $barcodes_B | wc -l)
        
        if [ $nA -eq $nB ];
	then
		if [ $npixel -eq $nA ] && [ $npixel -eq $nA ];
		then
			barcode_files=${outdir}/barcodes.txt
			if [ ! -f $barcode_files ];
			then
				cat $barcodes_B | while read bB;
				do
					cat $barcodes_A | while read bA;
					do
						echo $bB$bA >> $barcode_files
					done
				done
			fi
		else
			echo "You set npixel is $npixel, but your barcodes_A is $nA, barcodes_B is $nB!"
			exit
		fi
	else
		echo "Your barcodes files are note the same length, your barcodes_A is $nA, barcodes_B is $nB!"
		exit
	fi
else
        echo "Please check the existence of barcodes files!"
        exit
fi

###########################################################################################################
####--------------------------------------Start Running Pipeline---------------------------------------####
###########################################################################################################

if [ $run_fastqc_every_step == "run" ];
then
	if [ ! -d $outdir/fastqc_res ];
	then
		mkdir -p $outdir/fastqc_res
        fi
        
        cd $outdir/fastqc_res
        
        echo "** fastqc raw fq start ** @ "`date`  >> $outdir/${sample}.checked.txt
        
        time $fastqc --java $java_jre -t $ppn -o ./ $fq1 $fq2 --noextract 1>>${outdir}/${sample}_log.mark 2>&1
        
        echo "** fastqc raw fq done ** @ "`date`  >> $outdir/${sample}.checked.txt
fi

###########################################################################################################
####-------------------------------------------Check Linkers-------------------------------------------####
###########################################################################################################

if [ $run_check_linker == "run" ];
then
	if [ ! -d $outdir/1_check_linker ];
	then
		mkdir -p $outdir/1_check_linker
        fi
        
	cd $outdir/1_check_linker
	
	echo "** check linker start ** @ "`date`  >> $outdir/${sample}.checked.txt
	
	echo "Starting to process $fq2:" >> check_linker.log
	echo "Raw reads in total:" >> check_linker.log
	zcat $fq2 | awk '{if (NR%4==2){print $0}}' - | wc -l >> check_linker.log
	# check the barcode-linker
	# Attention: you may need to adjust the linker cut site!
	echo "Proper barcode2-linker:" >> check_linker.log
	zcat $fq2 | awk '{if (NR%4==2){print $0}}' - | cut -c 9-23 | sort - | uniq -c | sed 's/^[ ]*//g' | sort -k1nr | head >> check_linker.log
	echo "Proper barcode1-linker:" >> check_linker.log
	zcat $fq2 | awk '{if (NR%4==2){print $0}}' - | cut -c 24-38,47-61 | sort - | uniq -c | sed 's/^[ ]*//g' | sort -k1nr | head >> check_linker.log
	# check multiple-same-linkers
	echo "Multiple barcode2-linker:" >> check_linker.log
	zcat $fq2 | grep "${barcode2_linker}[ATCG]\{22,24\}${barcode2_linker}" | wc -l >> check_linker.log
	echo "Multiple barcode1-linker:" >> check_linker.log
	zcat $fq2 | grep "${barcode1_linker}[ATCG]\{22,24\}${barcode1_linker}" | wc -l >> check_linker.log
	echo "Multiple barcode2-linker and misreplaced by ME:" >> check_linker.log
	zcat $fq2 | grep "${barcode2_linker}[ATCG]\{26,28\}${barcode2_linker}" | wc -l >> check_linker.log
	
	zcat $fq2 | grep "${barcode2_linker}[ATCG]\{22,24\}${barcode2_linker}" > R2_multi_barcode2_linker.txt
	zcat $fq2 | grep "${barcode1_linker}[ATCG]\{22,24\}${barcode1_linker}" > R2_multi_barcode1_linker.txt
	zcat $fq2 | grep "${barcode2_linker}[ATCG]\{27,29\}${barcode2_linker}" > R2_multi_barcode2_linker_mis_ME.txt
	
	echo "** check linker done ** @ "`date`  >> $outdir/${sample}.checked.txt
fi

###########################################################################################################
####---------------------------------------------cutadapt----------------------------------------------####
###########################################################################################################

if [ $run_cutadapt == "run" ];
then
	if [ ! -d $outdir/2_cutadapt ];
    	then
    		mkdir -p $outdir/2_cutadapt
	fi

	cd $outdir/2_cutadapt
	
	echo "** cutadapt fq start ** @ "`date`  >> $outdir/${sample}.checked.txt
	
	$cutadapt \
	-j $ppn \
	-a $ME \
	-m 37 \
	--trim-n \
	--pair-filter=first \
	-o ${sample}.cutadapt.R1.fastq.gz \
	-p ${sample}.cutadapt.R2.fastq.gz \
	--too-short-output ${sample}.cutadapt.short.R1.fastq.gz \
	--too-short-paired-output ${sample}.cutadapt.short.R2.fastq.gz \
	$fq1 $fq2 1>>${outdir}/${sample}_log.mark 2>&1
	
	echo "** cutadapt fq done ** @ "`date`  >> $outdir/${sample}.checked.txt
	
	if [ $run_fastqc_every_step == "run" ];
	then
	echo "** fastqc cutadapt fq start ** @ "`date`  >> $outdir/${sample}.checked.txt
	
	time $fastqc --java $java_jre -t $ppn -o $outdir/fastqc_res ${sample}.cutadapt.R1.fq.gz ${sample}.cutadapt.R2.fq.gz --noextract 1>>$outdir/${sample}_log.mark 2>&1
	
	echo "** fastqc cutadapt fq done ** @ "`date`  >> $outdir/${sample}.checked.txt
	fi
fi
          
if [ -f "${outdir}/2_cutadapt/${sample}.cutadapt.R1.fastq.gz" ];
then
	fq1=${outdir}/2_cutadapt/${sample}.cutadapt.R1.fastq.gz
	fq2=${outdir}/2_cutadapt/${sample}.cutadapt.R2.fastq.gz
fi

###########################################################################################################
####---------------------------------------------debarcode---------------------------------------------####
###########################################################################################################

if [ $run_debarcode == "run" ];
then
	if [ ! -d $outdir/3_umi_tools ];
        then
        	mkdir -p $outdir/3_umi_tools
	fi
	
	cd $outdir/3_umi_tools
	
	echo "** umi_tools fq start ** @ "`date`  >> $outdir/${sample}.checked.txt
	
	$umi_tools extract \
	--extract-method=regex \
	--bc-pattern2="^.{8}(?P<umi_1>.)(?P<discard_1>${debarcode_discard1}){e<=3}.{8}(?P<discard_2>${debarcode_discard2}){e<=3}.{10}(?P<discard_3>.*)$" \
	-I $fq1 -S ${sample}.extract.R1.fastq.gz \
	--read2-in=$fq2 --read2-out=${sample}.extract.R2.fastq.gz \
	-L ${sample}.umitools.debarcode.log
	
	echo "** umi_tools fq done ** @ "`date`  >> $outdir/${sample}.checked.txt
	
	if [ $run_fastqc_every_step == "run" ];
	then
	echo "** fastqc umi_tools fq start ** @ "`date`  >> $outdir/${sample}.checked.txt
	
	time $fastqc --java $java_jre -t $ppn -o $outdir/fastqc_res ${sample}.extract.R1.fastq.gz ${sample}.extract.R2.fastq.gz --noextract 1>>$outdir/${sample}_log.mark 2>&1
	
	echo "** fastqc umi_tools fq done ** @ "`date`  >> $outdir/${sample}.checked.txt
	fi
fi
	
if [ -f "${outdir}/3_umi_tools/${sample}.extract.R1.fastq.gz" ];
then
	fq1=${outdir}/3_umi_tools/${sample}.extract.R1.fastq.gz
	fq2=${outdir}/3_umi_tools/${sample}.extract.R2.fastq.gz
fi

###########################################################################################################
####------------------------------------Visualize Debarcoded Result------------------------------------####
###########################################################################################################

if [ $run_debarcoded_stat_visualize == "run" ];
then
	echo "** stat/visualize debarcoded results start ** @ "`date`  >> $outdir/${sample}.checked.txt
	
	zcat $fq2 | awk 'NR%4==2{print substr($1,1,16)}' | sed 's/^\([A-Z]\{8\}\)\([A-Z]\{8\}\).*$/\1_\2/g' - > ${sample}.extract.R2.barcodes.txt
	# Attention: Actually debarcode produces all the pixels within a permit-rule for downstream-analyzing, but we only visualize the exactly correct pixels!
	$python_exec $script_path/check_barcode_visualize.py \
	--bc=${sample}.extract.R2.barcodes.txt \
	--refA=$barcodes_A \
	--refB=$barcodes_B \
	--fig=${sample}.debarcoded_reads_heatmap.png \
	--map_stat=${sample}.debarcoded_passed_reads_stat.csv \
	--unmap_stat=${sample}.debarcoded_failed_reads_stat.csv \
	--barcodeNum=$npixel 1>>${outdir}/${sample}_log.mark 2>&1
	
	$r_exec $script_path/draw_svg.R \
	-f ${sample}.debarcoded_passed_reads_stat.csv \
	-s ${sample} \
	-p ./ 1>>${outdir}/${sample}_log.mark 2>&1
	
	echo "** stat/visualize debarcoded results done ** @ "`date`  >> $outdir/${sample}.checked.txt
fi

###########################################################################################################
####------------------------------------------------zUMI-----------------------------------------------####
###########################################################################################################

if [ $run_zumi == "run" ];
then
	if [ ! -d $outdir/4_zUMIs ];
    	then
    		mkdir -p $outdir/4_zUMIs
	fi
	
	cd ${outdir}/4_zUMIs
	
	echo "** zUMIs start ** @ "`date`  >> $outdir/${sample}.checked.txt
	
	if [ $zUMIs_yaml == "NULL" ];
	then
		if [ $yaml_default == "run_0.3" ];
		then
			cp $zUMIs_yaml_default1 zUMIs_${sample}.yaml
		else
			cp $zUMIs_yaml_default zUMIs_${sample}.yaml
		fi
			
		sed -i "s/sample_temp/$sample/" zUMIs_${sample}.yaml # s/ A/B /, if match one then next line
		sed -i "s/ppn_N/$ppn/" zUMIs_${sample}.yaml
		sed -i "s/mem_N/$memory/" zUMIs_${sample}.yaml

		outpath=$outdir/4_zUMIs
		temp=${outpath////;} # change / to ;
		sed -i "s/out_path/$temp/" zUMIs_${sample}.yaml
		sed -i "s/;/\//g" zUMIs_${sample}.yaml # s/ ;/\/ /g, change ; to / in the whole file

		temp=${fq1////;}
		sed -i "s/shareseq_temp_fq1/$temp/" zUMIs_${sample}.yaml
		sed -i "s/;/\//g" zUMIs_${sample}.yaml
		
		temp=${fq2////;}
		sed -i "s/shareseq_temp_fq2/$temp/" zUMIs_${sample}.yaml
		sed -i "s/;/\//g" zUMIs_${sample}.yaml
		
		temp=${STAR_idx////;}
		sed -i "s/star_index_path/$temp/" zUMIs_${sample}.yaml
		sed -i "s/;/\//g" zUMIs_${sample}.yaml

		temp=${ref_gtf////;}
		sed -i "s/ref_gtf/$temp/" zUMIs_${sample}.yaml
		sed -i "s/;/\//g" zUMIs_${sample}.yaml
		
		temp=${barcode_files////;}
		sed -i "s/barcode_files/$temp/" zUMIs_${sample}.yaml
		sed -i "s/;/\//g" zUMIs_${sample}.yaml
		
		if [ $run_zumi_from == "Filtering" ] || [ $run_zumi_from == "Mapping" ] || [ $run_zumi_from == "Counting" ] || [ $run_zumi_from == "Summarising" ];
		then
			sed -i "s/which_Stage: Filtering/which_Stage: $run_zumi_from/" zUMIs_${sample}.yaml
		else
			echo "$run_zumi_from is not in (Mapping, Counting, Summarising)!"
			exit
		fi
		
		$zUMIs -y zUMIs_${sample}.yaml 1>>${outdir}/${sample}_log.mark 2>&1
	elif [ -f "$zUMIs_yaml" ];
	then
		$zUMIs -y $zUMIs_yaml 1>>${outdir}/${sample}_log.mark 2>&1
	else
		echo "$zUMIs_yaml does not exist!"
		exit
	fi
	
	echo "** zUMIs done ** @ "`date`  >> $outdir/${sample}.checked.txt
fi

###########################################################################################################
####---------------------------------------------Final Stat--------------------------------------------####
###########################################################################################################
if [[ $run_final_stat == "run" ]]
then
    if [ ! -d $outdir/5_stat ];
    	then
    		mkdir -p $outdir/5_stat
	fi
	
	cd ${outdir}/5_stat
	
    ### bam ###
    # <sample>.filtered.tagged.unmapped.bam: fq filtered, base quality & BC & UMI (STAR-mapping input)
    
    # <sample>.filtered.tagged.Aligned.out.bam: STAR-mapping-to-genome (FeatureCount-counting input)
    # <sample>.filtered.tagged.Aligned.toTranscriptome.out.bam: STAR-mapping-to-transcriptome
    
    # <sample>.filtered.Aligned.GeneTagged.sorted.bam: FeatureCount-counting and coordinate-sorting

	# <sample>.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam: UMI-correcting (UX: raw, UB: corrected)
	
	### stat ###
	# <sample>.genecounts.txt: nFeature per cell
	# <sample>.UMIcounts.txt: nCount per cell
	# <sample>.readspercell.txt: reads per cell
	# <sample>.reads_per_gene.txt: reads per gene
	
	echo "** stat start ** @ "`date`  >> $outdir/${sample}.checked.txt
	
	echo "## check linker ##" > ${sample}.stat_full_pipeline.txt
	cat $outdir/1_check_linker/check_linker.log >> ${sample}.stat_full_pipeline.txt
	
	echo "## cDNA-ME-trimmed from cutadapt ##" >> ${sample}.stat_full_pipeline.txt
	expr `zcat ${outdir}/2_cutadapt/${sample}.cutadapt.R1.fastq.gz | wc -l` / 4 >> ${sample}.stat_full_pipeline.txt
	
	echo "## debarcoded from umi_tools ##" >> ${sample}.stat_full_pipeline.txt
	cat ${outdir}/3_umi_tools/${sample}.extract.R2.barcodes.txt | wc -l >> ${sample}.stat_full_pipeline.txt
	
	echo "## mapped from STAR ##" >> ${sample}.stat_full_pipeline.txt
	samtools view -b -h -F 4 ${outdir}/4_zUMIs/${sample}.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam > ${sample}.mapped.bam
	samtools view -c ${sample}.mapped.bam >> ${sample}.stat_full_pipeline.txt
	
	echo "## unique mapped from STAR ##" >> ${sample}.stat_full_pipeline.txt
	samtools view -b -h -F 4 -q 255 ${sample}.mapped.bam > ${sample}.mapped.unique.bam
	samtools view -c ${sample}.mapped.unique.bam >> ${sample}.stat_full_pipeline.txt
	
	echo "## annotated to gtf from FeatureCount ##" >> ${sample}.stat_full_pipeline.txt
	echo "# total annotated reads #" >> ${sample}.stat_full_pipeline.txt
	samtools view ${sample}.mapped.bam | grep "Assigned3" | wc -l >> ${sample}.stat_full_pipeline.txt
	echo "# distribution of annotated reads #" >> ${sample}.stat_full_pipeline.txt
	# 两次FeatureCount的结果(exon+intron)加起来不等于Assigned reads的总数, 因为一条read同时注释到exon和intron上只+1
	samtools view ${sample}.mapped.bam | grep "ES:Z:.*IS:Z:\w*" -o | awk -F ' ' '{print $1"   "$NF}' | sort | uniq -c >> ${sample}.stat_full_pipeline.txt
	echo "## demultiplexed from taggd (fron matrix.RDS) ##" >> ${sample}.stat_full_pipeline.txt
	$r_exec $script_path/stat.R \
	-f ${outdir}/4_zUMIs/zUMIs_output/expression/${sample}.dgecounts.rds \
	-p $outdir/5_stat/${sample}.stat_full_pipeline.txt 1>${sample}_stat_R_log.txt 2>&1
	cat ${sample}_stat_R_log.txt >> $outdir/${sample}_log.mark
	#cat ${outdir}/4_zUMIs/zUMIs_output/stats/${sample}.readspercell.txt | sed '1d' | awk -F ' ' '{x[$3] += $2} END {for(i in x){print i, x[i]}}' >> ${sample}.stat_full_pipeline.txt
	#cat ${outdir}/4_zUMIs/zUMIs_output/stats/${sample}.UMIcounts.txt | sed '1d' | awk -F ' ' '{x[$3] += $1} END {for(i in x){print i, x[i]}}' >> ${sample}.stat_full_pipeline.txt
	
    echo "** stat done ** @ "`date`  >> $outdir/${sample}.checked.txt
fi
