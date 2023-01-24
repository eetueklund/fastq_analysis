#!/usr/bin/bash -l
#$ -cwd
#$ -q all.q
#$ -o stdoutput.out
#$ -e error.out
#$ -m eetu.eklund@maryland.gov
### add more options above to set node, increase memory...
# These options will be used if the program is run with a qsub command (which does not work)

################################################################################
# Script Name: master.sh
# Description: QC Analysis and Sequence Alignment Script.
#               Script currently works in the same directory as fastq read files
#               Performs QC, BWA Alignment, and ANI calculation on read files
# Author: Eetu Eklund
# Email: eetu.eklund@maryland.gov
# Execute Script: bash master.sh
################################################################################

# Future additions:
# Error handling

# 1. Reference Preparation and indexing
	# Done in a for loop through all fastq files. 
	# sendsketch matches current fastq file to a reference sequence and gets ANI
	# ncbi-genome-download downloads the reference
	# bwa index, samtool faidx, and picard dict are run on the reference sequences
	# bwa mem alignments are done using the current pair of read files
	# samtools used for BAM creation and sorting

# 2. Sample QC and cleaning 
	# fastqc is run on all fastq files and sorted bam files
	# confindr determines levels of contamination in read files
	# qualimap multi-bamqc is performed on all sorted bam files
	# multiqc combines QC metrics into an html
	# write_summary.py writes and excel file with key QC metrics for all fastq files
	# directory is cleaned and unnecesary/large files are deleted or zipped.


# check if condaenv exists. Activates condaenv if it does, creates condaenv and installs requirements if it does not.
if conda env list | grep ".*condaenv.*"; then
        source ~/miniconda3-2/etc/profile.d/conda.sh
        conda activate condaenv
else
        source ~/miniconda3-2/etc/profile.d/conda.sh
        conda create -n condaenv python=3.9 anaconda
	conda init bash
	conda activate condaenv
        conda install -c conda-forge bwa samtools biokit fastqc multiqc qualimap ncbi-genome-download picard xlsxwriter bcftools #seqkit seqkt mummer
fi


# Program takes in one argument: The current date to name directory where all files are put in 
# date (YY/MM/DD)
dir=$1
if [ -z "$dir" ]; then
	echo Please provide a directory name as an argument. 
	exit
fi
mkdir ./$dir
if [ -d "./scripts/" ]; then
        mv ./scripts/* .
fi
# creates directories to put output files in.
mkdir $dir/results
mkdir $dir/reference
mkdir $dir/results/alignments
mkdir $dir/results/sendsketch
if [ -d "./reads" ]; then
        mv ./reads/* . 
fi


LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/mdhsequencing/gsl/lib
export LD_LIBRARY_PATH
# for loop performs these operations for each pair of fastq file at a time
        # get sample id information and other read group info
        # sendsketch to match it to reference and get ANI
        # BWA and samtools index for each reference
        # Picard CreateSequenceDictionary
        # bwa mem for all fastq pairs and their reference (takes a while)
R1File=''
taxid_list=()
bothfiles=false
mkdir tmpRef
for file in *.fastq*; do
	
	# Gets Read Group information to include for BWA run. (ID,SM,LB,PU,PL)
	front="${file%%_*}"
        front2="${front#*-}"
        front3="${front2#*-}"
        back="${file#*_}"
        back2="${back#*_}"
        ID="${file%%-*}"
        SM="${front3#*-}"
        LB="${back%%_*}"
        PU="${back2%%_*}"
        PL=illumina

    	# checks if we have both read files (and that they are pairs)
    	if [ "$bothfiles" = true ]; then
		# sendsketch finds closest matching taxon and its taxid
		~/bbmap/sendsketch.sh $R1File reads=1m samplerate=0.5 minkeycount=2 > sendsketch.txt
		
		# extracts taxid for the species (only first single reference file)
		taxids=$(python3 extract_taxid2.py)
		# if taxid has not been used, reference file is downloaded and index/dict files created
		if ! echo ${taxid_list[@]} | grep -w -q $taxids; then
			# downloads reference file
			ncbi-genome-download --taxids $taxids --formats fasta --parallel 4 bacteria,viral
			taxid_list=(${taxid_list[@]} $taxids)
			refdirname=$(find ./refseq/bacteria/ -name 'GCF_*' | head -n 1)
                	mv $refdirname/*fna* ./tmpRef
                	gunzip ./tmpRef/*fna*
			mkdir $dir/reference/$taxids
			# gets name of reference sequence to be used in CreateSequenceDictionary
                	ref=$(basename ./tmpRef/*.fna)
                	# creates index files, .dict, and .fai.
                	bwa index ./tmpRef/*
                	samtools faidx ./tmpRef/*.fna
                	picard CreateSequenceDictionary R=./tmpRef/$ref >/dev/null 2>&1
			rm -r refseq
		else # if taxid has been used before, move reference and index files back to tmpRef
			mv $dir/reference/$taxids/* ./tmpRef/.
		fi

		# bwa alignment
		#bwa mem -v 2 -t 4 tmpRef/*.fna $R1File $file > results/alignments/$R1File.sam
		bwa mem -R '@RG\tID:'$ID'\tLB:'$LB'\tPL:'$PL'\tSM:'$SM'\tPU:'$PU -v 2 -t 8 tmpRef/*.fna $R1File $file > $dir/results/alignments/$R1File.sam

		# samtools indexing and making bam files
		samtools view -@ 8 $dir/results/alignments/$R1File.sam -o $dir/results/alignments/$R1File.bam
		samtools sort -@ 8 $dir/results/alignments/$R1File.bam -o $dir/results/alignments/$R1File.sorted.bam
		samtools index -@ 8 $dir/results/alignments/$R1File.sorted.bam

		# moves all reference files to reference directory once we are done with them
		mv sendsketch.txt $dir/results/sendsketch/$R1File.sendsketch.txt
		mv ./tmpRef/* $dir/reference/$taxids/.
		bothfiles=false
		# continue skips rest of for loop
		continue
	fi 
    	R1File=$file
    	bothfiles=true
done
rm -r ./tmpRef


# Quality Control with FastQc
fastqc *.fastq* -o $dir -t 12
fastqc $dir/results/alignments/*.sorted.bam -o $dir -t 12

# Include these in our path for confindr to use
export PATH=$PATH:/home/mdhsequencing/bbmap
export PATH=$PATH:/home/mdhsequencing/kma

# contamination finder
# Not easy to get working. had to change a line of code in two python files (database_setup.py, bbtools.py)
# :~/miniconda3-2/envs/condaenv/lib/python3.9/site-packages/confindr_src
# confindr_src/wrappers/bbtools.py ----> in bbduk_bait function, added --ignorejunk flag to bbduk.sh commands
# confindr_src/database_setup.py ----> line 209, added .encode() to the end of the line
confindr.py -i ./ -o $dir/results/confindr_out --rmlst -d ~/.confindr_db


# Qualimap mapping quality report
unset display # stops qualimap display from appearing
#qualimap bamqc -bam results/alignments/*.sorted.bam --java-mem-size=5G -outdir ./results/qualimap

# Python program to write bam_files.txt (file of bam file names use by multi-bamqc)
python3 write_bam_files2.py $dir
mv bam_files.txt $dir/results/alignments

# qualimap for multiple bam files
qualimap multi-bamqc -r -d $dir/results/alignments/bam_files.txt --java-mem-size=5G -outdir $dir/results/qualimap


# multiqc scans current directory and makes a report from fastqc outputs
multiqc $dir

# Move files to clean working directory
mkdir $dir/results/fastqc
mv multiqc_data $dir/results
mv multiqc_report.html $dir/results


if [ ! -f "$dir/*fastqc*" ]; then
        mv $dir/*fastqc* $dir/results/fastqc
fi

mv  $dir/results/alignments/*_stats $dir/results/qualimap


# unzip and move fastqc files for write_summary.py to be able to extract and calculate mean quality scores
mkdir $dir/results/fastqc/tmp
for file in $dir/results/fastqc/*fastqc.zip; do
	unzip $file -d $dir/results/fastqc/tmp >/dev/null 2>&1
done
rm -r $dir/results/fastqc/tmp/*.sorted_fastqc

# Writes all quality control statistics, Coverage, ANI... from read files to one excel file
python write_summary2.py $dir


# Remove large files (like sam files but keep bam) to save space
# zip fastq files and other large files
echo "Deleting, compressing, and moving files to clean directory and conserve storage space..."
rm $dir/results/alignments/*.sam
rm -r $dir/results/fastqc/tmp

ref_dirs=$(ls $dir/reference/)
for ref_dir in $ref_dirs; do
	if [ -f "$dir/reference/ref_dir/*.fna" ]; then
		gzip $dir/reference/*/*.fna
	fi
	if [ -f "$dir/reference/ref_dir/*.fa" ]; then
		gzip $dir/reference/*/*.fa
	fi
done

if [ ! -d "$dir/reads" ]; then
	mkdir $dir/reads
	mv *.fastq* $dir/reads
fi
if [ ! -d "./scripts/" ]; then
        mkdir scripts
fi
mv *.sh scripts
mv *.py scripts



conda deactivate
