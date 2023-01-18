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
# confindr - finds intraspecies contamination in raw illumina data

# 1. Reference Preparation and indexing
	# Done in a for loop through all fastq files. 
	# Kraken2 matches current fastq file to a reference sequence and downloads it
	# bwa index, samtool faidx, and picard dict are run on the reference sequences
	# bwa mem alignments are done using the current pair of read files
	# samtools used for BAM creation and sorting
	# bcftools to make a consensus sequence based on read files 
	# FastANI to calculate Average Nucleotide Identity

# 2. Sample QC and cleaning 
	# fastqc is run on all fastq files and sorted bam files
	# qualimap multi-bamqc is performed on all sorted bam files
	# multiqc combines QC metrics into an html
	# write_summary.py writes and excel file with key QC metrics for all fastq files
	# directory is cleaned and unnecesary/large files are deleted or zipped.


# check if condaenv exists. Activates condaenv if it does, creates condaenv and installs requirements if it does not.
if conda env list | grep ".*condaenv.*"; then
        source ~/miniconda3-2/etc/profile.d/conda.sh
	conda init bash
        conda activate condaenv
else
        source ~/miniconda3-2/etc/profile.d/conda.sh
        conda create -n condaenv python=3.9 anaconda
	conda init bash
	conda activate condaenv
        conda install -c conda-forge bwa samtools biokit fastqc multiqc faqcs qualimap ncbi-genome-download picard xlsxwriter bcftools #seqkit seqkt mummer
fi


# checks if directory exists, creates one if it does not
if [ ! -d "./results" ]; then
        mkdir results
fi
if [ ! -d "./reference/" ]; then
        mkdir reference
fi
if [ ! -d "./results/alignments" ]; then
        mkdir ./results/alignments
fi
if [ ! -d "./results/variants" ]; then
        mkdir ./results/variants
fi
if [ ! -d "./results/ANI" ]; then
        mkdir ./results/ANI
fi
if [ -d "./scripts/" ]; then
        mv ./scripts/* . 
fi
if [ ! -d "./results/K2" ]; then
        mkdir ./results/K2
fi
if [ -d "./reads" ]; then
        mv ./reads/* . 
fi

# FaQCs - Quality Control and trimming
#../FaQCs/FaQCs -p *.fastq* -d FaQCs_output
~/FaQCs/FaQCs -1 *R1*.fastq* -2 *R2*.fastq* -d results/FaQCs_output


if (file *fastq* | grep -q compressed ) ; then
     gunzip *fastq.gz
fi

# code to make fastANI work by finding gsl lib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/mdhsequencing/gsl/lib
export LD_LIBRARY_PATH

# for loop performs these operations for each pair of fastq file at a time
        # get sample id information and other read group info
        # Kraken2 to match it to reference
        # BWA index for each reference
        # Picard CreateSequenceDictionary
        # samtools index for each reference
        # bwa mem for all fastq pairs and their reference (takes a while)
        # bcftools used to create a consensus fasta to be used for calculating ANI
        # fastANI calculates Average Nucleotide Identity from reference and consensus sequences.
R1File=''
taxid_list=()
bothfiles=false
mkdir tmpRef
for file in *.fastq; do
	
	# Not the cleanest way to get Read Group information to include for BWA run.
    	ID1="${file%-*}"  # gives everything before the last -
        ID2="${ID1%-*}"
        ID3="${ID2%-*}"
        after="${file#*-}" # gives everything after the first -
        after2="${after#*-}"
        after3="${after2#*-}"
        after4="${after3#*-}"
        middle="${after3%%.*}" # gives everything before . and after -
        middle2="${middle%_*}"
        middle3="${middle2%_*}"
        middle4="${middle3%_*}"
	LB2="${middle3#*_}"

        ID="${ID3%-*}"
        SM="${middle4%_*}"
        PU="${middle4#*_}"
        LB="${LB2#*_}"
        PL=illumina
	
    	# checks if we have both read files (and that they are pairs)
    	if [ "$bothfiles" = true ]; then
		
		# Kraken2 finds which species the read pair is from
		kraken2 --db ~/kraken2/standard_db --paired --use-names --threads 4 --output k2_output.txt --report k2_report.txt --classified-out seqs#.fq $R1File $file
		
		# extracts taxid for the species (only first single reference file)
		taxids=$(python3 extract_taxid.py)
		# if taxid has not been used, reference file is downloaded and index/dict files created
		if ! echo $taxid_list | grep -w -q $taxids; then
			# downloads reference file
			ncbi-genome-download --taxids $taxids --formats fasta --parallel 4 bacteria,viral
			taxid_list=(${taxid_list[@]} $taxids)
			refdirname=$(find ./refseq/bacteria/ -name 'GCF_*' | head -n 1)
                	mv $refdirname/*fna* ./tmpRef
                	gunzip ./tmpRef/*fna*
			mkdir reference/$taxids
			# gets name of reference sequence to be used in CreateSequenceDictionary
                	ref=$(basename ./tmpRef/*.fna)
                	# creates index files, .dict, and .fai.
                	bwa index ./tmpRef/*
                	samtools faidx ./tmpRef/*.fna
                	picard CreateSequenceDictionary R=./tmpRef/$ref >/dev/null 2>&1
		else # if taxid has been used before, move reference and index files back to tmpRef
			ls reference/$taxids
			mv ./reference/$taxids/* ./tmpRef/.
		fi
		
		# Stores taxid in list. use ncbi-download only when taxid has not already been used.
                #taxid_list=(${taxid_list[@]} $taxids)
		
		# moves 1 reference file to temporary directory where bwa indexing will be performed
		#refdirname=$(find ./refseq/bacteria/ -name 'GCF_*' | head -n 1)
		#mv $refdirname/*fna* ./tmpRef
        	#gunzip ./tmpRef/*fna*

		# gets name of reference sequence to be used in CreateSequenceDictionary
		#ref=$(basename ./tmpRef/*.fna)
		
		# creates index files, .dict, and .fai. 
		#bwa index ./tmpRef/*
		#samtools faidx ./tmpRef/*.fna
		#picard CreateSequenceDictionary R=./tmpRef/$ref >/dev/null 2>&1

		# bwa alignment
		#bwa mem -v 2 -t 4 tmpRef/*.fna $R1File $file > results/alignments/$R1File.sam
		bwa mem -R '@RG\tID:'$ID'\tLB:'$LB'\tPL:'$PL'\tSM:'$SM'\tPU:'$PU -v 2 -t 4 tmpRef/*.fna $R1File $file > results/alignments/$R1File.sam

		# samtools indexing and making bam files
		samtools view ./results/alignments/$R1File.sam -o ./results/alignments/$R1File.bam
		samtools sort ./results/alignments/$R1File.bam -o ./results/alignments/$R1File.sorted.bam
		samtools index ./results/alignments/$R1File.sorted.bam

		# calculate ANI
		# make consensus fasta file from read files by running variant calling/filtering and 
		# making a consensus sequence
		bcftools mpileup -f ./tmpRef/$ref ./results/alignments/$R1File.sorted.bam | bcftools call -m --threads 4 -v -Ov -o ./results/variants/$R1File.variants.vcf.gz

		bcftools index ./results/variants/$R1File.variants.vcf.gz

		bcftools norm -f ./tmpRef/$ref ./results/variants/$R1File.variants.vcf.gz -Ob -o ./results/variants/$R1File.norm.bcf

		bcftools filter --IndelGap 5 ./results/variants/$R1File.norm.bcf -Ob -o ./results/variants/$R1File.norm.flt-indels.bcf

		cat ./tmpRef/$ref | bcftools consensus ./results/variants/$R1File.variants.vcf.gz > ./tmpRef/$R1File.consensus.fa

		# FastANI uses new fasta file and reference fasta to calculate average nucleotide identity
		~/FastANI/fastANI -q ./tmpRef/$R1File.consensus.fa -r ./tmpRef/$ref -o ./results/ANI/$R1File.fastani.out

		# moves all reference files to reference directory once we are done with them
		mv ./tmpRef/* ./reference/$taxids/.
		mv k2_report.txt ./results/K2/$R1File.k2_report.txt
		rm -r refseq
		bothfiles=false
		# continue skips rest of for loop
		continue
    	fi 
    	R1File=$file
    	bothfiles=true
done

rm -r ./tmpRef


# Include samtools coverage here? makes histogram or table of coverage
# also done in fastqc or qualimap
# samtools coverage bwa_output.sorted.bam


# Quality Control with FastQc
fastqc *.fastq
fastqc results/alignments/*.sorted.bam
rm seqs_*.fq

# Include these in our path for confindr to use
export PATH=$PATH:/home/mdhsequencing/bbmap
export PATH=$PATH:/home/mdhsequencing/kma

# contamination finder
# Not easy to get working. had to change a line of code in two python files (database_setup.py, bbtools.py)
# :~/miniconda3-2/envs/condaenv/lib/python3.9/site-packages/confindr_src
# confindr_src/wrappers/bbtools.py ----> in bbduk_bait function, added --ignorejunk flag to bbduk.sh commands
# confindr_src/database_setup.py ----> line 209, added .encode() to the end of the line
confindr.py -i ./ -o ./results/confindr_out --rmlst -d ~/.confindr_db


# Qualimap mapping quality report
unset display # stops qualimap display from appearing
#qualimap bamqc -bam results/alignments/*.sorted.bam --java-mem-size=5G -outdir ./results/qualimap

# Python program to write bam_files.txt (file of bam file names use by multi-bamqc)
python3 write_bam_files.py
mv bam_files.txt ./results/alignments

# qualimap for multiple bam files
qualimap multi-bamqc -r -d results/alignments/bam_files.txt --java-mem-size=5G -outdir ./results/qualimap


# multiqc scans current directory and makes a report from fastqc outputs
multiqc .

# Move files to clean working directory
mkdir ./results/fastqc
mv multiqc_data results
mv multiqc_report.html results


if [ ! -f "./*fastqc*" ]; then
        mv *fastqc* ./results/fastqc
fi
if [ ! -f "./results/alignments/*fastqc*" ]; then
        mv ./results/alignments/*fastqc* results/fastqc
fi
mv  results/alignments/*_stats results/qualimap


# unzip and move fastqc files for write_summary.py to be able to extract and calculate mean quality scores
mkdir results/fastqc/tmp
for file in results/fastqc/*fastqc.zip; do
	unzip $file -d ./results/fastqc/tmp >/dev/null 2>&1
done
rm -r results/fastqc/tmp/*.sorted_fastqc

# Writes all quality control statistics, Coverage, ANI... from read files to one excel file
python write_summary.py


# Remove large files (like sam files but keep bam) to save space
# zip fastq files and other large files
echo "Deleting, compressing, and moving files to clean directory and conserve storage space..."
rm results/alignments/*.sam
rm k2_output.txt
rm results/FaQCs_output/*.fastq*
rm -r results/fastqc/tmp
gzip reference/*/*.fna
gzip reference/*/*.fa
gzip *.fastq

if [ ! -d "./reads" ]; then
        mkdir reads
       mv *.fastq* ./reads
fi
if [ ! -d "./scripts/" ]; then
        mkdir scripts
        mv *.sh ./scripts
        mv *.py ./scripts
else
	mv *.sh ./scripts
        mv *.py ./scripts
fi

conda deactivate
