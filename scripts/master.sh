#!/usr/bin/bash -l
#$ -cwd
#$ -q all.q
#$ -o stdoutput.out
#$ -e error.out
#$ -m eetu.eklund@maryland.gov
### add more options above to set node, increase memory...
# These options will be used if the program is run with a qsub command (which is not available currently)

################################################################################
# Script Name: master.sh
# Description: QC Analysis and Sequence Alignment Script.
#               Script currently works in the same directory as fastq read files
#               Performs QC, BWA Alignment, ANI calculations, and more on read files
# Author: Eetu Eklund
# Email: eetu.eklund@maryland.gov
# Execute Script: bash master.sh YYMMDD runAMR
# Examples: bash scripts/master.sh 230101
# 	    bash scripts/master.sh 230101 runAMR  <--option to run AMRFinderPlus
# Works with read files in current working directory
################################################################################

# Future additions:
# docker containerization!!
# Error handling!!
# add argument flags!!
# kraken 2 as contamination check and double check for sendsketch. DONE
# phylogenetic trees? This should probably be done manually after with the fastq files and some Reference fasta files with outgroups
# BUSCO & MIDAS (Busco only works on .fasta)


#1. Reference Preparation and indexing
	# The following operations are done in a for loop through all fastq files.
	# fastp to filter read files and QC
	# sendsketch matches current fastq file to a reference sequence and gets ANI
	# ncbi-genome-download downloads the reference
	# bwa index, samtool faidx, and picard dict are run on the reference sequences
	# bwa mem alignments are done using the current pair of read files
	# samtools used for BAM creation and sorting
	# optionally runs AMRFinderPlus and creates .vcf and consensus fasta files

# 2. Sample QC and cleaning
	# fastqc is run on all fastq files and sorted bam files
	# confindr determines levels of contamination in read files
	# qualimap multi-bamqc is performed on all sorted bam files
	# multiqc combines QC metrics into an html
	# write_summary.py writes and excel file with key QC metrics for all fastq files
	# directory is cleaned and unnecesary/large files are deleted or zipped.


# check if condaenv exists. Activates condaenv if it does, creates condaenv and installs requirements if it does not.
if conda env list | grep ".*fastq_analysis.*"; then
        source ~/miniconda3-2/etc/profile.d/conda.sh
        conda activate fastq_analysis
else
        source ~/miniconda3-2/etc/profile.d/conda.sh
        conda create -n fastq_analysis python=3.9 anaconda
	conda init bash
	conda activate fastq_analysis
        conda install -c conda-forge bwa samtools biokit fastqc multiqc qualimap ncbi-genome-download picard xlsxwriter bcftools fastp=0.22.0 #confindr seqkit seqkt mummer
	pip install confindr #<-- in the conda env. A python file in the confindr has to be modified to work
fi


# These inputs should be changed to flags (-o, --output, -AMR, -BMGAP)
dir=$1
if [ -z "$dir" ]; then
	echo Please provide a directory name as an argument.
	exit
fi
mkdir ./$dir
mkdir ./$dir/results

runAMR=$2
runBMGAP=$3
BMGAP=False
if [[ $runAMR == "runAMR" ]]; then
	echo
	echo
        echo Running program with AMRFinderPlus to find antimicrobial resistance genes in samples.
        echo This will at least double the compute time...
	echo
	echo
	sleep 1s
	mkdir $dir/results/AMRFinder
	mkdir $dir/results/consensus_seqs
	mkdir $dir/results/vcf_files
	if [[ $runBMGAP == "runBMGAP" ]]; then
		mkdir $dir/results/BMGAP_output
		echo Running BMGAP to serotype Meningitidis sequences.
		echo This will only work with Meningitidis read files and will increase compute time.
		BMGAP=True
	fi
fi

#if [ -d "./scripts/" ]; then
#        mv ./scripts/* .
#fi
# creates directories to put output files in.
mkdir $dir/reference
mkdir $dir/reads
mkdir $dir/results/alignments
mkdir $dir/results/sendsketch
mkdir $dir/results/kraken2
if [ -d "./reads" ]; then
        mv ./reads/* .
fi

# export paths for MIDAS
export PYTHONPATH=/home/mdhsequencing/.local/lib/python3.8:~/7T/MIDAS
export PATH=$PATH:~/7T/MIDAS/scripts
export MIDAS_DB=~/7T/MIDAS/MIDAS/midas_db_v1.2


# export gsl path
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/mdhsequencing/gsl/lib
export LD_LIBRARY_PATH

# for loop performs these operations for each pair of fastq file at a time
        # get sample id information and other read group info
	# fastp performs QC and filtering
        # sendsketch to match it to reference and get ANI
	# ncbi-genome-download downloads reference sequence based on taxid from sendsketch
        # BWA and samtools index each reference
        # Picard CreateSequenceDictionary makes dict file for bwa
        # bwa mem for all fastq pairs and their reference
	# samtools creates sorted bam file
	# optionally: AMRFinderPlus is used to find antimicrobial resistance genes
	#	      This also creates a vcf (variant) file and a consensus sequence file based on the bam

# initialize variables and make directory to be used in for loop
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
		filename="${R1File%%.*}"
		filename2="${file%%.*}"

		# Add fastp here for quality control, filtering, and adapter trimming
		fastp -i $R1File -I $file -o "${filename}_filtered.fastq.gz" -O "${filename2}_filtered.fastq.gz" \
		-j $filename.json -h $filename.html
		# Remove fastq files once we have filtered fastq files
		rm $R1File
		rm $file

		# runs kraken2 on filtered fastq files and extracts top 5 results at species level
		bash run_classifier.sh "${filename}_filtered.fastq.gz" "${filename2}_filtered.fastq.gz"


		# sendsketch finds closest matching taxon and its taxid
		~/bbmap/sendsketch.sh "${filename}_filtered.fastq.gz" reads=1m samplerate=0.5 minkeycount=2 > sendsketch.txt


		# MIDAS
		conda deactivate
		conda activate py2
		run_midas.py species $dir/results/MIDAS/$filename -1 "${filename}_filtered.fastq.gz" \
		-2 "${filename2}_filtered.fastq.gz" -t 12 -d ~/7T/MIDAS/MIDAS/midas_db_v1.2 --remove_temp
		conda deactivate
		conda activate fastq_analysis


		# could use kraken2 as a double check on sendsketch and confindr
		# would tell us taxid of closest genus+species if ncbi-genome-download fails with taxid from sendsketch
		# would find % of other species' sequences found in sample and alert if that % is too high
		# Takes a while to load into memory, has to do it for each fastq pair

		# extracts taxid for the species (only first single reference file)
		taxids=$(python3 scripts/extract_taxid.py)

		# if taxid has not been used, reference file is downloaded and index/dict files created
		if ! echo ${taxid_list[@]} | grep -w -q $taxids; then
			# downloads reference file
			# uses second taxid if first does not have a downloadable file
			# reference sequences are not ideal but work fine
			ncbi-genome-download --taxids $taxids --formats fasta --parallel 4 bacteria,viral || \
			(sed '1d' sendsketch.txt > tmpfile; mv tmpfile sendsketch.txt;taxids=$(python3 \
			scripts/extract_taxid.py);ncbi-genome-download --taxids $taxids --formats fasta --parallel 4 \
			bacteria,viral; sed -i '1s/^/\n/' sendsketch.txt)

			# keeps list of taxids used and makes a directory named after taxid for reference
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
		bwa mem -R '@RG\tID:'$ID'\tLB:'$LB'\tPL:'$PL'\tSM:'$SM'\tPU:'$PU -v 2 -t 8 tmpRef/*.fna \
		"${filename}_filtered.fastq.gz" "${filename2}_filtered.fastq.gz" > $dir/results/alignments/$filename.sam

		# samtools indexing and making bam files
		samtools view -u -@ 8 $dir/results/alignments/$filename.sam | samtools sort -@ 8 \
		-o $dir/results/alignments/"${filename}_sorted.bam"
		samtools index -@ 8 $dir/results/alignments/"${filename}_sorted.bam"

		#Mark pcr duplicates in the sorted bam file before making vcf files for AMRFinder

		# if runAMR option is given, AMRFinderPlus is used to find antimicrobial resistance genes
		# This will increase runtime a good amount, so it is optional
		if [[ $runAMR == "runAMR" ]]; then
			bash runAMR.sh $dir/results/alignments/"${filename}_sorted.bam" tmpRef/*.fna $filename
			mv *.vcf.gz* $dir/results/vcf_files
			mv *_consensus.fasta $dir/results/consensus_seqs
			mv *_AMRout.txt $dir/results/AMRFinder
		fi

		# moves all reference files to reference directory once we are done with them
		mv sendsketch.txt $dir/results/sendsketch/"${filename}_sendsketch.txt"
		mv ./tmpRef/* $dir/reference/$taxids/.
		bothfiles=false
		# continue skips rest of for loop so bothfiles stays false
		continue
	fi
    	R1File=$file
    	bothfiles=true
done
rm -r ./tmpRef
mv top_kraken_results_*.txt $dir/results/kraken2
mv *k2_report.txt $dir/results/kraken2


# THIS FASTQC could be deleted when fastp has been tested and works well, but summary file needs fastqc results
# if we want to delete this fastqc, write_summary.py has to be edited to find QC metrics from fastp
# QC general statistics (pass/fail) are not in the results of fastp
# Quality Control with FastQc
fastqc *_filtered.fastq.gz -o $dir -t 12
#fastqc is set temporarily to a different directory so multiqc doesnt use it (so summary file gets written properly)
#fastqc $dir/results/alignments/*.sorted.bam -o $dir/results -t 12


# Include these in our path for confindr to use
#export PATH=$PATH:/home/mdhsequencing/bbmap
#export PATH=$PATH:/home/mdhsequencing/kma

# contamination finder
# Not easy to get working. had to change a line of code in two python files (database_setup.py, bbtools.py)
# :~/miniconda3-2/envs/condaenv/lib/python3.9/site-packages/confindr_src
# confindr_src/wrappers/bbtools.py ----> in bbduk_bait function, added --ignorejunk flag to bbduk.sh cmd = ... commands
# confindr_src/database_setup.py ----> line 209, added .encode() to the end of the line
#confindr.py -i ./ -o $dir/results/confindr_out --rmlst -d ~/.confindr_db


# Qualimap mapping quality report
unset display # stops qualimap display from appearing

# uses bamqc if only one pair of fastq files provided, else it uses multi-bamqc
if [[ $(ls -lR ./*.fastq.gz | wc -l) == 2 ]]; then
	qualimap bamqc -bam $dir/results/alignments/*.bam --java-mem-size=5G -outdir $dir/results/qualimap
else
	# Python program to write bam_files.txt (file of bam file names use by multi-bamqc)
	python3 scripts/write_bam_files.py $dir
	mv bam_files.txt $dir/results/alignments
	# qualimap for multiple bam files
	qualimap multi-bamqc -r -d $dir/results/alignments/bam_files.txt --java-mem-size=5G -outdir $dir/results/qualimap
fi

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


# Writes all quality control statistics, Coverage, ANI... from read files to one excel file
python scripts/write_summary.py $dir

# move bam fastqc file and others back after summary written
mv $dir/results/fastqc/tmp/* $dir/results/fastqc

# writes summary file from AMRFinderPlus results
if [[ $runAMR == "runAMR" ]]; then
	python3 scripts/AMR_write.py $dir
fi

# Remove large files (like sam files but keep bam) to save space
# zip fastq files and other large files
echo
echo "Deleting, compressing, and moving files to clean directory and conserve storage space..."
echo

rm $dir/results/alignments/*.sam
rm -r $dir/results/fastqc/tmp

ref_dirs=$(ls $dir/reference/)
for ref_dir in $ref_dirs; do
	gzip $dir/reference/$ref_dir/*.fna
	gzip $dir/reference/$ref_dir/*.fna.bwt
done
mv *.fastq* $dir/reads
if [ ! -d "./scripts/" ]; then
        mkdir scripts
fi
#mv *.sh scripts
#mv *.py scripts
mv *.html $dir/results/fastqc
mv *.json $dir/results/fastqc

### Delete unfiltered fastq files to conserve space
#shopt -s extglob
#rm $dir/reads/!(*filtered*)

conda deactivate

# BMGAP is only for serotyping meningitis sequences. this can be removed or added as an argument flag with runAMR
if [[ $BMGAP == True ]]; then
	source ~/7T/BMGAP/pipeline/bmgap/bin/activate
	echo
	echo Running BMGAP to serotype Meningitis sequences.
	echo
	python3 ~/7T/BMGAP/pipeline/PMGA/blast_pubmlst.py -d $dir/results/consensus_seqs -o $dir/results/BMGAP_output -t 4 -sg -fr -p
fi
if [ -d "$dir/results/consensus_seqs" ]; then
        gzip $dir/results/consensus_seqs/*
fi

echo
echo "Pipeline Completed Succesfully!"
