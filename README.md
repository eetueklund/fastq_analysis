# Allin1
paired fastq file all-in-one analysis script. Quality Control, Alignment, ANI, sendsketch, multiqc, confindr...

How to use: Have all fastq files either in current directory or ./reads directory.

        Command: bash master.sh
        
        # if scripts are in the ./scripts directory
        Command: bash scripts/master.sh

1. Reference Preparation and indexing

        Done in a for loop through all fastq files.
        sendsketch matches current fastq file to a reference sequence and gets ANI
        ncbi-genome-download downloads the reference
        bwa index, samtool faidx, and picard dict are run on the reference sequences
        bwa mem alignments are done using the current pair of read files
        samtools used for BAM creation and sorting

2. Sample QC and cleaning

        fastqc is run on all fastq files and sorted bam files
        qualimap multi-bamqc is performed on all sorted bam files
        multiqc combines QC metrics into an html
        write_summary.py writes and excel file with key QC metrics for all fastq files 
        directory is cleaned and unnecesary/large files are deleted or zipped.
