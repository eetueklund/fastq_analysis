# Allin1
## Paired fastq file all-in-one analysis script. Quality Control, BWA Alignment, Samtools, Kraken2, ANI, bcftools, multiqc, confindr...

How to use: Have all fastq files either in current directory or ./reads directory. Have all script files either in current directory or ./scripts directory. 
   
        Command: bash master.sh
        
        Command: bash scripts/master.sh

### 1. Reference Preparation and indexing

        - Done in a for loop through all fastq files.
        - Kraken2 matches current fastq file to a reference sequence and downloads it
        - bwa index, samtool faidx, and picard dict are run on the reference sequences
        - bwa mem alignments are done using the current pair of read files
        - samtools used for indexing, BAM creation, and sorting
        - bcftools to make a consensus sequence based on read files
        - FastANI to calculate Average Nucleotide Identity

### 2. Sample QC and cleaning

        - fastqc is run on all fastq files and sorted bam files
        - confindr determines levels of contamination in read files
        - qualimap multi-bamqc is performed on all sorted bam files
        - multiqc combines QC metrics into an html
        - write_summary.py writes and excel file with key QC metrics for all fastq files
        - directory is cleaned and unnecesary/large files are deleted or zipped.
