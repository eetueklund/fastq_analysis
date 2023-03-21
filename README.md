# Allin1
## Paired fastq file all-in-one analysis script. Quality Control, BWA Alignment, Samtools, sendsketch, bcftools, AMRFinder+, confindr, multiqc...

How to use: All fastq files should be located in current working directory. Have all script files either in current working directory or ./scripts directory. 
   
        Command: bash master.sh
        
        Command: bash scripts/master.sh

### 1. Reference Preparation and indexing

        - Done in a for loop through all fastq files.
        - sendsketch matches current fastq file to a reference sequence and downloads it using ncbi-genome-download
        - bwa index, samtool faidx, and picard dict are run on the reference sequences
        - bwa mem alignments are done using the current pair of read files
        - samtools used for indexing, BAM creation, and sorting
        - bcftools to make a consensus sequence based on read files
        - optionally AMRFinder+ finds Antimicrobial resistance genes
            - This also creates a consensus fasta sequence from the .bam alignment file

### 2. Sample QC and cleaning

        - fastqc is run on all fastq files and sorted bam files
        - confindr determines levels of contamination in read files
        - qualimap multi-bamqc is performed on all sorted bam files
        - multiqc combines QC metrics into an html
        - write_summary.py writes and excel file with key QC metrics for all fastq files
        - directory is cleaned and unnecesary/large files are deleted or zipped.
