# Allin1
## Paired fastq file all-in-one analysis script. Quality Control, BWA Alignment, Samtools, sendsketch, bcftools, AMRFinder+, confindr, multiqc...

### How to use: 
All fastq files should be located in current working directory. Have all script files either in current working directory or ./scripts directory. 

First argument should be the date in YYMMDD format. This is where all files will be moved to.

Optional second argument is the string "runAMR". This tells the program to run AMRFinder+. (Run time ~doubles)

Conda Environment is automatically created and downloads all required packages if the environment does not already exist. (conda install step sometimes fails. This step should be containerized in the future)
### EXAMPLES:
   
        Command: bash master.sh 230321
        
        Command: bash scripts/master.sh 230321 runAMR

### 1. Reference Preparation and indexing

        - Done in a for loop through all fastq files.
        - fastp preprocesses and filters raw reads
        - sendsketch matches current fastq file to a reference sequence and downloads it using ncbi-genome-download
        - bwa index, samtool faidx, and picard dict are run on the reference sequences
        - bwa mem alignment is performed
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
