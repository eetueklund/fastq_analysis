# fastq_analysis
## Paired fastq file analysis script. Quality Control, BWA Alignment, Samtools, bcftools, Kraken2, MIDAS, AMRFinder+, multiqc...

### How to use: 
All fastq files should be located in current working directory. Have all script files either in current working directory or ./scripts directory. 

First argument is the output directory name to be created. This is where all files will be moved to.

Optional second argument is the string "runAMR". This tells the program to run AMRFinder+ and also creates a consensus sequence. (Run time ~doubles)

Conda Environment is automatically activated. If it does not already exist, a fastq_analysis conda environment is created and downloads all required packages 
                                                                        (conda install step sometimes gets stuck. This step should be containerized in the future)

### EXAMPLES:
   
        Command: bash master.sh 230321
        
        Command: bash scripts/master.sh 230321_CRE runAMR

### 1. FASTQ Analysis Steps

        - fastp preprocesses and filters raw reads
        - top kraken results are extracted using run_classifier.sh as a taxonomy and contamination check
        - MIDAS is used as a double check for contamination
        - sendsketch from BBTools matches current fastq file to a reference sequence and downloads it using ncbi-genome-download
        - bwa index, samtool faidx, and picard dict are run on the reference sequences
        - bwa mem alignment is performed
        - samtools used for indexing, BAM creation, and sorting
        - optionally AMRFinder+ finds Antimicrobial resistance genes
            - This also creates a consensus fasta sequence from the .bam alignment file using bcftools
        - fastqc is run on all fastq files and sorted bam files
        - qualimap multi-bamqc is performed on all sorted bam files
        - multiqc combines QC metrics into an html
        - write_summary.py writes a summary excel file with key QC metrics and other results for all fastq files
        - directory is cleaned and unnecesary/large files are deleted or zipped
        - all results from the above tools will be in the results directory
