#!/usr/bin/bash -l

# Description: Used by master.sh to run AMRFinderPlus to find antimicrobial resistance genes from bam file
# inputs: sorted bam file, reference used to create bam file, and name of readfile
# outputs: vcf file, consensus fasta sequence, AMR resistances file

bam_path=$1
ref_path=$2
readfile=$3
readfile="${readfile%%_*}"

#Creates consensus fasta sequence from read files to be used in amrfinder

# call  variants
bcftools mpileup -f $ref_path $bam_path --threads 8 | bcftools call -mv -Oz -o $readfile.vcf.gz --threads 8
bcftools index $readfile.vcf.gz --threads 8

# normalize indels --- might be unnecessary
#bcftools norm -f $ref_path $readfile.vcf.gz -Ob -o $readfile.norm.bcf

# filter adjacent indels within 5bp --- might be unnecessary
#bcftools filter --IndelGap 5 calls.norm.bcf -Ob -o calls.norm.flt-indels.bcf

# apply variants to create consensus sequence
cat $ref_path | bcftools consensus $readfile.vcf.gz > "${readfile}_consensus.fasta"


#bcftools call -mv -Oz -o $readfile.vcf.gz tabix $readfile.vcf.gz cat $ref_path | \
#bcftools consensus $readfile.vcf.gz > $readfile.consensus.fa

# Second option for consensus generation
#samtools mpileup -d 1000 -A -Q 0 test.bam | ivar consensus -p consensus.out -q 20 -t 0

docker run --rm -v ${PWD}:/data ncbi/amr amrfinder -n "${readfile}_consensus.fasta" --threads 8 > "${readfile}_AMRout.txt"


# ABRicate command
#docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) -it staphb/abricate abricate * > abricate_out.tsv
#docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) staphb/mlst mlst *.fa > mlst_out.tsv
