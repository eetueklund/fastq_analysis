#!/usr/bin/bash -l


# THIS FOR PAIRED ALIGNED TO blaKPC2
#R1File=''
#bothfiles=false
#for file in *.fastq*; do
#	if [ "$bothfiles" = true ]; then
#		bothfiles=false
#		front="${file%%_*}"
#		bwa mem -v 2 -t 8 blaKPC-2.fasta $R1File $file > $front.sam
#		samtools view -u -@ 8 $front.sam | samtools sort -@ 8 -o "${front}_sorted.bam"
#		samtools index -@ 8 "${front}_sorted.bam"
#		qualimap bamqc -bam "${front}_sorted.bam" -outdir "${front}_filtered_KPC_qualimap"
		# continue skips rest of for loop so bothfiles stays false
#		continue
#	fi
#	R1File=$file
#	bothfiles=true
#done


# This for single fasta file (SRR7012993_OG.fastq.gz)
#bwa mem -v 2 -t 8 blaKPC-2.fasta SRR7012993_OG.fastq.gz > SRR7012993.sam
#samtools view -u -@ 8 SRR7012993.sam | samtools sort -@ 8 -o SRR7012993_sorted.bam
#samtools index -@ 8 SRR7012993_sorted.bam
#qualimap bamqc -bam SRR7012993_sorted.bam -outdir SRR7012993_qualimap


# THIS IS PAIRED ALIGNED TO REFERENCE STRAIN FASTA
R1File=''
bothfiles=false
for file in *.fastq*; do
       if [ "$bothfiles" = true ]; then
               bothfiles=false
               front="${file%%_*}"
               bwa mem -v 2 -t 8 GCF_003071325.1_ASM307132v1_plasmid_NZ_CP028958.1.fna $R1File $file > "${front}_plasmid.sam"
               samtools view -u -@ 8 "${front}_plasmid.sam" | samtools sort -@ 8 -o "${front}_plasmid_sorted.bam"
               samtools index -@ 8 "${front}_plasmid_sorted.bam"
               qualimap bamqc -bam "${front}_plasmid_sorted.bam" -outdir "${front}_plasmid_qualimap"
                # continue skips rest of for loop so bothfiles stays false
               continue
       fi
       R1File=$file
       bothfiles=true
done


# THIS IS PAIRED AILIGNED TO KPC2 FROM REFERENCE SEQ
#R1File=''
#bothfiles=false
#for file in *.fastq*; do
#       if [ "$bothfiles" = true ]; then
#               bothfiles=false
#               front="${file%%_*}"
#               bwa mem -v 2 -t 8 KPC2_gene_from_ref.fna $R1File $file > $front.sam
#               samtools view -u -@ 8 $front.sam | samtools sort -@ 8 -o "${front}_sorted.bam"
#               samtools index -@ 8 "${front}_sorted.bam"
#               qualimap bamqc -bam "${front}_sorted.bam" -outdir "${front}_filtered_refKPC2_qualimap"
                # continue skips rest of for loop so bothfiles stays false
#               continue
#       fi
#       R1File=$file
#       bothfiles=true
#done




#KPC2_gene_from_ref.fna
# This for single fastq file (SRR7012993_OG.fastq.gz)
#bwa mem -v 2 -t 8 KPC2_gene_from_ref.fna SRR7012993_OG.fastq.gz > SRR7012993_2.sam
#samtools view -u -@ 8 SRR7012993_2.sam | samtools sort -@ 8 -o SRR7012993_2_sorted.bam
#samtools index -@ 8 SRR7012993_2_sorted.bam
#qualimap bamqc -bam SRR7012993_2_sorted.bam -outdir SRR7012993_2_qualimap
