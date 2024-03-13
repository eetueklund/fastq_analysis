#!/bin/bash


fastqR1=''
for fastqR2 in *.fastq*; do
	if [ "$bothfiles" = true ]; then

		ID=$(echo $fastqR1 | sed 's:.*/::')

		ID=$(echo $ID | cut -f1,2 -d'_')
		echo $ID
		#ID=${ID%".fastq.gz"}

		# Bacterial DB: ~/7T/kraken2/bacterial_db
		# Minikraken DB: ~/7T/databases/minikraken2_v2_8GB_201904_UPDATE
		# standard DB: ~/7T/kraken2/standard_db
		# viral DB: ~/7T/kraken2/viral_db
		# additional options: --minimum-hit-groups 4 --confidence 0.1

		# Run kraken instead so we dont have to gunzip and gzip every fastq file
		~/7T/kraken2/kraken2 --db ~/7T/kraken2/standard_db --paired --use-names --gzip-compressed --threads 16 \
		--report "${ID}_k2_standard_report.txt" $fastqR1 $fastqR2 >/dev/null 2>&1

		#~/7T/kraken2/kraken2 --db ~/7T/kraken2/viral_db --paired --use-names --gzip-compressed --threads 8 \
                #--report "${ID}_k2_viral_report.txt" $fastqR1 $fastqR2 >/dev/null 2>&1

		bothfiles=false
                # continue skips rest of for loop so bothfiles stays false
                continue
        fi
        fastqR1=$fastqR2
        bothfiles=true
done
