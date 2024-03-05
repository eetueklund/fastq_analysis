#!/usr/bin/bash -l


source ~/miniconda3-2/etc/profile.d/conda.sh

conda activate nullarbor
# Build all assemblies
R1File=''
bothfiles=false
for file in *.fastq*; do

	if [ "$bothfiles" = true ]; then
		bothfiles=false
		front="${file%%_*}"
		shovill --outdir "${front}_shovill_out" --R1 $R1File --R2 $file
		# continue skips rest of for loop so bothfiles stays false
		continue
	fi
	R1File=$file
	bothfiles=true

done
conda deactivate
# run hicap on all assemblies
conda activate hicap
for folder in *_shovill_out; do
	sample="${folder%%_*}"
	mkdir -p "${sample}_hicap_out"
	hicap --query_fp "${folder}/spades.fasta" --output_dir "${sample}_hicap_out"
done
conda deactivate
