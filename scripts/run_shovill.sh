#!/usr/bin/bash -l



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
