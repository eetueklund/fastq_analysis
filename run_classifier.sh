#!/bin/bash


# takes in a set of paired fastq files and extracts the top 5 kraken2 results on the species level
now=$(date +"%Y-%m-%d")
touch top_kraken_results_"${now}".txt

fastqR1=$1
fastqR2=$2

#i=0

top4=(0.0)
#ID_PATH="${R1_list[i]}" #%%-*}"

# this only works for certain specific paths. need to fix to always get the first thing after the last /

ID2=$(echo $fastqR1 | sed 's:.*/::')

#ID="${ID2%%-*}"

# 2023GO (CRE) samples have their ID after the first dash
#if [[ $ID == "2023GO" || $ID == "GCWGS" ]]; then
ID=$(echo $ID2 | cut -f1,2 -d'-')
#fi

ID=${ID%".fastq.gz"}

# Bacterial DB: ~/7T/kraken2/bacterial_db
# Minikraken DB: ~/7T/databases/minikraken2_v2_8GB_201904_UPDATE
# additional options: --minimum-hit-groups 4 --confidence 0.1

# Run kraken instead so we dont have to gunzip and gzip every fastq file
~/7T/kraken2/kraken2 --db ~/7T/kraken2/bacterial_db --paired --use-names --gzip-compressed --threads 8 \
--report "${ID}_k2_report.txt" $fastqR1 $fastqR2 >/dev/null 2>&1

# for some reason when there is a comma in the kraken2 results file it seems to mess up the splitting of taxonomy
# this removes all commas from kraken results and creates a tmp_file
sed 's/,//g'  "${ID}_k2_report.txt"  > tmp_file.txt

# extract top 4 kraken matches into a tsv
# get kraken results
percentage_match=( $(awk '{print $1}' "${ID}_k2_report.txt") )
species=( $(awk '{print $4}' "${ID}_k2_report.txt") )
taxonomy=( $(awk '{print $6","}' "${ID}_k2_report.txt") ) #tmp_file.txt) )
taxonomy2=( $(awk '{print $7","}' "${ID}_k2_report.txt") )
taxonomy3=( $(awk '{print $8","}' "${ID}_k2_report.txt") )
IFS=$'\n' tax_array=(${taxonomy[@]})
index_lst=()

# taxonomy separates each taxonomy as its own index (species level and subspecies level are 2 different indices)
# can these be set as one index, or use both at the same time somehow
# Above IFS is solution to split taxonomy based on new line

# goes through kraken2 results file and extracts the top 4-5 species
for j in "${!species[@]}"; do
	if [[ "${species[$j]}" == "S" || "${species[$j]}" == "U" ]]; then
		# call python to get min of top4 list
		min=$(python3 -c 'import sys; print(min(sys.argv[1:]))' "${top4[@]}")

		# if current % is larger than the minimum, add it to the list
		if (( $(echo ${percentage_match[j]} '>' $min | bc -l) )); then
			# if length is less than 5, append
			if (( "${#top4[@]}"<5 )); then
				top4=(${top4[@]} "${percentage_match[j]}")
				#IFS=$'\n' top4=($(sort -n <<<"${top4[*]}"))
				index_lst=(${index_lst[@]} $j)

				# if third index is empty, then the species name is just two words
				# else, use include third index in the name (Neisseria sp. Marseille-Q6792)
				if [[ "${taxonomy3[j]}" == "," ]]; then
					taxonomy_lst=(${taxonomy_lst[@]} $(echo "${tax_array[j]}" | tr -d ',') "${taxonomy2[j]}" )
				else
					taxonomy_lst=(${taxonomy_lst[@]} $(echo "${tax_array[j]}" | tr -d ',') $(echo "${taxonomy2[j]}" | tr -d ',') "${taxonomy3[j]}" )
				fi
				#taxonomy_lst=(${taxonomy_lst[@]} $(echo "${tax_array[j]}" | tr -d ',') "${taxonomy2[j]}" )

			# if length is 4, delete smallest value, append, and sort
			elif (( "${#top4[@]}"==5 )); then
				# only remove 0.0 from list once 4 are filled so the list has a minimum
				if (( $(echo $top4 '==' 0.0 |bc -l) )) ;then
					unset top4[0]
				else
					unset top4[4]
				fi
				top4=(${top4[@]} "${percentage_match[j]}")
				index_lst=(${index_lst[@]} $j)

				# if third index is empty, then the species name is just two words
                                # else, use include third index in the name (Neisseria sp. Marseille-Q6792)
				if [[ "${taxonomy3[j]}" == "," ]]; then
                                        taxonomy_lst=(${taxonomy_lst[@]} $(echo "${tax_array[j]}" | tr -d ',') "${taxonomy2[j]}" )
                                else
                                        taxonomy_lst=(${taxonomy_lst[@]} $(echo "${tax_array[j]}" | tr -d ',') $(echo "${taxonomy2[j]}" | tr -d ',') "${taxonomy3[j]}" )
                                fi
				#taxonomy_lst=(${taxonomy_lst[@]} $(echo "${tax_array[j]}" | tr -d ',') "${taxonomy2[j]}" )
			fi
		fi
	fi
done

# add ID to top of kraken results
echo $ID >> top_kraken_results_"${now}".txt

#tax_list=$(python3 -c 'import sys; tax=sys.argv[1:]; tax2[:] = [x for x in tax if x]; print(" ".join([str(elem) for elem in tax2]));' "${taxonomy_lst[@]}")
#tax_list=$(python3 -c 'import sys; tax=sys.argv[1:];  string=" ".join(tax); tax2=string.split(","); tax2[:] = [x for x in tax2 if x]; print(" ".join([str(elem) for elem in tax2]));' "${taxonomy_lst[@]}")
#tax_list=($(python3 -c 'import sys; tax=sys.argv[1:];  string=" ".join(tax); tax2=string.split(","); tax2[:] = [x for x in tax2 if x]; print(tax2);' "${taxonomy_lst[@]}" | tr -d '[],'))

# split taxonomy_lst by comma to get a list of full taxonomy names
IFS=',' read -r -a list <<< "${taxonomy_lst[@]}"
#IFS=$','; taxonomy_lst2=(${taxonomy_lst[@]})


# can the top4 list be sorted by % and then sort the taxonomy in the same way so the values stay consistent?
# yes with python
#top4_sorted=$(python3 -c 'import sys; top4=sys.argv[1:6]; tax=sys.argv[6:11]; top4, tax=(list(t) for t in zip(*sorted(zip(top4, tax)))); print(top4); print(tax);' "${top4[@]}" "${taxonomy_lst[@]}")
top4_sorted=$(python3 -c 'import sys; top4=sys.argv[1:6]; tax=sys.argv[6:11]; top4, tax=zip(*sorted(zip(top4, tax))); print(top4); print(tax);' "${top4[@]}" "${list[@]}")


tmp=0
# tmp gets % taxonomy, tmp2 gets taxonomy ID
# this splits up top4_sorted list that contains both % and tax ID
for item in $top4_sorted; do
	if [[ "$tmp" =~ ^[0-9]+$ ]]; then
		tmp=$item
	else
		tmp2=$item
	fi
done

# removes () from lists
IFS=',()' read -r -a array1 <<<$tmp
IFS=',()' read -r -a array2 <<<$tmp2

# xargs removes quotes from % tax and taxonomy lists
# This appends top5 kraken results into a file
for k in "${!array1[@]}"; do
	if (( $k!=5 )); then
		echo -n -e "${array1[5-$k]}" : "${array2[5-$k]}" | xargs >> top_kraken_results_"${now}".txt
       	fi
done

echo $'\n' >> top_kraken_results_"${now}".txt
rm tmp_file.txt
x=$((i+1))
i=$x


#~/7T/kraken2/kraken2 --db ~/7T/databases/minikraken2_v2_8GB_201904_UPDATE --paired --use-names --gzip-compressed \
#--threads 8 --minimum-hit-groups 4 --report 1613_k2_report.txt --output 1613_k2_output.txt \
#/home/mdhsequencing/7T/fastq_analysis/test/231031_LFFM/reads/1613-MD-M04394-231018_S34_L001_R1_001_filtered.fastq.gz \
#/home/mdhsequencing/7T/fastq_analysis/test/231031_LFFM/reads/1613-MD-M04394-231018_S34_L001_R2_001_filtered.fastq.gz



