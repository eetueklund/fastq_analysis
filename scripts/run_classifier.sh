#!/bin/bash


# takes in a set of paired fastq files and extracts the top 5 kraken2 results on the species level
now=$(date +"%Y-%m-%d")
touch top_kraken_results_"${now}".txt

fastqR1=$1
fastqR2=$2

#i=0

top4=()
#ID_PATH="${R1_list[i]}" #%%-*}"

ID2=$(echo $fastqR1 | sed 's:.*/::')

# 2023GO (CRE) samples have their ID after the first dash
#if [[ $ID == "2023GO" || $ID == "GCWGS" ]]; then
ID=$(echo $ID2 | cut -f1,2 -d'_')



# Bacterial DB: ~/7T/kraken2/bacterial_db
# Minikraken DB: ~/7T/databases/minikraken2_v2_8GB_201904_UPDATE
# standard DB: ~/7T/kraken2/standard_db
# viral DB: ~/7T/kraken2/viral_db
# additional options: --minimum-hit-groups 4 --confidence 0.1

# Run kraken instead so we dont have to gunzip and gzip every fastq file
~/7T/kraken2/kraken2 --db ~/7T/kraken2/standard_db --paired --use-names --gzip-compressed --threads 8 \
--report "${ID}_k2_report.txt" $fastqR1 $fastqR2 >/dev/null 2>&1


# get kraken results
percentage_match=( $(awk '{print $1}' "${ID}_k2_report.txt") )
species=( $(awk '{print $4}' "${ID}_k2_report.txt") )
taxonomy=( $(awk '{print $6","}' "${ID}_k2_report.txt") ) #tmp_file.txt) )
taxonomy2=( $(awk '{print $7","}' "${ID}_k2_report.txt") )
taxonomy3=( $(awk '{print $8","}' "${ID}_k2_report.txt") )
IFS=$'\n' tax_array=(${taxonomy[@]})
index_lst=()


# goes through kraken2 results file and extracts all species
for j in "${!species[@]}"; do
	if [[ "${species[$j]}" == "S" || "${species[$j]}" == "U" ]]; then

		top4=(${top4[@]} "${percentage_match[j]}")

		# if third index is empty, then the species name is just two words
		# else, use include third index in the name (Neisseria sp. Marseille-Q6792)
		if [[ "${taxonomy3[j]}" == "," ]]; then
			taxonomy_lst=(${taxonomy_lst[@]} $(echo "${tax_array[j]}" | tr -d ',') "${taxonomy2[j]}" )
		else
			taxonomy_lst=(${taxonomy_lst[@]} $(echo "${tax_array[j]}" | tr -d ',') $(echo "${taxonomy2[j]}" | tr -d ',') "${taxonomy3[j]}" )
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


# can the top4 list be sorted by % and then sort the taxonomy in the same way so the values stay consistent?
# yes with python
#top4_sorted=$(python3 -c 'import sys; top4=sys.argv[1:6]; tax=sys.argv[6:11]; top4, tax=zip(*sorted(zip(top4, tax))); print(top4); print(tax);' "${top4[@]}" "${list[@]}")


# this sorts list of % matches and sorts taxonomy list in the same order
# has to convert % matches list from strings to floats to sort properly
END_IDX="${#list[@]}"
sorted=$(python3 -c 'import sys; end=int(sys.argv[1]); top=sys.argv[2:end+2]; tax=sys.argv[end+2:]; \
top_float = [float(ele) for ele in top]; top5, tax=zip(*sorted(zip(top_float, tax))); print(top5); print(tax);' \
"${END_IDX[@]}" "${top4[@]}" "${list[@]}") # input arguments to python command


# this has to be initialized for some reason
tmp=0
# tmp gets % taxonomy, tmp2 gets taxonomy ID
# this splits up top4_sorted list that contains both % and tax ID
for item in $sorted; do
	if [[ "$tmp" =~ ^[0-9]+$ ]]; then
		tmp=$item
	else
		tmp2=$item
	fi
done

# Removes () and '' from the two lists
tmp=( $(echo $tmp | sed 's/[()]//g') )
tmp2=( $(echo $tmp2 | sed 's/[()'']//g') )

# list is delimited with commas
# This sparates the values in lists using , and makes the arrays bash iterable
IFS=',' read -r -a array1 <<<"$tmp"
IFS=',' read -r -a array2 <<<"$tmp2"

# xargs removes quotes from % tax and taxonomy lists
# This appends top5 kraken results into a file
for k in {1..5}; do
	#echo "${array1[k]}"
	echo -n -e "${array1[$END_IDX-$k]}" : "${array2[$END_IDX-$k]}" | xargs >> top_kraken_results_"${now}".txt
done
# add new line to the end of line
echo $'\n' >> top_kraken_results_"${now}".txt


#~/7T/kraken2/kraken2 --db ~/7T/databases/minikraken2_v2_8GB_201904_UPDATE --paired --use-names --gzip-compressed \
#--threads 8 --minimum-hit-groups 4 --report 1613_k2_report.txt --output 1613_k2_output.txt \
#/home/mdhsequencing/7T/fastq_analysis/test/231031_LFFM/reads/1613-MD-M04394-231018_S34_L001_R1_001_filtered.fastq.gz \
#/home/mdhsequencing/7T/fastq_analysis/test/231031_LFFM/reads/1613-MD-M04394-231018_S34_L001_R2_001_filtered.fastq.gz



