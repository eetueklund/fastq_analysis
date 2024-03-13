# Writes QC Summary in Excel Sheet
# Eetu Eklund
# Eetu.Eklund@Maryland.gov
# 12/28/2022

import xlsxwriter
import re
import os
import sys
from datetime import datetime


now = datetime.now()
year = now.strftime("%Y")
day = now.strftime("%d")
month = now.strftime("%m")

dir = sys.argv[1]
workbook = xlsxwriter.Workbook(dir+'/summary_'+str(month)+'-'+str(day)+'-'+str(year[2:])+'.xlsx')
worksheet = workbook.add_worksheet('QC Summary')

headers = ["Sample", "Basic Statistics", "Per Base Sequence Quality", "Per Sequence Quality Scores",
           "Per Base Sequence Content", "Per Sequence GC Content", "Per Base N Content",
           "Sequence Duplication Levels", "Adapter Content", "Total Reads", "Coverage",
           "Mapping Quality", "% Aligned", "ANI", "Per Base Quality", "Taxonomy Prediction",
           "Top 2 Kraken Species", "% Match"]
row = 0
col = 0

# Writes column headers
for item in headers:
    worksheet.write(row, col, item)
    col += 1

col = 0
row += 1

# Could color cells below to green, yellow, red based on pass, warn, fail.
# cell_format.set_font_color('red')
# worksheet.write(0, 0, 'Wheelbarrow', cell_format)

# Extracts fastqc general statistics
f = open(dir+"/results/multiqc_data/multiqc_fastqc.txt", "r")
f.readline()
sample_list=[]
for line in f:
    line = line.split("\t")
    sample = line[0]
    basic = line[12]
    per_base_qual = line[13]
    per_sequence_qual = line[15]
    per_base_sequence_content = line[16]
    per_sequence_GC = line[17]
    per_base_N = line[18]
    sequence_dup = line[20]
    adapter_content = line[22]

    worksheet.write(row, col, sample)
    sample_list.append(sample)
    worksheet.write(row, col+1, basic)
    worksheet.write(row, col+2, per_base_qual)
    worksheet.write(row, col+3, per_sequence_qual)
    worksheet.write(row, col+4, per_base_sequence_content)
    worksheet.write(row, col+5, per_sequence_GC)
    worksheet.write(row, col+6, per_base_N)
    worksheet.write(row, col + 7, sequence_dup)
    worksheet.write(row, col + 8, adapter_content)
    row += 1
f.close()

# writes mean coverage, quality, and % aligned to every other row (one for each pair of fastq files)
f2 = open(dir+"/results/multiqc_data/multiqc_qualimap_bamqc_genome_results.txt", "r")
f2.readline()
row = 1
for l in f2:
    l = l.split("\t")
    sample = l[0]
    total_reads = l[2]
    map_qual = l[8]
    mean_cov = l[10]
    percentage_aligned = l[11]

    # worksheet.write(row, col, sample)
    worksheet.write(row, col + 9, int(float(total_reads)))
    worksheet.write(row, col+10, float(mean_cov))
    worksheet.write(row, col+11, float(map_qual))
    worksheet.write(row, col+12, float(percentage_aligned))

    row += 2

# Extracts ANI numbers
files = sorted(os.listdir(dir+"/results/sendsketch"))
row = 1

for file in files:
    # Extracts 1 ANI and Taxonomy name from each file and writes it to summary
    f = open(dir+"/results/sendsketch/"+file, 'r')
    f.readline()
    f.readline()
    f.readline()
    l = f.readline()
    l = l.split("\t")
    worksheet.write(row, col + 13, l[2])
    # removes ANSI escape sequences from Taxonomy name, might only work for one color. make this work for all ANSI
    taxonomy = l[11].split()
    #width = len(l[11].replace("\x1b[0m\n", ""))
    #worksheet.set_column(row, col+17, width)
    worksheet.write(row, col + 15, l[11].replace("\x1b[0m\n", ""))
    row += 2
    f.close()
f2.close()

row = 1
sum = 0.0
score = 0.0
n = 0.0
calculate = False
# Extracts and calculates average quality score
dirs = sorted(os.listdir(dir+"/results/fastqc/tmp"))
for dir2 in dirs:
    file2 = open(dir+"/results/fastqc/tmp/"+dir2+"/fastqc_data.txt", 'r')
    for line in file2:
        if ">>Per base sequence quality" in line:
            file2.readline()
            line = file2.readline()
            calculate = True
        if ">>END_MODULE" in line:
            calculate = False
        if calculate:
            score = line.split("\t")[1]
            sum += float(score)
            n += 1
    if n != 0.0:
        worksheet.write(row, col + 14, sum/n)
        row += 1
    sum = 0.0
    score = 0.0
    n = 0.0
    file2.close()

"""# extracts and writes confindr results
f3 = open(dir+"/results/confindr_out/confindr_report.csv", "r")
f3.readline()
row = 1
for l3 in f3:
    l3 = l3.split(",")
    ContamStatus = l3[3]
    PercentContam = l3[4]
    worksheet.write(row, col + 15, ContamStatus)
    worksheet.write(row, col + 16, PercentContam)
    row += 2
f3.close()
"""

# writes top kraken results into summary file
kraken_results = os.listdir(dir+"/results/kraken2/")
for result in kraken_results:
    if "top_kraken_results" in result:
        kraken_file = result


cell_format = workbook.add_format({'font_color': 'red', 'border_color': 'red'})
cell_format2 = workbook.add_format({'font_color': 'orange'})

combined_samples = ''.join(sample_list)

f4 = open(dir+"/results/kraken2/"+kraken_file, "r")
iter=0
new_sample=True
both=False
for line in f4:
    # if line is not empty (two empty lines between samples)
    if line != "" and new_sample==True:
        # lines with : contain percent match and taxonomy
        # only get the first (top) match using new_sample variable set to true.
        # new_sample is changed to false after top taxonomy is written
        if ":" in line:
            percent = line.split(":")[0]
            taxonomy = line.split(":")[1]
            both=True
        # lines that are empty and have no : are the lines with sample IDs
        else:
            ID=line.strip()
            both=False

        # write percent and taxonomy on rows matching with ID
        # iterate through the worksheet while in this for loop, or add another for loop to loop through with each sample
        if both==True and ID in combined_samples and ID !='':

            if percent != "":
                if float(percent) <= 80:
                    worksheet.write(iter+1, 16, taxonomy.lstrip(' '), cell_format2)
                    worksheet.write(iter+1, 17, percent, cell_format2)
                elif float(percent) <= 65:
                    worksheet.write(iter+1, 16, taxonomy.lstrip(' '), cell_format)
                    worksheet.write(iter+1, 17, percent, cell_format)
                else:
                    worksheet.write(iter+1, 16, taxonomy.lstrip(' '))
                    worksheet.write(iter+1, 17, percent)

            # gets second top kraken match
            line=f4.readline()
            percent = line.split(":")[0]
            taxonomy = line.split(":")[1]

            if percent != "":
                if float(percent) >= 10:
                    worksheet.write(iter+2, 16, taxonomy.lstrip(' '), cell_format)
                    worksheet.write(iter+2, 17, percent, cell_format)
                elif float(percent) >= 5:
                    worksheet.write(iter+1, 16, taxonomy.lstrip(' '), cell_format2)
                    worksheet.write(iter+1, 17, percent, cell_format2)
                else:
                    worksheet.write(iter+2, 16, taxonomy.lstrip(' '))
                    worksheet.write(iter+2, 17, percent)

            new_sample=False
            both=False
            iter+=2
            ID = ""
            if iter>=len(sample_list):
                break
    else:
        new_sample=True

f4.close()

# Create new sheet for MIDAS results
worksheet2 = workbook.add_worksheet("MIDAS")
headers2=["Sample", "Species", "count_reads", "coverage", "relative_abundance"]
col=0
for item in headers2:
    worksheet2.write(0, col, item)
    col += 1


# Extract MIDAS results
midas_results = os.listdir(dir+"/results/MIDAS/")
iter2=1
sample_iter=0
second_species=False
for folder in midas_results:
    f5 = open(dir+"/results/MIDAS/"+folder+"/species/species_profile.txt", "r")
    f5.readline()
    for line in f5:
        line_list = line.split('\t')
        if float(line_list[3])==0.0:
            second_species=False
            f5.close()
            sample_iter+=2
            iter2+=1
            break
        else:
            species=line_list[0]
            read_count=line_list[1]
            coverage=line_list[2]
            abundance=line_list[3]

            if second_species==True:
                worksheet2.write(iter2, 0, folder, cell_format2)
                worksheet2.write(iter2, 1, species, cell_format2)
                worksheet2.write(iter2, 2, int(read_count), cell_format2)
                worksheet2.write(iter2, 3, float(coverage), cell_format2)
                worksheet2.write(iter2, 4, float(abundance), cell_format2)
            else:
                second_species=True
                worksheet2.write(iter2, 0, folder)
                worksheet2.write(iter2, 1, species)
                worksheet2.write(iter2, 2, int(read_count))
                worksheet2.write(iter2, 3, float(coverage))
                worksheet2.write(iter2, 4, float(abundance))

            iter2+=1

workbook.close()
f5.close()
