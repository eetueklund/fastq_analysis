# Writes QC Summary in Excel Sheet
# Eetu Eklund
# Eetu.Eklund@Maryland.gov
# 12/28/2022

import xlsxwriter
import re
import os
import sys

dir = sys.argv[1]
workbook = xlsxwriter.Workbook(dir+'/summary.xlsx')
worksheet = workbook.add_worksheet()

headers = ["Sample", "Basic Statistics", "Per Base Sequence Quality", "Per Sequence Quality Scores", "Per Base Sequence Content",
          "Per Sequence GC Content", "Per Base N Content", "Sequence Duplication Levels", "Adapter Content", "Coverage",
           "Mean Mapping Quality", "% Aligned", "ANI", "Per Base Quality", "Contamination", "% Contamination", "Taxonomy"]
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
for line in f:
    line = line.split("\t")
    sample = line[0]
    basic = line[10]
    per_base_qual = line[11]
    per_sequence_qual = line[13]
    per_base_sequence_content = line[14]
    per_sequence_GC = line[15]
    per_base_N = line[16]
    sequence_dup = line[18]
    adapter_content = line[20]

    worksheet.write(row, col, sample)
    worksheet.write(row, col+1, basic)
    worksheet.write(row, col+2, per_base_qual)
    worksheet.write(row, col+3, per_sequence_qual)
    worksheet.write(row, col+4, per_base_sequence_content)
    worksheet.write(row, col+5, per_sequence_GC)
    worksheet.write(row, col+6, per_base_N)
    worksheet.write(row, col + 7, sequence_dup)
    worksheet.write(row, col + 8, adapter_content)
    # worksheet.write(row, col+7, )  coverage
    # worksheet.write(row, col+8, )  Mean Quality
    row += 1
f.close()

# writes mean coverage, quality, and % aligned to every other row (one for each pair of fastq files)
f2 = open(dir+"/results/multiqc_data/multiqc_qualimap_bamqc_genome_results.txt", "r")
f2.readline()
row = 1
for l in f2:
    l = l.split("\t")
    sample = l[0]
    map_qual = l[8]
    mean_cov = l[10]
    percentage_aligned = l[11]

    # worksheet.write(row, col, sample)
    worksheet.write(row, col+9, float(mean_cov))
    worksheet.write(row, col+10, float(map_qual))
    worksheet.write(row, col+11, float(percentage_aligned))

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
    worksheet.write(row, col + 12, l[2])
    # removes ANSI escape sequences from Taxonomy name, might only work for one color. make this work for all ANSI
    taxonomy = l[11].split()
    worksheet.write(row, col + 16, l[11].replace("\x1b[0m\n", ""))
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
        worksheet.write(row, col + 13, sum/n)
        row += 1
    sum = 0.0
    score = 0.0
    n = 0.0
    file2.close()

f3 = open(dir+"/results/confindr_out/confindr_report.csv", "r")
f3.readline()
row = 1
for l3 in f3:
    l3 = l3.split(",")
    ContamStatus = l3[3]
    PercentContam = l3[4]
    worksheet.write(row, col + 14, ContamStatus)
    worksheet.write(row, col + 15, PercentContam)
    row += 2
f3.close()
workbook.close()
