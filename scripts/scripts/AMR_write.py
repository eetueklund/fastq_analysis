# Writes AMR Summary in Excel Sheet
# Eetu Eklund
# Eetu.Eklund@Maryland.gov
# 08/02/2023

import os
import sys
import xlsxwriter

dir = sys.argv[1]

workbook = xlsxwriter.Workbook(dir+'/AMR_summary.xlsx')
worksheet = workbook.add_worksheet()

headers = ["Sample", "Contig ID", "Gene symbol", "Sequence name", "Element type", "Element subtype", "Class", "Subclass", "Method",
           "% Identity to Reference", "Accession of Closest Protein", "Name of Closest Protein"]
row, col = 0, 0
for item in headers:
    worksheet.write(row, col, item)
    col += 1
row, col = 1, 0

contig_ID = []
gene_symbol = []
seq_name = []
type = []
subtype = []
Class = []
subclass = []
method = []
identity = []
accession = []
name = []

files = sorted(os.listdir(dir+"/results/AMRFinder"))
for file in files:
    f = open(dir+"/results/AMRFinder/"+file, "r")
    f.readline()  # Skip header line
    for line in f:
        line = line.split("\t")

        contig_ID.append(line[1])
        gene_symbol.append(line[5])
        seq_name.append(line[6])
        type.append(line[8])
        subtype.append(line[9])
        Class.append(line[10])
        subclass.append(line[11])
        method.append(line[12])
        identity.append(line[16])
        accession.append(line[18])
        name.append(line[19])


    if len(contig_ID)!=0:
        # writes name of fastq file used for AMR in first column
        # writes key AMR results in a list per cell
        worksheet.write(row, col, file.split(".")[0])
        worksheet.write(row, col + 1, ','.join(contig_ID))
        worksheet.write(row, col + 2, ','.join(gene_symbol))
        worksheet.write(row, col + 3, ','.join(seq_name))
        worksheet.write(row, col + 4, ','.join(type))
        worksheet.write(row, col + 5, ','.join(subtype))
        worksheet.write(row, col + 6, ','.join(Class))
        worksheet.write(row, col + 7, ','.join(subclass))
        worksheet.write(row, col + 8, ','.join(method))
        worksheet.write(row, col + 9, ','.join(identity))
        worksheet.write(row, col + 10, ','.join(accession))
        worksheet.write(row, col + 11, ','.join(name))
    elif len(contig_ID) == 0:
        worksheet.write(row, col, file.split(".")[0])
        worksheet.write(row, col + 1, "No Antimicrobial Resistance Genes Identified by AMRFinderPlus")

    row += 1
    f.close()
    contig_ID = []
    gene_symbol = []
    seq_name = []
    type = []
    subtype = []
    Class = []
    subclass = []
    method = []
    identity = []
    accession = []
    name = []

workbook.close()
