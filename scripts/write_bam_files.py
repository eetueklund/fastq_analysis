# Writes bam_files.txt file to be used by qualimap multibamqc
# Eetu Eklund
# Eetu.Eklund@Maryland.gov
# 1/5/2023

import os

files = sorted(os.listdir("./results/alignments"))
f = open("bam_files.txt", "w")

for file in files:
    if file.endswith(".sorted.bam"):
        f.write(file[:-11] + "\t" + file + "\n")

f.close()
