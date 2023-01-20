# Returns list of taxids from k2_report file created by kraken2
# Eetu Eklund
# Eetu.Eklund@Maryland.gov
# 12/28/2022

def extract():
    file = open("sendsketch.txt", "r")
    taxids = []
    id = ''
    first = True
    file.readline()
    file.readline()
    file.readline()
    line = file.readline()
    line = line.split("\t")
    id = line[8]
    file.close()
    return id

# prints and returns a comma separated taxid string
print(extract())
