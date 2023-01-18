# Returns list of taxids from k2_report file created by kraken2
# Eetu Eklund
# Eetu.Eklund@Maryland.gov
# 12/28/2022

def extract():
    file = open("./k2_report.txt")
    taxids = []
    id = ''
    first = True
    for line in file:
        line = line.split("\t")

        # sets taxid as first S
        if line[3] == 'S' and first is True:
            id = line[4]
            break
    file.close()

    return id


# prints and returns a comma separated taxid string
print(extract())
