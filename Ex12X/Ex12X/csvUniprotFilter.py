import csv

def csvUniprotFilter(path, colName="uniprot"):
    """ Filter uniprot values from a csv file.
    THE COLUMN NAME MUST BE UNIPROT, IF NOT ,ADD SECOND ARGUMENT """

    csvF = open(path, "r", newline="\n")
    csvFR = csv.DictReader(csvF)
    uniprot = []
    for row in csvFR:
        uniprot.append(row[colName])
    csvF.close()
    return uniprot