
""" write a module which imports all following modules and solves the task 
modul1:
- scrapes uniprot to get fasta sequences and protein name with a list input of uniprot accesion numbers
- remove everything with is not part of sequence itself 
- create a dictonair for eacht uniprot ID with the corresponding [fasta sequence and protein name], check if there are dublicates, if yes remove them 
- save this dict in a pickle and json file for later use
modul2
- Store each entry in a file in a given directory with the uniprot ID as filename
- check if file is already in directory, if yes pass if no then create """


from csvUniprotFilter import (csvUniprotFilter)
from uniprotScraper import (uniprotScraper, saveJson)
from pandas import DataFrame

uniprot = csvUniprotFilter(r"C:\Users\kevin\OneDrive - ZHAW\KEVIN STUFF\ZHAW\_PYTHON_R\R\Train_data\vulcano.csv")

scrapeDic = uniprotScraper(uniprot,save="y")


df = DataFrame(scrapeDic).to_csv("uniprotNamen.csv")

type(df.iloc[1,1])

from bs4 import BeautifulSoup as bs
from urllib import request
import gzip
import lxml
seq=[]
re = request.urlopen('https://www.uniprot.org/uniprot/' + "P01876" + '.xml.gz').read() #load compressed uniprot xml site
bsRe = bs(gzip.decompress(re), 'lxml') #from binary to readably and transform to beautifulsoup format

bsRe.find_all("sequence")[-1].string

bsRe.find_all("sequence")[-1].string

if bsRe.uniprot.entry.sequence.string == None:
    print("ne")

bsRe.find('div', class_='section', id='sequences')
bsRe.find('table', class_='databaseTable SEQUENCE')

################



