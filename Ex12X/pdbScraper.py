""" write a module which imports all following modules and solves the task 
modul1:
- scrapes uniprot to get fasta sequences and protein name with a list input of uniprot accesion numbers
- remove everything with is not part of sequence itself 
- create a dictonair for eacht uniprot ID with the corresponding [fasta sequence and protein name], check if there are dublicates, if yes remove them 
- save this dict in a pickle and json file for later use
modul2
- Store each entry in a file in a given directory with the uniprot ID as filename
- check if file is already in directory, if yes pass if no then create """

# chooser working directoriy
import os
path = r"C:\Users\kevin\OneDrive - ZHAW\KEVIN STUFF\ZHAW\MASTER\V5_1_Programming Algorithms and Data-Structures\Ex12X"
os.chdir(path)
#Get internet stuff
from urllib import request
from bs4 import BeautifulSoup as bs
import gzip
#tocheck
import io
from datetime import datetime
#to pickle / jason
import json
""" import pickle """ 


#load xml file from uniprot

def uniprotScraper(protList, save="n"):

    if isinstance(protList, io.IOBase):
        pro = []
        for i in protList:
            i.replace(",","")
            i.replace("\n", "")
            pro.append(i)
        protList = pro

    print(protList)
    input("Passt das so?")
        
    if isinstance(protList, list):
        pass
    else:
        print("The input is neither a list nor a open() object")

    for el in protList:     #check dublicates
        protList.count(el) > 1
        print(el + " is at least a dublicate")
        protList.remove(el) #because remove only removes 1 element

    dic = {}
    name = []
    orga = []
    seq = []

    for nr, p in enumerate(protList):
        re = request.urlopen('https://files.rcsb.org/download/' + p + '.pdb.gz').read() #load compressed uniprot xml site

        bsRe = bs(gzip.decompress(re),'html.parser') #from binary to readably and transform to beautifulsoup format
        #get protein name
        try:
            name.append(bsRe.uniprot.entry.protein.fullname.string)
        except:
            print("No Protein Name found")
            name.append("-")
            pass
        try:
            orga.append(bsRe.uniprot.entry.organism.select("[type~=common]")[0])
        except:
            print("No Protein Organism found")
            orga.append("-")
            pass
        try:
            seq.append(bsRe.uniprot.entry.sequence.string.lower())
        except:
            print("No Protein Sequence for "+str(protList[nr])+ " found")
            seq.append("-")
            pass
 

    dic["uniprotID"] = protList
    dic["proteinName"] = name
    dic["organism"] = orga
    dic["fastaSequence"] = seq

    if save == "y":
        saveJson(dic)
        
    return dic

def saveJson(element):
    #save json 
    JSON = open(os.path.join(path, datetime.now().strftime("%Y_%m_%d__%H:%M")+"_uniprotScraperDict.json"), "w")
    json.dumb(element, JSON)
    JSON.close()

