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

from numpy import isin
path = r"C:\Users\Kevin.Yar\Projects\biognosys-research\Kevin_Yar\_DataPrep_Simultan\menu3"
os.chdir(path)
#Get internet stuff
from urllib import request
from bs4 import BeautifulSoup as bs
import gzip
import lxml
#tocheck
import io
from datetime import datetime
from time import sleep
#to pickle / jason
import json


import pandas as pd
import numpy as np

import re


########################################################################################################
########################################################################################################
#load xml file from uniprot

def uniprotScraper(protList, save="n", path=""):
    """ this function scrapes through uniprot and gets the most important information like, proteinname, organism and AA sequence
    Just insert protList a list with uniprot ID's and select if you want to save the JSON file with save to be either n = no oder y= yes.
    Furthermore add an output FOLDER path. If you want to save the JSON file in a panda just pandas.read_json() this thing"""

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
    elif isinstance(protList, pd.Series) or isinstance(protList, np.ndarray):
        protList = list(protList)
    else:
        print("The input is neither a list nor a open() object")

    for el in protList:     #check dublicates
        if protList.count(el) > 1:
            print(el + " is at least a dublicate")
            protList.remove(el) #because remove only removes 1 element

    dic = {}
    name = []
    gene = []
    orga = []
    seq = []
    liLen = len(protList)
    animation = ["[■□□□□□□□□□]","[■■□□□□□□□□]", "[■■■□□□□□□□]", "[■■■■□□□□□□]", "[■■■■■□□□□□]", "[■■■■■■□□□□]", "[■■■■■■■□□□]", "[■■■■■■■■□□]", "[■■■■■■■■■□]", "[■■■■■■■■■■]"]
    for nr, p in enumerate(protList):

        print(f"Process: {nr} / {len(protList)}", end='\r')
        #countdown animation

            
        #get url
        request_ = request.urlopen('https://www.uniprot.org/uniprot/' + p + '.xml.gz').read() #load compressed uniprot xml site
        bsRe = bs(gzip.decompress(request_),'lxml') #from binary to readably and transform to beautifulsoup format
        #get protein name
        try:
            name.append(str(bsRe.uniprot.entry.protein.fullname.string))
        except:
            print("No Protein Name found")
            name.append("-")
            pass

            #get gene
        try:
            geneString = str(bsRe.uniprot.entry.gene.select("name"))
            geneString = re.search(r'>(.*?)<', geneString).group(1)
            gene.append(geneString)
        except:
            gene.append("-")
            pass
        #get organism
        try:
            organismString = str(bsRe.uniprot.entry.organism.select("name"))
            orga.append(re.search(r'>(.*?)<',organismString).group(1))
        except:
            orga.append("-")
            pass

        #get sequence
        try:
            if not bsRe.uniprot.entry.sequence.string == None:
                seq.append(str(bsRe.uniprot.entry.sequence.string))
            else:
                raise AttributeError
                
        except:
            try:
                seq.append(str(bsRe.find_all("sequence")[-1].string)) #searches the last element from multiple sequence entry
            except:
                print("No Protein Sequence for  found")
                seq.append("-")
                pass
            

    dic["uniprot"] = protList
    dic["proteinName"] = name
    dic["gene"] = gene
    dic["organism"] = orga
    dic["fastaSequence"] = seq

    if save == "y":
        saveJson(dic, path=path)

    print(""" |￣￣￣￣￣￣￣￣￣￣￣|
        
            Done!!

    |＿＿＿＿＿＿＿＿＿＿＿| 
        \ ( •◡• ) / 
           \    / 
            ---
            |   | """)   
    return dic

def saveJson(element, path=""):
    #save json 
    JSON = open(os.path.join(path, datetime.now().strftime("%Y_%m_%d__%H_%M")+"_uniprotScraperDict.json"), "w")
    json.dump(element, JSON)
    JSON.close()

########################################################################################################

source = r"C:\Users\Kevin.Yar\Projects\biognosys-research\Kevin_Yar\informationGathering\uniprotList.csv"
exit = r"C:\Users\Kevin.Yar\Projects\biognosys-research\Kevin_Yar\_DataPrep_Simultan\menu3"
names = pd.read_csv(source)

#get uniprot information
uniprotScraper(names.uniprot.unique(), save="y", path=exit)

#read in the created json file
yeah = pd.read_json(r"C:\Users\Kevin.Yar\Projects\biognosys-research\Kevin_Yar\_DataPrep_Simultan\menu3\2022_03_17__17_04_uniprotScraperDict.json")

#add the molecular properties to each protein
os.chdir(r"S:\Big_No_Backup\Kevin\Uniprot-fastaScraper-by-xml\Ex12X")
import biopython
biopython.proteinInfos(yeah, "y")

#query=P178616&fields=accession,id,entry name,reviewed,protein names,genes,organism,length