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
path = r"S:\Big_No_Backup\Kevin\biognosys-research\Kevin_Yar\Ex12X"
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
    orga = []
    seq = []
    liLen = len(protList)
    animation = ["[■□□□□□□□□□]","[■■□□□□□□□□]", "[■■■□□□□□□□]", "[■■■■□□□□□□]", "[■■■■■□□□□□]", "[■■■■■■□□□□]", "[■■■■■■■□□□]", "[■■■■■■■■□□]", "[■■■■■■■■■□]", "[■■■■■■■■■■]"]
    for nr, p in enumerate(protList):
        #countdown animation
        if nr < liLen * 0.1:
            print(animation[0], end='\r') 
        if  liLen * 0.2 < nr < liLen * 0.3:
            print(animation[1], end='\r')
        if liLen * 0.3 < nr < liLen * 0.4:
            print(animation[2], end='\r')
        if liLen * 0.4 < nr < liLen * 0.5:
            print(animation[3], end='\r')
        if liLen * 0.5 < nr < liLen * 0.6:
            print(animation[4], end='\r')
        if liLen * 0.6 < nr < liLen * 0.7:
            print(animation[5], end='\r')
        if liLen * 0.7 < nr < liLen * 0.8:
            print(animation[6], end='\r')
        if liLen * 0.8 < nr < liLen * 0.9:
            print(animation[7], end='\r')
        if liLen * 0.9 < nr < liLen:
            print(animation[8], end='\r')
        if nr == liLen:
            print(animation[9], end='\r')
            sleep(1)
            
        #get url
        re = request.urlopen('https://www.uniprot.org/uniprot/' + p + '.xml.gz').read() #load compressed uniprot xml site
        bsRe = bs(gzip.decompress(re),'lxml') #from binary to readably and transform to beautifulsoup format
        #get protein name
        try:
            name.append(str(bsRe.uniprot.entry.protein.fullname.string))
        except:
            print("No Protein Name found")
            name.append("-")
            pass
        #get organism
        try:
            orga.append(str(bsRe.uniprot.entry.organism.select("[type~=common]")[0]))
        except:
            print("No Protein Organism found")
            orga.append("-")
            pass


        #get sequence
        try:
            if not bsRe.uniprot.entry.sequence.string == None:
                seq.append(str(bsRe.uniprot.entry.sequence.string))
            else:
                raise AttributeError
        except AttributeError:
            seq.append(str(bsRe.find_all("sequence")[-1].string)) #searches the last element from multiple sequence entry
        except:
            print("No Protein Sequence for "+str(protList[nr])+ " found")
            seq.append("-")
            pass
 

    dic["uniprot"] = protList
    dic["proteinName"] = name
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
########################################################################################################
source = r"C:\Users\Kevin.Yar\Projects\biognosys-research\Kevin_Yar\_DataPrep_Simultan\menu3\missingProteinNames.csv"
exit = r"C:\Users\Kevin.Yar\Projects\biognosys-research\Kevin_Yar\_DataPrep_Simultan\menu3"
names = pd.read_csv(source)

uniprotScraper(names.uniprot.unique(), save="y", path=exit)



yeah = pd.read_json(r"C:\Users\Kevin.Yar\Projects\biognosys-research\Kevin_Yar\_DataPrep_Simultan\menu3\2022_02_16__15_07_uniprotScraperDict.json")
yeah.to_csv(r"uniprotNamenExtended.csv")

p = "Q9Y547"