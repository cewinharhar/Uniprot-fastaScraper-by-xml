from urllib import request
from bs4 import BeautifulSoup as bs
import gzip
import lxml
#tocheck
import io
from datetime import datetime
from time import sleep

import re
import pandas as pd
import numpy as np

#file
pathIrina = r"C:\Users\kevin\OneDrive - ZHAW\KEVIN STUFF\ZHAW\_PYTHON_R\_GITHUB\Uniprot-fastaScraper-by-xml\IrinaMonn.html"
pathFranco = r"C:\Users\kevin\OneDrive - ZHAW\KEVIN STUFF\ZHAW\_PYTHON_R\_GITHUB\Uniprot-fastaScraper-by-xml\francoCurschellas.html"
pathSimon = r"C:\Users\kevin\OneDrive - ZHAW\KEVIN STUFF\ZHAW\_PYTHON_R\_GITHUB\Uniprot-fastaScraper-by-xml\SimonSgier.html"


maLi = []
for i in [pathIrina, pathFranco, pathSimon]:
    with open(i, encoding="utf8") as re:
        bsRe = bs(re,'html.parser')   
    li = []
    className = "oajrlxb2 g5ia77u1 qu0x051f esr5mh6w e9989ue4 r7d6kgcz rq0escxv nhd2j8a9 nc684nl6 p7hjln8o kvgmc6g5 cxmmr5t8 oygrvhab hcukyx3x jb3vyjys rz4wbd8a qt6c0cv9 a8nywdso i1ao9s8h esuyzwwr f1sip0of lzcic4wl gmql0nx0 gpro0wi8"
    classNameA = "d2edcug0 hpfvmrgz qv66sw1b c1et5uql lr9zc1uh a8c37x1j keod5gw0 nxhoafnm aigsh9s9 d3f4x2em fe6kdd0r mau55g9w c8b282yb mdeji52x a5q79mjw g1cxx5fr lrazzd5p oo9gr5id"
    refs = bsRe.find_all("a", class_= className)

    for i in refs:
        li.append(i.string)

    maLi.append(li)


same = [x for x in maLi[1] if x in maLi[2] and x in maLi[0]]
#is the same as 
lu=[]
for x in maLi[1]:
    if x in maLi[2] and x in maLi[0]:
        lu.append(x)

#####################facebook mobile webseite
re =  request.urlopen("https://mobile.facebook.com/irina.monn/friends").read()  
bsRe = bs(re,'lxml') 
bsRe

className = "_52jh _5pxc _8yo0"

bsRe.find_all("div", class_= className)
bsRe.find_all(className)

li = []
for i in refs:
    li.append(i.string)
##################### facebook webseite
re =  request.urlopen("https://www.facebook.com/irina.monn/friends").read()  
bsRe = bs(re,'html.parser') 
bsRe






re.findall("\dir=\"auto\">", refs[0])
