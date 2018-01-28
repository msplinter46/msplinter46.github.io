# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 12:50:17 2018

@author: asus
"""

import requests


query=''
fh=open("information.txt", "r")             #abre o ficheiro com informaçoes dos genes

for line in fh.readlines():
    query+=  (line.split('\t')[0]) + " "
    
url= 'http://www.uniprot.org/mapping/'

params = {
'from':'P_ENTREZGENEID',
'to':'ACC',
'format':'tab',
'query':query}

r = requests.get(url,params)            #encontra ID's da Uniprot correspondentes aos id's dos genes

new_file= open("ncbiToUniprot.txt","w")
new_file.write(r.text)                  #grava num ficheiro os id's dos genes e proteinas

new_file.close()

    

