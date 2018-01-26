# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 12:50:17 2018

@author: asus
"""

import requests


query=''
fh=open("information.txt", "r")

for line in fh.readlines():
    query+=  (line.split('|')[0]) + " "
    
url= 'http://www.uniprot.org/mapping/'

params = {
'from':'P_ENTREZGENEID',
'to':'ACC',
'format':'tab',
'query':query}

r = requests.get(url,params)

new_file= open("ncbiToUniprot.txt","w")
new_file.write(r.text)

new_file.close()

    

