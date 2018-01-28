# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 10:21:45 2018

@author: asus
"""

from Bio import Entrez
from Bio import Medline
import time



Entrez.email = "pg35360@alunos.uminho.pt"
handle = Entrez.egquery(term = "Legionella  pneumophila")
record = Entrez.read(handle)

#Contar o numero total de artigos para preencher parâmetro retmax do Entrez.esearch
for row in record["eGQueryResult"]:
    if row["DbName"]=="pubmed":
        total = row["Count"]
        
        
#aceder aos identificadores pubmed de cada artigo
handle = Entrez.esearch(db = "pubmed", term = "Legionella  pneumophila Philadelphia 1", retmax=total)#restringir à estirpe Philadelphia 1
record = Entrez.read(handle)
idlist = record["IdList"]

#download da informação correspondente aos ids obtidos
handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
records = list(Medline.parse(handle))
print(records)


#Exportação da informação obtida para um ficheiro de texto 

#cabeçalho do ficheiro
results = open("Artigos PubMed.txt", 'w')
lang={'eng':'English', 'fre':'French', 'rus':'Russian', 'ger':"German", "spa":"Spanish", "pol":"Polish", "jpn":"Japanese"}
results.write ("Legionella pneumophila subsp. pneumophila str. Philadelphia 1\n\n")
results.write('Artigos Disponiveis em PubMed a '+ time.strftime("%d/%m/%Y")+'\n')
results.write('Total de artigos disponiveis:' + str(len(records)) + '\n\n')

#preenchimento das linhas com a informação respetiva a cada artigo
i=0
for record in records:
    i+=1
    lin1 ="[" + str(i) + "]\n" + "Title: " + str(record.get("TI", "?")) 
    lin2 ="\nDate:" + str(record.get("DP", "?"))
    lin3 = "\nAuthors: " + str(record.get("AU", "?")) 
    lin4 = "\nLanguage:" + lang[(record.get("LA","?"))[0]]
    lin5 = "\nSource: " + str(record.get("SO", "?")) 
    results.write(lin1)
    results.write(lin2)
    results.write(lin3)
    results.write(lin4)
    results.write(lin5)
    results.write("\n\n")
results.close()


#Mesmo procedimento que o anterior mas com temo de pesquisa referente à estirpe e via metabolica em estudo
handle2 = Entrez.esearch(db = "pubmed", term = "Legionella  pneumophila Philadelphia 1 threonine", retmax=total)#restringir à estirpe Philadelphia 1
record2 = Entrez.read(handle2)
idlist2 = record2["IdList"]

handle2 = Entrez.efetch(db="pubmed", id=idlist2, rettype="medline", retmode="text")
records2 = list(Medline.parse(handle2))


record_results2 = open("Artigos PubMed_treonina.txt", 'w')
lang={'eng':'English', 'fre':'French', 'rus':'Russian', 'ger':"German", "spa":"Spanish", "pol":"Polish", "jpn":"Japanese"}
record_results2.write ("Legionella pneumophila subsp. pneumophila str. Philadelphia 1\n\n")
record_results2.write('Artigos Disponiveis em PubMed a '+ time.strftime("%d/%m/%Y")+'\n')
record_results2.write('Total de artigos disponiveis:' + str(len(records2)) + '\n\n')

i=0
for record in records2:
    i+=1
    lin1 ="[" + str(i) + "]\n" + "Title: " + str(record.get("TI", "?")) 
    lin2 ="\nDate:" + str(record.get("DP", "?"))
    lin3 = "\nAuthors: " + str(record.get("AU", "?")) 
    #lin4 = "\nLanguage:" + lang[(record.get("LA","?"))[0]]
    lin5 = "\nSource: " + str(record.get("SO", "?")) 
    record_results2.write(lin1)
    record_results2.write(lin2)
    record_results2.write(lin3)
    #record_results2.write(lin4)
    record_results2.write(lin5)
    record_results2.write("\n\n")
record_results2.close()



















