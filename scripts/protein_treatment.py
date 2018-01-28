# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 10:54:51 2018

@author: asus
"""
from re import search
from Bio import ExPASy
from Bio import SwissProt


def ncbi_to_uniprot():       #retorna um dicionario em que as chaves são os id's dos genes e os valores são os id's das proteinas na uniprot
    w = open("ncbiTOuniprot.txt","r")
    d = {}
    for l in w.readlines():
        d[l.split("\t")[0]] = l.split("\t")[1].replace("\n","")
    del d["From"]
    w.close()
    return d

def gerar_protein_info(ids):        #procura a proteina na uniprot e gera várias informações que são depois adicionadas a uma lista
    try:
        handle = ExPASy.get_sprot_raw(ids)
        swiss_record = SwissProt.read(handle) #Cria objeto SwissProt.Record
        entry_name=str(swiss_record.entry_name)
        length=str(swiss_record.sequence_length)
        description=str(swiss_record.description).strip("RecName: Full=")
        features=str(swiss_record.features).strip("']")
        features=features.strip("['")
        keywords=str(swiss_record.keywords).strip("']")
        keywords=keywords.strip("['")
        accessions=str(swiss_record.accessions).strip("']")
        accessions=accessions.strip("['")
        comments=str(swiss_record.comments).strip("']")
        comments=comments.strip("['")
        cross_refs=str(swiss_record.cross_references).strip("']")
        cross_refs=cross_refs.strip("['")
        dataclass=str(swiss_record.data_class)      
        handle.close()
    except Exception as e:
       pass
    regex_location= "SUBCELLULAR LOCATION:\s(.*)\s{.*\."
    regex_ec= "EC=(.*){.*;"
    l=[]
    if search(regex_location, comments):            #encontra a localização da proteina, e depois o ec number de forma automatica
        location=search(regex_location, comments).group(1)
    else:
        location="Unknown"
    if search(regex_ec, description):
        ec=search(regex_ec, description).group(1)
    else:
        ec="Unfound"
    l.append(accessions+"\t"+entry_name+"\t"+dataclass+"\t"+description+"\t"+ length+"\t"+location+"\t"
             +ec+"\t"+comments +"\t"+ features +"\t"+ keywords+"\t"+cross_refs+"\t")
    protein_info=''.join(str(e) for e in l)   
    
    return protein_info
    
    
def add_to_table():         #adiciona as informações obtidas na função anterior à tabela relativa aos genes
 
    fh= open("information.txt", "r")
    final=open("final_table.txt", "w")
    lines=[]
    final.write(fh.readline())
    for line in fh:
        line=line.replace("\n","")
        lines.append(line)
    d=ncbi_to_uniprot() 
    for key in d.keys():
        for line in lines:
            if search(key, line):
                line+=gerar_protein_info(d[key])
                final.write(line+"\n")
                  
    fh.close()
    final.close()

add_to_table()
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
