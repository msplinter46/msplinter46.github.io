# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 10:54:51 2018

@author: asus
"""
from re import search
from re import finditer
from Bio import ExPASy
from Bio import SwissProt


def ncbi_to_uniprot():
    w = open("ncbiTOuniprot.txt","r")
    d = {}
    for l in w.readlines():
        d[l.split("\t")[0]] = l.split("\t")[1].replace("\n","")
    del d["From"]
    w.close()
    return d


def swissprot(geneid, d):
    res=[]
    try:
        handle = ExPASy.get_sprot_raw(geneid)
        swiss_record = SwissProt.read(handle) #Cria objeto SwissProt.Record
        entry_name="Entry name: " + str(swiss_record.entry_name)+"\n"
        length="Protein length: "+ str(swiss_record.sequence_length)+"\n"
        description="Description: "+ str(swiss_record.description).strip("RecName: Full=")+"\n"
        features="Features: "+  str(swiss_record.features)+"\n"
        keywords="Keywords: "+ str(swiss_record.keywords)+"\n"
        annotations= "Annotations: " + str(swiss_record.annotation_update)+"\n"
        res.append(entry_name + length + description + features + keywords + annotations)      
        handle.close()
        return ''.join(str(e) for e in res)
    except Exception as e:
       pass

def gerar_protein_info(filename):
    fh = open(filename+".txt", "r")
    regex_id= "AC \s* (.*);"
    regex_revision= "ID \s (.*) \s (.*);"
    regex_pname= "Full=(.*)\s{"
    regex_location= "SUBCELLULAR LOCATION:\s(.*)\s{.*\."
    regex_ec= "EC=(.*){.*;"
    regex_function= "FUNCTION:\s*(.*)((\s(.*))*)"
    l=[]
    lines=fh.readlines()
    string=''.join(str(e) for e in lines)

    if search(regex_id, string):
        ids=search(regex_id, string).group(1)
    if search(regex_revision, string):
        revision= search(regex_revision, string).group(2)
    if search(regex_pname, string):
        pname=search(regex_pname, string).group(1)
    else:
        regex_pname="Full=(.*);"
        pname=pname=search(regex_pname, string).group(1)
    if search(regex_location, string):
        location=search(regex_location, string).group(1)
    else:
        location="Unknown"
    if search(regex_ec, string):
        ec=search(regex_ec, string).group(1)
    else:
        ec="Not available"
    if search(regex_function, string):
        function= search(regex_function, string).group()
        func_aux=finditer(r"ECO:.*\s.*CATALYTIC ACTIVITY:", function)
        for x in func_aux:
            fim=x.span()[0]
            function=function[10:fim]
        function=function.replace("\n","")
        function=function.replace("CC","")
        function=function.replace("{","")
        function=function.replace("   ", "")
        function=function.rstrip()        
    else:
        function="Unknown"
    l.append(ids+"|"+ "|"+revision+"|"+pname+"|"+ "|"+location+"|"+ec+"|"+function+"|")
    protein_info=''.join(str(e) for e in l)
    return protein_info
    
    
def add_to_table():
    fh= open("information.txt", "r")
    final=open("final_table.txt", "w")
    file=open("swissprot.txt","w")
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
                file.write(swissprot(d[key],d)+"\n")  
    fh.close()
    final.close()



add_to_table()
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
