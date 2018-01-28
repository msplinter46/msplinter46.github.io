# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 18:13:45 2018

@author: asus
"""
from re import search
from protein_treatment import ncbi_to_uniprot
from Bio import ExPASy
from Bio import SwissProt


def cria_table():
    final_blast=open("Blast_table.txt","w")
    ids=ncbi_to_uniprot() 
    final_blast.write("Protein ID (NCBI)"+"\t"+"Protein ID (UniProt)"+"\t"+"Protein name"+"\t"+"Molecular Function"+"\t"+"Biological Process"+"\t"+"CDD result"+"\n")
    
    for i in ids.values():       
        handle = ExPASy.get_sprot_raw(i)
        swiss_record = SwissProt.read(handle)
        pname=str(swiss_record.description).strip("RecName: Full=")  
        name=find_name(pname)
        ncbi_id=str(uniprot_to_ncbi(i).keys())
        ncbi_id=search("(YP.{9})",ncbi_id).group(1)
        final_blast.write(ncbi_id+"\t"+i+"\t"+name +"\n")
        

def find_name(description):
    regex="{ECO:(.*)}"
    regex2="{ECO:(.*)};"
    regex_ec="EC=(.*)"
    regex_short="Short=(.*)"
    l=description.split()
    for i in l:
        if search(regex_ec, i):
            l.remove(i)
    for j in l:
        if search(regex, j):
            l.remove(j)   
    for j in l:
        if search(regex2, j):
            l.remove(j)
    for k in l:
        if search(regex_short, k):
            l.remove(k)
    string=' '.join(str(e) for e in l)
    string=string.strip("SubName: Full=")
    string=string.strip("Short"+"(.*)")
    if string[len(string)-1]=="s":
        string=string+"e"
    return string
    


def uniprot_to_ncbi(unip):
    handle = ExPASy.get_sprot_raw(unip)
    swiss_record = SwissProt.read(handle)
    sr=str(swiss_record.cross_references)
    regex="(YP.{9})"
    d={}
    if search(regex,sr):
        d[(str(search(regex,sr).group(1)))]=unip
    return d
        
        
        
        
        
    


cria_table()