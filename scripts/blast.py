# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 11:47:41 2018

@author: asus
"""

from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from re import search
from Bio import ExPASy
from Bio import SwissProt


record = SeqIO.parse(open("sequencias.fasta"),"fasta")  #abre e faz parse das sequencias proteicas obtidas anteriormente 
p = open("blast_results.xml","w")   #abre um novo ficheiro para guardar resultados
for rec in record:
    records = NCBIWWW.qblast("blastp","swissprot", rec.format("fasta"))    #faz um blast para cada proteina na swissprot
    for r in records.readlines():
        p.write(r)
p.close()


handle=open("blast_results.xml","r")        
records = NCBIXML.parse(handle)
regex="LEGP."        #esta regex vai ser usada para remover resultaods do blast da espécie Legionella pneumophila
e_thresh=0.05        #e-value maximo que irá ser aceite

    
dump = open("dump.txt","w")     #ficheiro onde serão guardados resultados



regex_query=">YP"
regex_id=".*\|(.*)\..*\|.*"
regex_func="FUNCTION:(.*)CATALYTIC"
regex_func2="FUNCTION:(.*)COFACTOR"
regex_name="(.*);"



def get_uniprot(ids):           #procura na uniprot as proteinas e retorn a função, se disponivel
    handle=  ExPASy.get_sprot_raw(ids)
    record=SwissProt.read(handle)
    com=str(record.comments)
    if search(regex_func, com):
        func=str(search(regex_func,com).group(1))
        return func
    elif search(regex_func2, com):
        func=str(search(regex_func2,com).group(1))
        return func
    else:
        return "Not available"
    

def blast_table(line):
    if search(regex_id, line):
        uni_id=str(search(regex_id, line).group(1))    
        if uni_id!="UniprotID":
            record= get_uniprot(uni_id)
            return record


def grava_tabela():
    for re in records:     #para cada record serão guardados num novo ficheiro partes relevantes dos resultados filtrados de 
                           #forma a remover resultados da espécie em estudo e e-values demasiado altos
        dump.write(">"+re.query+"\n")
        dump.write("gi|UniprotID|Entry_name"+"\t"+"Title"+"\t"+"e-value"+"\t"+"Score"+"\t"+"Align length"+"\t"+"Function"+"\n")
        for alignment in re.alignments:
            for hsp in alignment.hsps:                     
                if hsp.expect<e_thresh:      
                    if not search(regex, alignment.hit_id):
                        hit_id=alignment.hit_id.strip("gi")
                        if search("sp", hit_id):
                            hit_id1=str(search("(.*)sp(.*)", hit_id).group(1)).strip("|\s")
                            hit_id2=str(search("(.*)sp(.*)", hit_id).group(2)) 
                            function=blast_table(hit_id2).strip("', '")
                            name=alignment.hit_def.strip("RecName: Full=")
                            if search(regex_name, name):
                                name=search(regex_name, name).group(1)
                                if search(regex_name, name):
                                    name=search(regex_name, name).group(1)                         
                            dump.write(hit_id1+hit_id2+"\t"+name +"\t" 
                                       +str(hsp.expect)+"\t"+str(hsp.score)+"\t"+str(hsp.align_length)+"\t"+function+"\n")
        
                   
        dump.write('\n'*3)
    
    dump.close()    

grava_tabela()

















