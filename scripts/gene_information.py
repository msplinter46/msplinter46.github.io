# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 10:32:29 2018

@author: asus
"""

from Bio import SeqIO


def ids():              #função que retorna a lista de ID's dos genes de interesse
    genesid=['19834097','19831637','19833099','19834533', '19833742', '19833208', '19833867', '19833377']
    return genesid

genesid=ids()


def parse_ficheiro(filename):               #parse do ficheiro obtido na página do genoma completo, accession = NC_002942.5
    records=SeqIO.parse(open(filename), "genbank")
    return records
    
def gerar_info(records):            #gera uma tabela com informações do gene e proteina e grava a tabela num ficheiro.
                                    #Grava também a sequencia de aminoacidos num outro ficheiro
    res=[]
    res.append("GeneID"+"\t"+ "Accession Number"+"\t"+"Locus Tag"+"\t"+"Gene"+"\t"+"Strand"+"\t"+
               "Protein accession (NCBI)"+"\t"+"Protein ID (UniProt)"+"\t"+"Entry Name"+"\t"+"Revision"+"\t"+"Protein Name"+
               "\t"+"Protein length"+"\t"+"Celular location"+ "\t"+"EC Number/TC Number"+"\t"+"Comments" +"\t"+
                "Features" +"\t"+ "GO Terms" +"\t"+"Cross References"+"\n" )
    seqs=open("sequencias.fasta","w")
    for r in records:
        for a in r.features:
            if a.type == "CDS" and a.qualifiers["db_xref"][0][7::] in genesid:
                gene_id=str(a.qualifiers["db_xref"][0][7::])
                accession="NC_002942.5"
                locus=str(a.qualifiers['locus_tag'][0])
                if "gene" in a.qualifiers:              
                    gene_name=str(a.qualifiers['gene'][0])
                else:
                    gene_name= "null"            
                strand=str(a.strand)
                protein_id=str(a.qualifiers['protein_id'][0])
                res.append(gene_id+"\t"+ accession+"\t"+ locus +"\t"+ gene_name +"\t"+ strand +"\t"+ protein_id +"\t"+ "\n") 
                seqs.write(">"+str(protein_id+"\n"))
                seqs.write(str(a.qualifiers['translation'][0]))
                seqs.write("\n")   
    f= open("information.txt", "w")
    for lin in res:
        f.write(lin)  
    f.close()
    seqs.close()


def teste():
    f=parse_ficheiro("sequence.gb")
    gerar_info(f)

teste()


