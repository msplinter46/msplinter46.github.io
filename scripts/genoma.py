# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 10:32:29 2018

@author: asus
"""

from Bio import SeqIO


genesid=['19834097','19831637','19833099','19834533', '19833742', '19833208', '19833867', '19833377']


def parse_ficheiro(filename):
    records=SeqIO.parse(open(filename), "genbank")
    return records
    
def gerar_info(records):
    res=[]
    res.append("GeneID|"+"Accession Number|"+"Locus Tag|"+"Gene|"+"Strand|"+"Protein accession (NCBI)|"
               + "Protein ID (UniProt)|"+"UniProt entry name|"+"Revision|"+"Protein Name|"+"Protein length|"
               +"Celular location|"+ "EC Number/TC Number|"+"Function|"+ "Description|"+ "Features|"+ "Keywords|"+
               "Annotations|"+"Comments|" + "\n" )
    seqs=open("sequencias.txt","w")
    for r in records:
        for a in r.features:
            if a.type == "CDS" and a.qualifiers["db_xref"][0][7::] in genesid:
                gene_id=str(a.qualifiers["db_xref"][0][7::])
                accession="|NC_002942.5|"
                locus=str(a.qualifiers['locus_tag'][0])
                if "gene" in a.qualifiers:              
                    gene_name=str(a.qualifiers['gene'][0])
                else:
                    gene_name= "null"            
                strand=str(a.strand)
                protein_id=str(a.qualifiers['protein_id'][0])
                res.append(gene_id+ accession+ locus + "|" + gene_name+"|" + strand + "|" + protein_id+ "|" + "\n") 
                seqs.write(">"+str(gene_id+"\n"))
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


