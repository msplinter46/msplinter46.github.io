# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 14:20:24 2018

@author: asus
"""
from Bio.Blast import NCBIXML
from re import search
from Bio import Entrez
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import Phylo

handle=open("blast_results.xml","r")        
records = NCBIXML.parse(handle)
regex="LEGPH"        #esta regex vai ser usada para remover resultaods do blast da espécie Legionella pneumophila
e_thresh=0.05

Entrez.email = 'pg35360@alunos.uminho.pt'

l=[]
l.append("Protein ID"+ "\t"+"Start/end position"+"\t"+"Concensus"+"\t"+"Conserved Percentage"+"\n")

def cria_ficheiros_seqs():
    for re in records:
        seqs=open(re.query+".fasta","w")
        handle = Entrez.efetch(db="nucleotide", id=re.query, rettype="fasta", retmode="text")
        rec=handle.read()
        seqs.write(rec.strip("\n"))
        seqs.write("\n")
        for align in re.alignments[0:10]:
            for hsp in align.hsps:                     
                    if hsp.expect<e_thresh:      
                        if not search(regex, align.hit_id):
                            seqs.write(">"+align.hit_id+"\n")
                            seqs.write(hsp.sbjct+"\n")
                                
        seqs.close()
        mult_align(str(re.query))
        tree(str(re.query))

def mult_align(ids):
    
    l.append(ids+"\t")
    alClustal = AlignIO.read(ids+".clustal", "clustal") 
    summary_align = AlignInfo.SummaryInfo(alClustal)
    consensus=open(ids+"_consensus.txt", "w")
    consensus_al = summary_align.dumb_consensus()
    consensus.write("Consensus:  " + str(consensus_al)+"\n")
#    print("Consensus:  ", consensus_al)
    numcols = alClustal.get_alignment_length()
    numseqs = len(list(alClustal))
    conserved = []
    for c in range(numcols):
        col = alClustal[:,c]
        fc = col[0]
        if col.count(fc)==numseqs:
            conserved.append(c)      
#    print (conserved)

    perc = len(conserved) / numcols * 100
    consensus.write("Percentage: " + str(perc)+"\n")
    ini = 0
    end = 0
    best_ini = 0
    best_end = 0
    maxs = 0
    pos = 1
    while pos < len(conserved):
        if conserved[pos] - conserved[pos-1] > 1:
            ini = conserved[pos]
        end = conserved[pos]
        if end - ini > maxs: 
            maxs = end - ini
            best_ini = ini
            best_end = end
        pos = pos + 1

    consensus.write("Initial: " + str(best_ini) + " " + str(best_end)+"\n")
    l.append(str(best_ini) + " " + str(best_end)+"\t")
    consensus.write("Conserved alignment:")
    consensus.write(str(alClustal[:,best_ini:best_end+1])+"\n")
    cons = consensus_al[best_ini:best_end+1]
    consensus.write("Consensus:"  + str(cons)+"\n")
    l.append(str(cons)+"\t")
    p=round(perc,2)
    l.append(str(p)+"%"+"\n")

def tree(ids):
    tree = Phylo.read(ids+".txt", "newick")
    
    
    Phylo.draw_ascii(tree)
    
def create_table():
    file=open("Multiple_align.txt","w")
    file.write(''.join(str(e) for e in l))

def test():
    cria_ficheiros_seqs()
    create_table()
    

test()


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
