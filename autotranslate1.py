import sys

import os.path
f = "hie"
d = () #blank tuple
c = ''

def read_fasta(x):
    os.path.dirname(x)
    a = '' #blank string
     
    m = open(x,'r')  #open file
    a = m.readline()  #read the first line of FASTA file
    b = '' #blank string
    global c
    global d
    c = a[1:-1]       #remove the '>' sign
       
    for line in m:   # prepare DNA sequence
        b = b + line.rstrip()   # string that contains DNA sequence
        d = (c,b)
    return b
    
#########################################################################3

from string import maketrans
def reverse_complement(x):
    k = ''
    nucleo = "ACGT"
    change = "TGCA"
    comple = maketrans(nucleo,change)

    for i in x:
        k = i + k
    return k.translate(comple)
   
##########################################################################        

             #####################################################             
def get_orf(x,y):
    d = -1
    k = () 
    for i in range(len(x)):
        if x[i:i+3] == "ATG":
            for a in range(i+3,len(x),3):
                if x[a:a+3] == 'TAA' or x[a:a+3] == 'TGA' or x[a:a+3] == 'TAG':  #if TTG is a stop codon
                    d = i%3
                    seq = x
                    break 
    
    if d == -1:
        
        for i in range(len(y)):
            #seq = "no start codon"
            if y[i:i+3] == "ATG":
                #seq = "no stop codon found"
                for a in range(i+3,len(y),3):
                    if y[a:a+3] == 'TAA' or y[a:a+3] == 'TGA' or y[a:a+3] == 'TAG':  #if TTG is a stop codon
                        d = i%3
                        seq = y
                        break            
    if d == -1:
        seq = "no ORF found in this sequence"
    
    k = (d, seq)
    return k[0], k[1]
                    
##########################################################################
def get_gene_by_ORF(x,y):  #x = sequence; y = ORF on that sequence
    for i in range(y,len(x),3):
        if x[i:i+3] == "ATG":
            for a in range(i+3,len(x),3):
                if x[a:a+3] == 'TTG' or x[a:a+3] == 'TAA' or x[a:a+3] == 'TGA' or x[a:a+3] == 'TAG':
                    return x[i:a+3]
                             
##########################################################################
def translate(a):
    aa = ""
    codontable = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        }    
    for n in range(0,len(a),3):
        aa = aa + codontable[a[n:n+3]]
    return aa
        
##########################################################################
def get_fasta(a,b):
    d = b[0]
    for i in range(1,len(b)):
        width = 60
        if (i+1)%width == 0:
            d = d + b[i] + '\n'
        else:
            d = d + b[i]
    
    c = ">" + a + "\n" + d + '\n'       
    return c
    

##########################################################################

DNA_seq = read_fasta(sys.argv[1])  # 1st processing
#print DNA_seq

DNA_revcomp = reverse_complement(DNA_seq)
#print DNA_revcomp

ORF = get_orf(DNA_seq,DNA_revcomp)

if int(ORF[0]) == -1:
    print "SEQUENCE ERROR: No ORF contained in sequence!!!!!"
else:
    GENE = get_gene_by_ORF(ORF[1],int(ORF[0]))  #ORF[0] is the frame number ORF[1] is the sequence either positive or rev_comp

PROT_seq = translate(GENE)

PROT_fasta = get_fasta('protein_seq001',PROT_seq)

##########################################################################

newfasta = open(sys.argv[2],'w+')
newfasta.writelines(PROT_fasta)
newfasta.close()
