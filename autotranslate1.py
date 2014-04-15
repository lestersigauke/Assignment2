import module
import sys

DNA_seq = module.read_fasta(sys.argv[1])  # 1st processing
#print DNA_seq

DNA_revcomp = module.reverse_complement(DNA_seq)
#print DNA_revcomp

ORF = module.get_orf(DNA_seq,DNA_revcomp)

if int(ORF[0]) == -1:
    print "SEQUENCE ERROR: No ORF contained in sequence!!!!!"
else:
    GENE = module.get_gene_by_ORF(ORF[1],int(ORF[0]))  #ORF[0] is the frame number ORF[1] is the sequence either positive or rev_comp

PROT_seq = module.translate(GENE)

PROT_fasta = module.get_fasta('protein_seq001',PROT_seq)

##########################################################################

newfasta = open(sys.argv[2],'w+')
newfasta.writelines(PROT_fasta)
newfasta.close()
