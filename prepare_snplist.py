from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from random import shuffle
import sys



snpref=sys.argv[1]

records=list(SeqIO.parse(sys.argv[2],'fasta'))
seq_dic={}
for r in records:
	seq=r.seq
	newseq=seq.tomutable()
	seq_dic[r.id]=newseq


snppos={}
snpid={}

nucs=['A','C','G','T']

f=open(snpref)
for l in f:
	split=l.split()
	alleles=list(set(split[4:]))
	if '0' in alleles:
		alleles.remove('0')
	if 'N' in alleles:
		alleles.remove('N')
	if len(alleles)!=2:
		continue
	c=split[0]
	p=int(split[3])
	if seq_dic[c][p-1] in alleles:
		if seq_dic[c][p-1]==alleles[0]:
			ref=alleles[0]
			alt=alleles[1]
		elif seq_dic[c][p-1]==alleles[1]:
			ref=alleles[1]
			alt=alleles[0]

		print c,p,ref,alt

f.close()

