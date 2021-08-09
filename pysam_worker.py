
from math import log,exp
import pysam
import sys
import os.path

err_dict={}

def init_err(maxq=101):
	for q in range(maxq):
		err_dict[q]=10**(-1*q/10)



def calc_GL(bases,quals,refbase,altbase,rl=0.5):

	if rl==0.0 or rl==1.0 or len(bases)==0:
		return([1/3,1/3,1/3])

	#GL for hom ref
	GL0=0.0
	for i in range(len(quals)):
		if bases[i]==refbase:
			GL0+=log((1-err_dict[quals[i]]))
		elif bases[i]==altbase:
			GL0+=log(err_dict[quals[i]]/3)
			
	
	#GL for hom alt
	GL2=0.0
	for i in range(len(quals)):
		if bases[i]==altbase:
			GL2+=log(( 1-err_dict[quals[i]]))
		elif bases[i]==refbase:
			GL2+=log(err_dict[quals[i]]/3)	
	
	
	#GL for het
	GL1=0.0
	for i in range(len(quals)):
		if bases[i]==refbase:
			GL1+=log(rl * (1-err_dict[quals[i]]) + (1-rl) * err_dict[quals[i]]/3)
		elif bases[i]==altbase:
			GL1+=log(rl * err_dict[quals[i]]/3  + (1-rl) * (1-err_dict[quals[i]]))
	
	

	GL0=exp(GL0)	
	GL1=exp(GL1)
	GL2=exp(GL2)

	S=sum([GL0,GL1,GL2])


	GL0/=S
	GL1/=S
	GL2/=S
	GLs=[GL0,GL1,GL2]
	return(GLs)
	




####test

init_err()

bamfile=sys.argv[1]
sitesfile=sys.argv[2]
RL_file=sys.argv[3]
refgenome=sys.argv[4]
outfile=sys.argv[5]
CORRECT=True


samfile = pysam.AlignmentFile(bamfile, "rb")

genome_fasta_open = pysam.Fastafile(refgenome)

refallele_dic={}
altallele_dic={}
RBdic={}

f=open(RL_file)

for line in f:
	split=line.split()
	refallele_dic[split[5]]=split[2]
	altallele_dic[split[5]]=split[3]
	RBdic[split[5]]=float(split[4])
f.close()



if os.path.isfile(outfile):
	#outfile exist, prep lines for appending info
	f=open(outfile)
	outlines=f.read().splitlines()
	f.close()
else:
	#outfile does not exist, prep lines
	f=open(sitesfile)
	outlines=[]
	for l in f:
		split=l.rstrip().split()
		if split[0].find('chr')>=0:
			marker="%s_%s" %(split[0],split[1])
		else:
			marker="chr%s_%s" %(split[0],split[1])
		outlines.append("%s\t%s\t%s" %(marker,split[2],split[3]))
	f.close()
	
	

f=open(sitesfile)
lines=f.readlines()
f.close()

assert len(lines)==len(outlines)

for i,l in enumerate(lines):
	split=l.split()
	c=split[0]
	pos=int(split[1])
	a1=split[2]
	a2=split[3]

	out_new='%s\t%s\t%s\n' %(1/3,1/3,1/3)

	for pileupcolumn in samfile.pileup(c, pos-1, pos,min_mapping_quality=30,truncate=True):

		refbase=genome_fasta_open.fetch(reference=c,start=pos-1,end=pos)
		bases=pileupcolumn.get_query_sequences()
		quals=pileupcolumn.get_query_qualities()
		bases=''.join(bases).upper()

		if refbase==a1:
			altbase=a2
		else:
			altbase=a1

		if len(bases)==len(quals):
			GLs=calc_GL(bases,quals,refbase,altbase)
		else:
			GLs=[1/3,1/3,1/3]

		
		out_new="%s\t%s\t%s\n" %(GLs[0],GLs[1],GLs[2])
	outlines[i]='%s\t%s' %(outlines[i],out_new)

samfile.close()

f=open(outfile,'w')
f.writelines(outlines)
f.close()

print("GL calculation done")

if CORRECT:
	print("staring correction")
	samfile = pysam.AlignmentFile(bamfile, "rb")

	genome_fasta_open = pysam.Fastafile(refgenome)
	outfile=outfile+'.corrected'

	if os.path.isfile(outfile):
		#outfile exist, prep lines for appending info
		f=open(outfile)
		outlines=f.read().splitlines()
		f.close()
	else:
		#outfile does not exist, prep lines
		f=open(sitesfile)
		outlines=[]
		for l in f:
			split=l.rstrip().split()
			if split[0].find('chr')>=0:
				marker="%s_%s" %(split[0],split[1])
			else:
				marker="chr%s_%s" %(split[0],split[1])
			outlines.append("%s\t%s\t%s" %(marker,split[2],split[3]))
		f.close()
		
		

	f=open(sitesfile)
	lines=f.readlines()
	f.close()

	assert len(lines)==len(outlines)

	for i,l in enumerate(lines):
		split=l.split()
		c=split[0]
		pos=int(split[1])
		a1=split[2]
		a2=split[3]
		if split[0].find('chr')>=0:
			marker="%s_%s" %(split[0],split[1])
		else:
			marker="chr%s_%s" %(split[0],split[1])

		out_new='%s\t%s\t%s\n' %(1/3,1/3,1/3)

		for pileupcolumn in samfile.pileup(c, pos-1, pos,min_mapping_quality=30,truncate=True):

			refbase=genome_fasta_open.fetch(reference=c,start=pos-1,end=pos)
			bases=pileupcolumn.get_query_sequences()
			quals=pileupcolumn.get_query_qualities()
			bases=''.join(bases).upper()

			if refbase==a1:
				altbase=a2
			else:
				altbase=a1
			if marker in RBdic and len(bases)==len(quals):
				GLs=calc_GL(bases,quals,refbase,altbase,rl=RBdic[marker])
			else:
				GLs=[1/3,1/3,1/3]

			out_new="%s\t%s\t%s\n" %(GLs[0],GLs[1],GLs[2])
		outlines[i]='%s\t%s' %(outlines[i],out_new)

	samfile.close()

	f=open(outfile,'w')
	f.writelines(outlines)
	f.close()
