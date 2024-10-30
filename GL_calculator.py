
from math import log,exp
import sys
import os.path

USE_PYSAM=0
bedfile='tmp.bed'
nucs=['A','C','G','T']
pileup_char=nucs+[',']+['.']
GL_miss=[1/3,1/3,1/3]

MIN_MQ=30
MIN_BQ=30

err_dict={}

def init_err(maxq=101):
	for q in range(maxq):
		err_dict[q]=10**(-1*q/10)



def calc_GL(bases,quals,refbase,altbase,rl=0.5):

	if len(bases)==0:
		return(GL_miss)

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
	
	

	eGL0=exp(GL0)	
	eGL1=exp(GL1)
	eGL2=exp(GL2)

	S=sum([eGL0,eGL1,eGL2])

	if S==0.0:
		max_value = max([GL0,GL1,GL2])
		max_index = [i for i, x in enumerate([GL0,GL1,GL2]) if x == max_value][0]
		GLs=[0.0,0.0,0.0]
		GLs[max_index]=1.0		
	else:
		eGL0/=S
		eGL1/=S
		eGL2/=S
		GLs=[eGL0,eGL1,eGL2]

	return(GLs)
	




####test

init_err()

bamfile=sys.argv[1]
sitesfile=sys.argv[2]
RL_file=sys.argv[3]
refgenome=sys.argv[4]
outfile=sys.argv[5]
CORRECT=True

if len(sys.argv)>6:
	cores=2 #int(sys.argv[6])/2



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



outlines_index={}


if os.path.isfile(outfile):
	#outfile exist, prep lines for appending info
	f=open(outfile)
	outlines=f.read().splitlines()
	ind_count=int(len(outlines[0].split())/3 - 1)
	outlines[0]+="\tInd%s\tInd%s\tInd%s\n" %(ind_count,ind_count,ind_count)
	f.close()
	for i,l in enumerate(outlines):
		split=l.rstrip().split()
		marker=split[0]		
		outlines_index[marker]=i
else:
	#outfile does not exist, prep lines
	f=open(sitesfile)
	sitelines=f.readlines()
	f.close()
	outlines=['']*(len(sitelines)+1)
	outlines[0]="marker\tallele1\tallele2\tInd0\tInd0\tInd0\n"
	for i,l in enumerate(sitelines):
		split=l.rstrip().split()
		if split[0].find('chr')>=0:
			marker="%s_%s" %(split[0],split[1])
		else:
			marker="chr%s_%s" %(split[0],split[1])
		outlines[i+1]="%s\t%s\t%s" %(marker,split[2],split[3])
		outlines_index[marker]=i+1


if CORRECT:
	print("staring correction")
	#samfile = pysam.AlignmentFile(bamfile, "rb")

	#genome_fasta_open = pysam.Fastafile(refgenome)
	outfile_cor=outfile+'.corrected'

	if os.path.isfile(outfile_cor):
		#outfile exist, prep lines for appending info
		f_cor=open(outfile_cor)
		outlines_cor=f_cor.read().splitlines()
		ind_count=int(len(outlines_cor[0].split())/3 - 1)
		outlines_cor[0]+="\tInd%s\tInd%s\tInd%s\n" %(ind_count,ind_count,ind_count)
		f_cor.close()
	else:
		#outfile does not exist, prep lines
		f=open(sitesfile)
		sitelines=f.readlines()
		f.close()
		outlines_cor=['']*(len(sitelines)+1)
		outlines_cor[0]="marker\tallele1\tallele2\tInd0\tInd0\tInd0\n"

		for i,l in enumerate(sitelines):
			split=l.rstrip().split()
			if split[0].find('chr')>=0:
				marker="%s_%s" %(split[0],split[1])
			else:
				marker="chr%s_%s" %(split[0],split[1])
			outlines_cor[i+1]="%s\t%s\t%s" %(marker,split[2],split[3])
		f.close()
		
		
	

f=open(sitesfile)
lines=f.readlines()
f.close()

#assert len(lines)==len(outlines)

if USE_PYSAM:
	import pysam
	samfile = pysam.AlignmentFile(bamfile, "rb")
	genome_fasta_open = pysam.Fastafile(refgenome)

	for i,l in enumerate(lines): #this loop turns out to be very slow need to find alternatives to individual pileup calls for each site
		split=l.split()
		c=split[0]
		pos=int(split[1])
		a1=split[2]
		a2=split[3]

		out_new='%s\t%s\t%s\n' %(1/3,1/3,1/3)
		out_new_cor='%s\t%s\t%s\n' %(1/3,1/3,1/3)
		if split[0].find('chr')>=0:
			marker="%s_%s" %(split[0],split[1])
		else:
			marker="chr%s_%s" %(split[0],split[1])

		for pileupcolumn in samfile.pileup(c, pos-1, pos,min_mapping_quality=MIN_MQ,truncate=True):

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
				GLs=GL_miss

			
			out_new="%s\t%s\t%s\n" %(GLs[0],GLs[1],GLs[2])
			
			if CORRECT:

				if refbase==a1:
					altbase=a2
				else:
					altbase=a1
				if marker in RBdic and len(bases)==len(quals):
					GLs=calc_GL(bases,quals,refbase,altbase,rl=RBdic[marker])
				else:
					GLs=GL_miss

				out_new_cor="%s\t%s\t%s\n" %(GLs[0],GLs[1],GLs[2])
			
				
		outlines[i+1]='%s\t%s' %(outlines[i+1],out_new) #i+1 since the first line contains ind info
		if CORRECT:
			outlines_cor[i+1]='%s\t%s' %(outlines_cor[i+1],out_new_cor) #i+1 since the first line contains ind info



	samfile.close()

else: #use subprocess to call samtools instead of using pysam for accessing the bam files -- seems to be a lot faster


	import subprocess	
	from subprocess import Popen

	f=open(bedfile,'w')

	a1_dic={}
	a2_dic={}
	found_dic={}

	for l in lines:
		split=l.split()
		c=split[0]
		pos=int(split[1])
		f.write('%s\t%s\t%s\n'%(c,pos-1,pos))
		a1_dic[(c,pos)]=split[2]
		a2_dic[(c,pos)]=split[3]
		found_dic[(c,pos)]=0
	f.close()
	


	cmd='samtools mpileup -B -q %s -Q %s -f %s -l %s %s' %(MIN_MQ,MIN_BQ,refgenome,bedfile,bamfile)
	process = Popen(cmd, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	
	for l in iter(process.stdout.readline,b''):
		split=l.decode().split('\t')
		#print(split)
		if len(split)>1 and (split[0],int(split[1])) in a1_dic:

			refbase=split[2]
			a1=a1_dic[(split[0],int(split[1]))]
			a2=a2_dic[(split[0],int(split[1]))]

			out_new='%s\t%s\t%s\n' %(1/3,1/3,1/3)
			out_new_cor='%s\t%s\t%s\n' %(1/3,1/3,1/3)
			if split[0].find('chr')>=0:
				marker="%s_%s" %(split[0],split[1])
			else:
				marker="chr%s_%s" %(split[0],split[1])

			if refbase==a1:
				altbase=a2
			else:
				altbase=a1

			

			if '-' in split[4] or '+' in split[4]: #exclude sites with indels
				skip=1

			else:			
				skip=0				
				pileup=split[4].rstrip().upper()
		
				bases=[]
				quals=[]
				count=0
				#print(split)
				skip_next=0
				for i,x in enumerate(pileup):
					if skip_next:
						skip_next=0
						continue
					elif x in nucs:
						#print(split[5][count])
						#print(ord(split[5][count]))
						#print(x)
						bases.append(x)
						quals.append(ord(split[5][count])-33)
						count=count+1
					elif x in [',','.']:
						#print(split[5][count])
						#print(ord(split[5][count]))
						#print(x)
						bases.append(refbase)
						quals.append(ord(split[5][count])-33)
						count=count+1
					elif x=='^':
						skip_next=1
					

			if len(bases)==len(quals) and not skip: 
				GLs=calc_GL(bases,quals,refbase,altbase)
			else:
				GLs=GL_miss

			
			out_new="%f\t%f\t%f\n" %(GLs[0],GLs[1],GLs[2])
			
			if CORRECT:

				if refbase==a1:
					altbase=a2
				else:
					altbase=a1
				if marker in RBdic and len(bases)==len(quals) and not skip:
					GLs=calc_GL(bases,quals,refbase,altbase,rl=RBdic[marker])
				else:
					GLs=GL_miss

				out_new_cor="%f\t%f\t%f\n" %(GLs[0],GLs[1],GLs[2])
			
				
			outlines[outlines_index[marker]]='%s\t%s' %(outlines[outlines_index[marker]],out_new)
			if CORRECT:
				outlines_cor[outlines_index[marker]]='%s\t%s' %(outlines_cor[outlines_index[marker]],out_new_cor)
			found_dic[(split[0],int(split[1]))]=1


			
	
				
	for i,l in enumerate(lines):
		split=l.split()
		c=split[0]
		pos=int(split[1])

		if not found_dic[(c,pos)]:
			outlines[i+1]='%s\t%f\t%f\t%f\n' %(outlines[i+1],1/3,1/3,1/3)
			if CORRECT:
				outlines_cor[i+1]='%s\t%f\t%f\t%f\n' %(outlines_cor[i+1],1/3,1/3,1/3)
		


			
	

f=open(outfile,'w')
f.writelines(outlines)
f.close()

if CORRECT:
	f=open(outfile_cor,'w')
	f.writelines(outlines_cor)
	f.close()

print("GL calculation done")



