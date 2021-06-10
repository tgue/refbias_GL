
import sys

GLfile_name=sys.argv[1]
GLfile=open(GLfile_name)
RBfile=open(sys.argv[2])
ind_num=int(sys.argv[3])


refallele_dic={}
altallele_dic={}
RBdic={}
beagle_nuc={'0':'A','1':'C','2':'G','3':'T'}

for line in RBfile:
	split=line.split()
	refallele_dic[split[5]]=split[2]
	altallele_dic[split[5]]=split[3]
	RBdic[split[5]]=float(split[4])

fout=open(GLfile_name+'.corrected','w')
fout_uncorr=open(GLfile_name+'.filter','w')

l=GLfile.readline()
fout.write(l)
fout_uncorr.write(l)

for line in GLfile:
	split=line.split()
	if not RBdic.has_key(split[0]):
		fout_uncorr.write('\t'.join(split) + '\n')
		fout.write('\t'.join(split) + '\n')
		continue	
	elif not refallele_dic.has_key(split[0]) or not altallele_dic.has_key(split[0]):
		continue
	if beagle_nuc[split[1]]==refallele_dic[split[0]] and beagle_nuc[split[2]]==altallele_dic[split[0]]:
		ref=ind_num*3
		alt=ind_num*3+2
	elif beagle_nuc[split[2]]==refallele_dic[split[0]] and beagle_nuc[split[1]]==altallele_dic[split[0]]:
		ref=ind_num*3+2
		alt=ind_num*3
	else:
		continue

	r=float(RBdic[split[0]])

	s=(1-r)**2*float(split[ref]) + 2*r*(1-r)*float(split[1+3*ind_num]) + r**2*float(split[alt])

	

	if not s>0.0:
		continue	

	fout_uncorr.write('\t'.join(split) + '\n')

	split[ref]=str((1-r)**2*float(split[ref])/s)
	split[1+3*ind_num]=str(2*r*(1-r)*float(split[1+3*ind_num])/s)
	split[alt]=str(r**2*float(split[alt])/s)

	
	fout.write('\t'.join(split) + '\n')

fout.close()
fout_uncorr.close()
		


