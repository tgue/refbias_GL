import sys

pop=sys.argv[1]

#f=open("chip.selected.MAF0.2.nomiss.het.tfam")
f=open("pseudohaploid.all.tfam")
names=[]
indices=[]
for i,l in enumerate(f):
	if pop in l:
		names.append(l.split()[1])
		indices.append(i)
f.close()

all_alleles={}

for l in open('chip.selected.MAF0.2.nomiss.tped'):
	split=l.rstrip().split()
	cur_label = "%s_%s"%(split[0],split[3])
	alleles=list(set(split[4:]))
	all_alleles[cur_label]=alleles

#outfile=open('%s.chip.selected.MAF0.2.nomiss.gl.gz' %pop,'w')
outfile=open('%s.pseudohaploid.gl' %pop,'w')


outfile.write("marker\tallele1\tallele2\t")
header = '\t'.join(['\t'.join([name]*3) for name in names])
outfile.write("%s\n"%header)

for l in open('pseudohaploid.all.tped'): #"chip.selected.MAF0.2.nomiss.het.tped"
#for l in open('chip.selected.MAF0.2.nomiss.tped'):
	
	split=l.rstrip().split()
	cur_label = "%s_%s"%(split[0],split[3])

	if not all_alleles.has_key(cur_label):
		continue

	GL=[]
	alleles=list(set(split[4:]))
	if '0' in alleles:
		alleles.remove('0')
	if len(alleles)!=2:
		continue
	if len(alleles)==2:
 
		if (alleles[0]=='A' and alleles[1]=='T') or (alleles[0]=='T' and alleles[1]=='A') or (alleles[0]=='G' and alleles[1]=='C') or (alleles[0]=='C' and alleles[1]=='G'):
			continue

		if alleles[0]==all_alleles[cur_label][0] and alleles[1]==all_alleles[cur_label][1]:
			pass
		elif alleles[1]==all_alleles[cur_label][0] and alleles[0]==all_alleles[cur_label][1]:
			pass
		else:
			continue
	
	for i in indices:
		index=i*2+4
		het=1
		missing=0
		if split[index]==split[index+1] and split[index]!='0':
			het=0
		elif split[index]==split[index+1] and split[index]=='0':
			missing=1
		if not het and not missing:
			if split[index]==alleles[0]:
				GL.append('1.0')
				GL.append('0.0')
				GL.append('0.0')
	
			else:
				GL.append('0.0')
				GL.append('0.0')
				GL.append('1.0')
		elif het and not missing:
			GL.append('0.0')
			GL.append('1.0')
			GL.append('0.0')
		else:
			GL.append('0.33333333333')
			GL.append('0.33333333333')
			GL.append('0.33333333333')			
			
	outfile.write("%s\t%s\t%s\t"%(cur_label,alleles[0],alleles[1]))
	outfile.write("%s\n"%'\t'.join(GL))			
outfile.close()
