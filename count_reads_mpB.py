import sys,subprocess
from subprocess import Popen

reffile=sys.argv[2]
snpref=sys.argv[3]
bedfile=sys.argv[4]
bamfile=sys.argv[1]

mq=30

snppos={}
snpid={}
bedfilter={}

f=open(bedfile)
for l in f:
	split=l.split()
	bedfilter[(split[0],split[2])]=1
f.close()

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
	if not (split[0],split[3]) in bedfilter:
		continue
	snppos[(split[0],split[3])]=alleles
	snpid[(split[0],split[3])]=split[1]
f.close()



cmd='samtools mpileup -B -q %s -f %s -l %s %s' %(mq,reffile,bedfile,bamfile)
process = Popen(cmd, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
stdout,stderr = process.communicate()

for l in stdout.decode().split('\n'):
	split=l.split('\t')
	if len(split)>1 and (split[0],split[1]) in snppos:
		if len(split[4]):
			alleles=snppos[(split[0],split[1])]
			pileup=split[4].upper()
			ref_allele=split[2]
			ref_count=pileup.count('.')+pileup.count(',')
			if ref_allele==alleles[0]:
				alt_allele=alleles[1]
			elif ref_allele==alleles[1]:
				alt_allele=alleles[0]
			else:
				continue
			alt_count=pileup.count(alt_allele)
			if (ref_count+alt_count)==0:
				continue
			prop=float(ref_count)/(ref_count+alt_count)
			print(ref_count, alt_count, ref_allele, alt_allele, prop, snpid[(split[0],split[1])])


