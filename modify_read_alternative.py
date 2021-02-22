#!/usr/bin/env python

import sys,re

from bisect import bisect_right

ref_dict={}
alt_dict={}

snpfile=open(sys.argv[1])
for l in snpfile:
	split=l.split()
	if ref_dict.has_key(split[0]):
		ref_dict[split[0]][int(split[1])]=split[2]
		alt_dict[split[0]][int(split[1])]=split[3]
	else:
		ref_dict[split[0]]={}
		alt_dict[split[0]]={}
		ref_dict[split[0]][int(split[1])]=split[2]
		alt_dict[split[0]][int(split[1])]=split[3]

snpfile.close()

chrs=ref_dict.keys()
chrs.sort()
key_dict={}

for c in chrs:
	keys=ref_dict[c].keys()
	keys.sort()
	key_dict[c]=keys
		

for l in sys.stdin:
	if l[0] == '@':
		print l.rstrip()
		continue
	split=l.rstrip().split('\t')
	read=split[9]
	readlen=len(read)
	ch=split[2]
	bp=int(split[3])

	md=l.rstrip().split('MD:Z:')[1].split('\t')[0]
	if md.find('^')>0:
		continue
	matches=re.split('\D+',md)
	mismatches=re.split('\d+',md)
	s=md
	ref=''
	i=0
	while len(s)!=0:
		i=len(ref)
		if re.findall('\A\D+',s):
			x=re.findall('\A\D+',s)[0]
			ref+=x
			s=re.sub('\A\D+','',s)
				
		elif re.findall('\A\d+',s):
			x=int(re.findall('\A\d+',s)[0])
			s=re.sub('\A\d+','',s)
			ref+=read[i:(i+x)]


	if len(ref)==len(read):

		start=bp
		end=bp+readlen
		keys=key_dict[ch]
		read=list(read)
		index=bisect_right(keys,start)-1
		split[0]="mod_"+split[0]

		for i in xrange(index,len(keys)):
			k=keys[i]
			if k<start:
				continue
			elif k>=start and k<end:
				snp_index=k-start		
				if read[snp_index]==ref_dict[ch][k]:
					read[snp_index]=alt_dict[ch][k]
				elif read[snp_index]==alt_dict[ch][k]:
					read[snp_index]=ref_dict[ch][k]
							
			elif k>end:
				break
		split[9]=''.join(read)
		print '\t'.join(split)
