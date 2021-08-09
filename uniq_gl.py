import sys

f=open(sys.argv[1])
fo=open(sys.argv[1]+'.uniq','w')
ind=int(sys.argv[2])
sitesfile=open(sys.argv[3])

sites_dic={}

sites_dic['marker']=1

for l in sitesfile:
	split=l.split()
	if split[4]!='NA':
		sites_dic[split[3]]=1

last_id=''

for l in f:
	split=l.split()
	if split[0]!=last_id and split[0] in sites_dic:
		fo.write("\t".join(split[:3]+split[ind*3:(ind*3+3)])+'\n')
	last_id=split[0]

f.close()
fo.close()

