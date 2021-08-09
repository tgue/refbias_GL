import sys

include_ids={}

f=open(sys.argv[1])

for l in f:
	include_ids[l.rstrip()]=1

f.close()

for line in sys.stdin:
	split=line.split()
	if include_ids.has_key(split[0]):
		print line.rstrip()


