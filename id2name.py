import os,sys
sys.path.append('.')
import readBed

filename = sys.argv[1]
refname = sys.argv[2]

reader = readBed.BEDReader(filename)
convert = {}
f = open(refname)
for r in f:
	if r.startswith('#'):
		continue
	tokens = r.strip().split('\t')
	convert[tokens[0]] = tokens[1].replace(' ','_')

f.close()

out = readBed.BEDWriter(sys.argv[3])
for g in reader:
	g['name'] = convert[g['name']]
	out.writerow(g)

