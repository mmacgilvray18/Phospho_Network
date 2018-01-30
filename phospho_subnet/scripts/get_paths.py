"""
Extracts paths at a given conf value from the path conf file.
Sif edges are last column in the path file.

pid\tconf\tnodes|...|...\t....\tedges|...|

Optionally, extracts only paths containing at least one of
specified genes (other args)

Optionally, don't include targets (nodes starting with "G_").
notarget as argument.

Prints sif file

"""
import sys
import csv

def main(argv):
	if len(argv) < 2:
		print "USAGE: python get_paths.py pathfile_with_confs.tab [min_conf ORF1 ... ] > yourfilename.sif"
		return 2

	targets=True
	temp=[]
	for a in argv:
		if a=="notarget":
			targets=False
		else:
			temp.append(a)
	argv=temp

	fn=argv[1]
	conf=float(argv[2])
	
	genes=[]
	if len(argv) > 3:
		genes=set(argv[3:])
	
	(pconf, pnode, pedge) = read_paths(fn)
	
	# first just get paths with conf
	above=[p for (p,c) in pconf.items() if c >= conf ]
	
	keep=[]
	if len(genes) == 0:
		keep=above
	else:
		for p in above:
			intersect = set.intersection(genes, pnode[p])
			if len(intersect) > 0:
				#print >> sys.stderr, p
				keep.append(p)
				
	
	edges=set()
	for p in keep:
		edges=set.union(edges, pedge[p])
		
	print >> sys.stderr, "paths >= %f with at least one of [%s]" % (conf, ", ".join(genes))
	for e in edges:
		sp=e.split()
		if not targets:
			if sp[2].find("G_")==0 or (sp[0].find("G_")==0 and "u" in sp[1]):
				continue
		print e.replace("(","").replace(")","")

def read_paths(fn):
	"""
	Read in path info file
	return (pconf, pnode, pedge)
	pconf as { p : conf }
	pnode and pedge in format: { p : set(ids) }
	"""
	pconf={}
	pnode={}
	pedge={}	
	
	with open(fn) as cf:
		cf=csv.reader(cf, delimiter="\t")		
		for row in cf:
			p=row[0]
			if "#" in p:
				continue
			pconf[p]=float(row[1])		
			pnode[p]=set(row[2].split("|"))
			pedge[p]=set(row[-1].split("|"))	
			
	return (pconf, pnode, pedge)
	
	
if __name__=="__main__":
	sys.exit(main(sys.argv))
