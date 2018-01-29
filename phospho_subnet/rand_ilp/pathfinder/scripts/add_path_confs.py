"""
Given a path file and a sigma conf file, produce a combined file that adds the
sigma conf for each path.

Usage: add_path_confs.py paths.tab sigma.tab > path_sigma.tab
"""
import sys, csv

if len(sys.argv) != 3:
	print "Usage: add_path_confs.py paths.tab sigma.tab > path_sigma.tab"
	sys.exit()
	
# read confs
confs={}
with open(sys.argv[2]) as f:
	reader=csv.reader(f, delimiter="\t")
	for row in reader:
		if "#" in row[0]:
			continue
		p=row[0].replace('"',"")
		conf=float(row[1])
		confs[p]=conf

# go through path file
with open(sys.argv[1]) as f:
	reader=csv.reader(f, delimiter="\t")
	for row in reader:
		if "#" in row[0]:
			# add column header for confs			
			print "%s\tconf\t%s" % (row[0], "\t".join(row[1:]))
			continue
		p=row[0].replace("'","")
		#if p not in confs:
			#print p
		#	continue
		conf=confs.get(p,0.0)
		newrow= "%s\t%f\t%s" % (p, conf, "\t".join(row[1:]))
		newrow=newrow.replace("'","")
		print newrow
		
		
