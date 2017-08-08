#!usr/bin/env python
# Converts module names in specified columns
# Usage: map_module_names.py map_file target_file col1 [... colN]
import sys
import csv

if len(sys.argv) < 3:
	print "USAGE: map_module_names.py map_file target_file col1 [... colN]"
	sys.exit(2)
	
argv=sys.argv


# read map file
mmap={}	# { old : new }
with open(argv[1]) as f:
	reader=csv.reader(f,delimiter="\t")
	for row in reader:
		if row[2] in mmap:
			print >> sys.stderr, "duplicate module name", row
			sys.exit()
		mmap[row[2]]=row[1]
		
# read input file
cols=[ int(x) for x in argv[3:] ]
with open(argv[2]) as f:
	reader=csv.reader(f,delimiter="\t")
	for row in reader:
		newrow=list(row)
		for c in cols:
			if row[c] in mmap:
				newrow[c]=mmap.get(newrow[c],newrow[c])
		print "\t".join(newrow)


