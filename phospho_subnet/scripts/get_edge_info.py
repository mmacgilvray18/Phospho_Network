"""
Gets source and PMID info for edges.
For directed edges, only gets info for sif edges A (d) B
For undirected edges, gets info for both A (u) B and B (u) A
"""
import sys, csv

# replace MKK1MKK2 ?
repmap={'YPL140C':'MKK1MKK2', 'YOR231W':'MKK1MKK2'}
dorep=True

with open(sys.argv[1]) as f:
	reader=csv.DictReader(f, delimiter="\t")
	headers=reader.fieldnames
	print "edge\tsources\tpmids"
	#print headers
	for row in reader:
		a=row["GENEA"]
		b=row["GENEB"]
		isdir=(row["dir"]=="1")
		source=row["source"]
		pmids=row["pmid"]
		
		if dorep:
			for (r,v) in repmap.items():
				a=a.replace(r,v)
				b=b.replace(r,v)
				
		if isdir:
			print "%s (d) %s\t%s\t%s" % (a,b, source, pmids)
		else:
			print "%s (u) %s\t%s\t%s" % (a,b, source, pmids)
			print "%s (u) %s\t%s\t%s" % (b,a, source, pmids)
			
			
