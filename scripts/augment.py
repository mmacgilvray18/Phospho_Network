#!/usr/bin/env python
"""
Gets BG edges between path-nodes that aren't already in the paths.
Will not duplicate undirected edges
If third argument, will also bring in edge attributes from file specified in
third argument.

Prints into tab format.

Usage: python augment.py paths.sif  bg.sif [edge_attrs] > augmented.tab


"""
import sys
import csv

if len(sys.argv) < 3:
	print "Usage: python augment.py paths.sif  bg.sif [edge_attrs] > augmented.tab"
	print "Produces a table-style network file."
	print "See header comments for more info"
	sys.exit(2)

# read path edges and nodes
nodes=[]
edges=[]

# read edge feats?
feats={}
otherfields=""
if len(sys.argv)==4:
	with open(sys.argv[3]) as f:
		reader=csv.DictReader(f, delimiter="\t")
		otherfields=[ h for h in reader.fieldnames if h != "sif_style" ]
		for row in reader:
			edge=row["sif_style"]
			esp=edge.split()
			if len(esp) != 3:
				continue
			
			# self edges - skip 'em
			if esp[0]==esp[2]:
				continue
			
			if edge not in feats:
				feats[edge]=dict( [(k,set([v])) for (k,v) in row.items() if k!="sif_style" ])
			else:
				# unite values
				for r in feats[edge]:
					feats[edge][r] = set.union(feats[edge][r], set([row[r]]))

# print incoming SIF as a table.
headers="node1\tinteraction\tnode2\tedge_status"
if len(otherfields)>0:
	headers="%s\t%s" % (headers, "\t".join(otherfields))

print headers

with open(sys.argv[1]) as f:
	for line in f:
		edges.append(line.strip())			
		edge=line.strip().split()
		nodes.append(edge[0])
		nodes.append(edge[2])
		
		# also add in reverse of undirected edges so we don't take both
		if "u" in edge[1]:
			edges.append(" ".join([edge[2], edge[1], edge[0] ]))
			
		efeats=["chosen"]
		if len(otherfields)>0 and "%s (%s) %s" % tuple(edge) in feats:
			efeats=efeats + [ "|".join(list(feats["%s (%s) %s" % tuple(edge) ][k])) for k in otherfields ]
		print "%s\t%s\t%s\t%s" % ( edge[0], edge[1], edge[2], "\t".join(efeats) ) 
		
nodes=set(nodes)
edges=set(edges)
		
# read BG and get edges between nodes that aren't already in the edge set.
with open(sys.argv[2]) as f:
	for line in f:
		edge=line.strip()
		esp=edge.split()
		if esp[0] in nodes and esp[2] in nodes and edge not in edges:
			efeats=["bg_only"]
			if len(otherfields)>0:
				efeats=efeats + [ "|".join(list(feats["%s (%s) %s" % tuple(esp) ][k])) for k in otherfields ]
			print "%s\t%s\t%s\t%s" % ( esp[0], esp[1], esp[2], "\t".join(efeats) ) 
