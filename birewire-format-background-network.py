#!/home/mplace/anaconda3/bin/python
"""
Program: birewire-format-background-network.py

Purpose: Convert a background network SIF file which looks like:

#nodeA  interaction     nodeB   is_directed
YCR066W ubiquitination  YBR088C 1
YJL030W prenylation_both        YOR370C 0
YJL030W prenylation_both        YPR176C 0
YOR370C prenylation_both        YPR176C 0
YBL016W kinase_substrate        YLR452C 1

to a bipartite graph table.

Input  : network file

Output: Bipartite Graph Table for input to BiRewire R package.      

Created Dec 2017

@author: mplace
"""
import argparse
import os
import re
import sys
from collections import defaultdict 


class rewire( object ):
    """
    Handle the procssing and conversion of the SIF file to required BiRewire input format.
    Also enable the conversion of BiRewire output i.e. the rewired network to SIF format.    
    """
    
    def __init__(self, inFile ):
        """
        inFile = input vcf  file
        geneList is a unique list of gene names
        targetList is a unique list of gene names which are connected to those in geneList
        table is a dict of dicts where 
            key = gene name from geneList, 
            value = dict of gene names from targetlist value = 0 or 1, 
            1 means there is a connection to the gene in the dict above. 
        """
        self.inFile = inFile  
        self.geneList = []
        self.targetList = []
        self.table  = defaultdict(dict)
    
    def makeList(self):
        """
        Open and parse sif file for gene information to create a unique set of
        gene names.  Used to create dictionary.
        """        
        with open( self.inFile, 'r' ) as file:
            for line in file:
                if line.startswith('#'):
                    continue
                row = line.split()
                self.geneList.append(row[0])
                self.targetList.append(row[2])
        
        self.geneList   = sorted(list(set(self.geneList)))    # unique list of genes
        self.targetList = sorted(list(set(self.targetList)))  # unique list of targets
               
    def createDict(self):
        """
        Create the dictionary to hold the bipartite graph table.
        """
        for g in self.geneList:        # rows in final table
            for i in self.targetList:  # columns in final table
                self.table[g][i] = 0
                
    def fillGraph(self):
        """
        Open background network file (again) and set interactions to 1 in self.table
        If no interaction found it will be zero.  Only zero and one are valid values.
        Check output, only values of 0,1 are valid, i.e. 2 should not be in the final table
        """
        with open(self.inFile, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                row = line.split()
                self.table[row[0]][row[2]] += 1
                
    def writeSifNetwork(self):
        """
        Read in BiRewire output file and convert to SIF network file.        
        """
        with open(self.inFile, 'r') as file, open('temp-rewired.tmp', 'w') as out:
            targets = file.readline().split()
            for line in file:
                row = line.split()
                gene = row.pop(0)
                indices = [ i for i, x in enumerate(row) if x == '1']
                for ind in indices:
                    out.write("%s %s\n" %(gene, targets[ind]))
                    
    def mergeNtwkInfo(self, ntwk):
        """
        Open temp-rewired.tmp, and original Background Network SIF file, add interaction and direction
        from original BG Nework file, write results.
        """
        outName = re.sub( '.tab', '.csv', self.inFile)
        with open(outName,'w') as out:   
            out.write('#nodeA interaction nodeB is_directed\n')
            with open(ntwk,'r') as file, open('temp-rewired.tmp', 'r') as tmp:
                file.readline()
                for n,t in zip(file,tmp):
                    n = n.rstrip()
                    infoN = n.split()
                    t = t.rstrip()
                    infoT = t.split()
                    out.write('%s %s %s %s\n' %(infoT[0], infoN[1], infoT[1], infoN[3]) )
        
                
    def graphToFile(self):
        """
        Write bipartite graph to text file
        """
        header = '\t' + '\t'.join(self.targetList)
        
        with open('bipartiteGraph.tab', 'w') as out:
            out.write(header + '\n')
            for row in self.geneList:
                edges = [row]
                for name in self.targetList:
                    edges.append(str(self.table[row][name]))
                out.write('\t'.join(edges) + '\n')     

def main():
    cmdparser = argparse.ArgumentParser(description="Create a Bipartie table for input to BiRewire R package.",
                                        usage='%(prog)s -f file.sif -t <b|n>' , prog='birewire-format-background-network.py'  )                                  
    cmdparser.add_argument('-f', '--File', action='store', dest='FILE', help='INPUT file', metavar='')
    cmdparser.add_argument('-i', '--info', action='store_true', dest='INFO', help='Print detailed description of program.')
    cmdparser.add_argument('-t', '--type', action='store', dest='TYPE', 
                           help='Type of file to generate, BiRewire input(b), or Background network from BiRewire(n)')
    cmdparser.add_argument('-n', '--network', action='store', dest='NETWORK', help='Original SIF Network file.' )
    cmdResults = vars(cmdparser.parse_args())
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)
    
    if cmdResults['INFO']:
        print("\n  birewire-format-background-network.py -f file.sif ")
        print("\n  Purpose: Create a bipartite graph table from sif file.")
        print("\t\tOR  Create a SIF file from BiRewire output file.")
        print("\n  Input  : Standard SIF file for a protein network")
        print("\t\tOR the output table from BiRewire.")
        print("\n  Output : Tab delimited text file.")
        print("\n  For help see Mike Place\n")
        sys.exit(1)
    
    if cmdResults['FILE']:
        inFile   = cmdResults['FILE']
    
    if not os.path.exists(inFile):
        print("\n\t-f input file does not exist.\n")
        cmdparser.print_help()
        sys.exit(1)
    
    if cmdResults['TYPE'] == 'b':
        print('Convert SIF network file to BiRewire input file.')
        dat = rewire(inFile)
        dat.makeList()
        dat.createDict()
        dat.fillGraph()
        dat.graphToFile()
    elif cmdResults['TYPE'] == 'n':
        if cmdResults['NETWORK']:
            ntwkFile = cmdResults['NETWORK']
        else:
            ntwkFile = 'phospho_v4_bgnet_siflike_withdirections_fix.cvs'
        sif = rewire(inFile)
        sif.writeSifNetwork()
        sif.mergeNtwkInfo(ntwkFile)
        
        
     
if __name__ == "__main__":
    main()