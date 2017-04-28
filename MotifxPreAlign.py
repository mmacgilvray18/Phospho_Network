#!/home/GLBRCORG/mplace/anaconda3/bin/python
"""
Program: MotifxPreAlign.py

Purpose: Use a list of phosphorylated peptide information to create a pre-aligned
          peptide file for the motifx website. 

          Program is called from Motifx.py, but can also be run standalone from
          the command line, output from command line is a text file: motifx_pre_align_results.txt
          
          If called from another python program, the data will be in a dictionary
          called pepInfo.
          
Input: A peptide information file and a orf fasta file (optional)

    peptide info file:
        
        YBR162W-A_T5,Induced,AVQtPR,AVQT*PR
        YBL005W-B_S960,Induced,AVsPTDSTPPSTHTEDSKR,AVS*PTDSTPPSTHTEDSKR
        YBL005W-B_T962,Induced,AVSPtDSTPPSTHTEDSKR,AVSPT*DSTPPSTHTEDSKR

    orf fasta file:  
        The default is SGD R64-2-1, but user may optionally use their own.
        
    width: default with to center the phosphorylated peptide on is 13, user
        may optionally change this to any odd number.
        
output: A list of peptide sequences aligned to a specific orf, all the same
         length.  
         * Only unique peptides will be used, any duplicates will be 
         thrown out.  
         
         * Any peptide which cannot be extended because it is too close to the
         end, will be discarded.  For instance if the width is 13 but the 
         central residue is at the end of the gene and thus the sequence cannot
         be extended it will be thrown out.
         
         Removed peptides will be listed in Peptides-removed.txt
         
Dependencies: Python 3
               BioPython      
               
Output: gene name with amino acid sequence and motif
author: Mike Place
Date:   12/09/2016
version: 1.0
"""
from Bio import SeqIO
from collections import defaultdict 
import argparse	
import sys

class MotifxPreAlign(object):
    """ An exact peptide sequence matching object.
    
    Attributes:
        genes   - Gene set, taken from the peptide info file, used to limit orf search space by excluding genes not in set
        orfRec  - Fasta file to search for exact match, i.e. the orf fasta file
        pepInfo - dictionary, key is the gene name 
                  example: 
                   'YGL076C_T8_S11', {'group': 'Induced', 
                   'motif': '......S.SP...', 
                   'origPep' : 'AAEKILtPEsQLKK',
                   'phosPep' : 'AAEKILT*PES*QLKK', 
                   'extended': [ Seq('AAEKILTPESQLK',SingleLetterAlphabet())]}  # extended is a list of Seq objects
                   
                   where 'extended' is a list of peptides with the centered residue
                   having 6 peptides on either side, (default peptide length of 13)
        window  - specified width of peptide string, defaults to 13
        """

    def __init__(self, fasta, pep, window):
        """ Contruct exactFastaMatch object """
        self.window        = window                       # set width of peptide string
        self.wing          = int((window -1)/2)           # number of amino acids on either side of phospho site
        self.genes         = set()                        # gene name set, for filtering orf fasta
        self.pepInfo       = self.getPepInfo(pep)         # get peptide information
        self.orfRec        = self.filterFasta(fasta)      # make a reduced list of seqIO objects (fasta sequences)
        self.genesToRemove = []                           # list of genes to remove from final results, i.e. no match or its been excluded
        self.duplicatePeps = set()                        # used to find duplicate peptides                

    def getPepInfo(self, pep):
        """ parse and return a dict of peptide information & sequence """
        peplst = defaultdict(dict)
        
        with open(pep, 'r') as file:
            for ln in file:
                if ln.startswith('Ppep'):
                    continue
                ln = ln.rstrip()
                row = ln.split(',')
                peplst[row[0]]['group']   = row[1]
                peplst[row[0]]['origPep'] = row[2].upper()  # capitalize peptide
                peplst[row[0]]['phosPep'] = row[3]
                peplst[row[0]]['extended']= []
                peplst[row[0]]['motif'] = ''                # rhis is blank by intention, expected to be filled by calling program.
                
                gene = row[0].split('_')
                self.genes.add(gene[0])             #Get gene name 
                
        return peplst
       
    def filterFasta(self,fasta):
        """ select only those fasta sequences that match the names in the gene list """
        # open and create a list of all input fasta sequences
        record      = list(SeqIO.parse( fasta, 'fasta'))
        filteredRec = []
        # loop through all sequence records and match name to gene name
        # if names match keep record
        for rec in record:
            if rec.name in self.genes:
                # if the * appears at the end of the peptide remove it
                if rec[-1] == '*':
                    rec = rec[0:-1]
                filteredRec.append(rec)
        
        return filteredRec 
        
    def matchSeq(self):
        """ Attempt to match the sequences to the reduced gene list"""
        # loop through short sequence to match to fasta records
        for gene, info in self.pepInfo.items():
            # locate the * index and subtracts 1 to get the phosphorylated site, will handle multiple sites w/in a peptide
            idx = [ i-1 for i, ltr in enumerate(info['phosPep']) if ltr == '*']
            count = 0
            for i in idx:
                # match the peptide to a gene
                for subject in self.orfRec:            
                    start = subject.seq.find(info['origPep'])       # the index of the start position
                    if start >= 0:                                  # if there is an exact match
                        orfLength = len(subject.seq)                # get the length of gene
                        end       = start + i + self.wing    # add start index to the full peptide length
                        # exclude peptides which cannot be extened the appropriate length
                        # peptide on right end
                        if end >= orfLength:
                            with open('Peptides-removed.txt', 'a') as rm:
                                removed = gene + ',' + info['phosPep'] + ',' + str(idx) + ' '  +  'Right, past end\n'
                                rm.write(removed)
                            rm.close()
                            self.genesToRemove.append(gene)
                            continue
                        elif (start + i) - self.wing < 0:                # exclude peptides which are unable to extend at left end
                            with open('Peptides-removed.txt', 'a') as rm:
                                removed = gene + ',' + info['phosPep'] + ',' + str(idx) + 'Left, past start\n'
                                rm.write(removed)
                            rm.close()
                            self.genesToRemove.append(gene)
                            continue
                        ## extend peptide around phosphorylated amino acid site
                        phosPepPos = start + i - count
                        count += 1
                        extendPep = subject.seq[(phosPepPos - self.wing):(phosPepPos +1 + self.wing)]
                        # check for duplicate peptides from the same gene and remove duplicate, allow one to pass
                        # for example:  YOR051C_S33 and YOR051C_S33_T34 -- produce 3 peptides in total
                        # we remove one of the YOR051C_S33, so it is not given extra weight by motifx
                        name = gene.split('_')                         # example split gene name YOR051C_S33 
                        dupName = ('%s\t%s' %(name[0] ,str(extendPep) ))    #  'YOR051C  DEPSRESTPVRSQ'
                        # check if this pattern has already been seen
                        if dupName not in self.duplicatePeps:
                            info['extended'].append(extendPep)
                            self.duplicatePeps.add(dupName)
                        else:
                            with open('Peptides-removed.txt', 'a') as rm:
                                removed = ('%s\t DUPLICATE PEPTIDE FROM THE SAME GENE' %(dupName))
                                rm.write('%s\n' %(removed))        
                            rm.close()                          
                        
    def cleanResult(self):
        """ Remove gene names from self.pepInfo, for genes that have been excluded"""
        for gn in self.genesToRemove:
            if gn in self.pepInfo:
                del self.pepInfo[gn]

    def writeToFile(self):
        """ Write the extended peptides to a file, as a table, example:
        
        Gene    Group   Motif   originalPep     phosphoPep      Motifx_extended_pep
        YNL103W_T302    Induced         TTHTPNR TTHT*PNR        LHRTTHTPNRSSP
        YJL020C_S893_T894       Induced         RSTTHDVGEISNNVK RS*T*THDVGEISNNVK       VASIRRSTTHDVG,ASIRRSTTHDVGE

        """
        with open('motifx_pre_align_results.txt', 'w') as out:
            # write table header
            out.write("Gene\tGroup\tMotif\toriginalPep\tphosphoPep\tMotifx_extended_pep\n")
            for gene,info in self.pepInfo.items():
                zpep = []
                # merge the new extened peptides into a single comma delimited string
                for i in info['extended']:
                    zpep.append(str(i))
                extd = ','.join(zpep)
                out.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(gene, info['group'], info['motif'], info['origPep'], info['phosPep'], extd))
                    
                                   
def main():
    """
    Main, parse cmd line args and call writeRefFasta. 
    """
#******************************************************************************
# Command line args
#******************************************************************************    
    cmdparser = argparse.ArgumentParser(description="Pre-align peptides to an orf fasta file for use with Motifx website.",
                                        usage='%(prog)s -f <fasta file> -p <peptide info file>' ,prog='MotifxPreAlign.py' )
    cmdparser.add_argument('-f', '--fasta',action='store',     dest='FASTA',help='Optional: orf fasta input file (default SGD R64-2-1)', metavar='')
    cmdparser.add_argument('-i', '--info', action='store_true',dest='INFO', help='Print a more detailed description')
    cmdparser.add_argument('-p', '--pep',  action='store',     dest='PEP',  help='Required: Peptide info file', metavar='') 
    cmdparser.add_argument('-w', '--width',action='store',     dest='WID',  help='Optional: width of peptide, (13)', metavar='')            
    cmdResults = vars(cmdparser.parse_args())
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)

#******************************************************************************
# Check and define input variables
#******************************************************************************    
    # Print detailed program information to screen
    if cmdResults['INFO']:
        print("\n\tMotifxPreAlign.py\n")
        print("\tPurpose: Pre-align peptides to an orf fasta file for use with Motifx website.")
        print("\n\tRequired parameters:")
        print("\n\t-p peptide info file,")
        print("\n\t   Example, one per line of the following:")
        print("\n\t   YBR162W-A_T5,Induced,AVQtPR,AVQT*PR ")
        print("\t   YGL076C_T8_S11,Induced,AAEKILtPEsQLKK,AAEKILT*PES*QLKK")
        print("\n\tOptional parameters:\n")
        print("\t-f fasta file, expected to be an orf file.")
        print("\t-w width of peptides, all peptides will be the same length.")
        print("\n\tOutut: A list of peptides centered on the phosphorylated amino acid.")
        print("\tList is meant to be used as input for the motifx website.")
        print("\n\tTo see Python Docs and get a better explaination of the program:")
        print("\n\tOpen python console and enter")
        print("\timport sys")
        print("\tsys.path.append('/full/path/to/script')")
        print("\timport MotifxPreAlign")
        print("\thelp(MotifxPreAlign)")
        print("\n\tSee Mike Place for help or suggestions.\n")
        sys.exit(1)        
    
    # get fasta sequence file to query
    if cmdResults['FASTA'] is not None:
        fastaFile = cmdResults['FASTA']
        
    # get short sequences to find in the fasta sequences.
    if cmdResults['PEP'] is not None:
        pepFile = cmdResults['PEP']
    
    # check for alternate window size
    if cmdResults['WID'] is not None:
        window = int(cmdResults['WID'])
    else:
        window = 13
        
    pepd = MotifxPreAlign(fastaFile, pepFile, window)       # create object
    pepd.matchSeq()
    pepd.cleanResult()
    pepd.writeToFile()

    
if __name__ == "__main__":
    main()

