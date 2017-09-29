import Bio
import os
from Bio.Seq import Seq
from Bio import motifs
from Bio.Alphabet import IUPAC

''' Script uses BioPython to generate position weight matrices from a directory containing Fasta files for
each modules phospho-peptides. Note: Duplicate amino acid sequences should be removed from the Fasta files before running this script, if they exist, to prevent overweigting
the matrix'''


alphabet = IUPAC.protein                                                           # use protein alphabet
instances = []
os.chdir("/Users/mmacgilvray18/Desktop/NaCl_Module_FASTA_Files_Dupes_Removed/")    # user defined directory containing Fasta files

def CreatePWM():
    ''' Function creates PWMs for each Module '''
    instances = []
    for x in os.listdir():                                                         # Iterate through the Fasta files in the directory
        if x.endswith('.txt'):
            with open(x, "r") as f:
                for line in f:
                    if line.startswith('>'):                                       
                        continue
                    line = line.rstrip()                                                 
                    instances.append(Seq(line, IUPAC.protein))                     # add amino acid sequence to instances
                m = motifs.create(instances)
                pwm = m.counts.normalize(pseudocounts = 1)                         # Add a +1 pseudocount
                instances = []
                
                
    
       
CreatePWM()


