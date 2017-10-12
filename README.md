# Computational phospho-proteomic network inference pipeline

This pipeline turns a list of _S. cerevesiae_ phospho-peptides that exhibit stress responsive
abundance changes, as measured by mass spectrometry, into a hierarchical signaling
network, connecting upstream kinases and phosphatases to their downstream targets. Our
computational pipeline is based on the premise that kinases and phosphatases recognize
target substrates through specific amino acid sequences at the phosphorylated residue,
called phosphorylation motifs. This pipeline groups phospho-peptides with similar
abundance changes and the same phosphorylation motif into modules. Modules are
partitioned into smaller groups, called submodules, based on differences in phospho-
peptide abundance in mutant strain(s) (sources). Candidate submodule regulators, called
shared interactors, are identified through enrichment analysis using a protein interaction
network in yeast (Chasman et al., 2014). Shared interactor-submodule pairs serve as
inputs for a previously developed Integer Programming (IP) approach that connects the
sources to their downstream target submodules (Chasman et al., 2014).

Please see our _bioRxiv_ preprint for additional information:

[Network inference reveals novel connections in pathways regulating growth and defense in the yeast salt response](https://doi.org/10.1101/176230).
Matthew E. MacGilvray<sup>+</sup>, Evgenia Shishkova<sup>+</sup>, Deborah Chasman, Michael Place, Anthony Gitter, Joshua J. Coon, Audrey P. Gasch.
_bioRxiv_ 2017. doi:10.1101/176230

<sup>+</sup>These authors contributed equally to the work.

## Prerequisites
The user should define differentially changing phospho-peptides in the "WT" or "Parent"
strain using their own criteria (eg; fold-change, p-value, etc.), followed by
grouping/clustering phospho-peptides based on similar directionality of abundance
change.

# Code Prerequisites:
 * This code has been run successfully on Mac OS X, Ubuntu 14
 * Python 3 (version 3.4 although any version 3 should work)
   Anaconda is an easy way to install python, https://www.anaconda.com/download/
 * python libraries required:  Pandas, numpy, Biopython, jupyter, beautifulsoup
   
   Code has successfully executed using these versions:
  > python                    3.4.4
  > biopython                 1.68               
  > jupyter                   1.0.0             
  > jupyter_core              4.1.0            
  > numpy                     1.11.0          
  > pandas                    0.19.2         
  > beautifulsoup4            4.4.1  (required for motfix.py)
   
   once anaconda has been downloaded the required libraries may be installed using:

    conda install -c anaconda biopython=1.68 

    conda install panda=0.19.2
   
    conda install numpy=1.11.0
    
    conda install jupyter=1.0.0

    conda install -c anaconda beautifulsoup4 


## To run pipeline requires that you clone the git repository
  Assuming you are running linux: 
  Open a terminal and type:  git clone https://github.com/mmacgilvray18/Phospho_Network.git
  This will copy the git repository into a folder called Phospho_Network.

1. Use the provided ipython notebook, Phospho_Network.ipynb.  This notebook was
   generated using jupyter notebook.  

   cd to the repository, the notebook is setup to be run from that directory.

2. Using the terminal run 'jupyter notebook' in that directory. If you installed anaconda python
   you may have to run ~/anaconda/bin/jupyter notebook.
   This will open a web browser w/ a page showing the directory contents.  
   Click on Phospho_Network.ipynb which opens another page showing the script.  

3. Follow the script and change any input and output names as desired. Gene Names will definitely
   need to be changed in the Identify Modules & Submodules step. All output directories and files
   will be created in your current directory.

   The code will create all required intermediate files.

   Directory outputs:
   
   >FastaFiles_Modules  -- peptide input files for 'Create PWMs From Module Fasta' step<br>
   >ClassA_NoShuffle_KL -- Kullback-Leibler Module to Each Kinase results<br>
   >Shuffle_KL          -- random shuffle of  Kullback-Leibler Module to Each Kinase used
                           to generate a FDR (False Discovery Rate)<br>

   Primary File outputs:

    >Network_Submodule_Nodes_background_Network.csv -- All identified submodules<br>
    >Final_enriched.csv          -- Filtered set of identified submodules, strict criteria<br>
    >SIs_submodule_Relationships_Define_ClassA_Network.csv -- Classify Shared Interactors Outputs<br>
    >position_weight_matrix.txt  -- position weight matrix created from fasta input<br>
    >FDR_Scores.csv              -- Calculate FDR Each Module step results<br>
    >FINAL_SIF.csv               -- SIF file for cytoscape<br>


## Scripts
************************************************************************
### MotifxPreAlign.py
```
Purpose: Pre-align peptides to an orf fasta file for use with Motifx website.

Required parameters:

        -p peptide info file

 Example, one per line of the following:

  YBR162W-A_T5,Induced,AVQtPR,AVQT*PR
  GL076C_T8_S11,Induced,AAEKILtPEsQLKK,AAEKILT*PES*QLKK

 Must set path to location where script is housed on your machine.

 Optional parameters:

        -f fasta file, expected to be an orf file.
        -w width of peptides, all peptides will be the same length.

 Output: A list of peptides centered on the phosphorylated amino acid.
 List is meant to be used as input for the motifx website
 ```
************************************************************************
### Motifx.py
```
Purpose: Automate submitting jobs to Motif-x website.

Input: Plain text file listing excel files to process, one excel file name per line.

Excel file format:
Ppep    Group   Localized_Sequence  Motif_X_Input_Peptide
YGL076C_T8_S11  Induced AAEKILtPEsQLKK  AAEKILT*PES*QLKK

Column order is unimportant, column names must match above.

usage: Motifx.py -f inputfiles
Must set path to location where script is housed on your machine.

 Optional Arguments:
  -o Minimum number of times each of your extracted motifs to occur in the data set (10)
  -s P-value threshold for the binomial probability (.000001)
  -u upload your own version of SGD proteome (orf_trans.fasta).
  -w Number of total characters in motif, (13)


 Output:
 A final text table named after the input file containing all the motifs matched to a gene.
 Given an input file named, motifx_sample.xlsx the final results file will be:
motifx_sample-Motifx-results.txt

 The other results are put in a directory.
 For instance if your input file is called motifx_sample.xlsx

 3 directories will be created one for each central character:

 motifx_sample_T    motifx_sample_S    motifx_sample_Y

These contain the LOGO pngs and the original html results page.
```
************************************************************************
## Primary Steps in Notebook Pipeline
### Identify Modules and Submodules
```
Purpose: Group phospho-peptides into modules and submodules based on stress-
dependent phosphorylation changes, motifs, and the presence mutant phenotypes.

Input: Plain csv text file

Csv format:
Ppep, Cluster, Motif, Peptide, FirstMutantNamePhenotype,
SecondMutantNamePhenotype, ThirdMutantNamePhenotype

YAL019W_T131, Induced, .sP., SPKKSETPPPVKK, Induced_Defective, empty,
Repressed_Amplified

Column order is unimportant, column names must match above.
The Cluster and MutantNamePhenotype columns are user defined and can be changed in
the script.

Output:
A final csv file in table format that contains on each row a phospho-peptide and the
module and submodule(s) it resides within.
```
************************************************************************
### Identify Shared Interactors
```
Purpose: Identify proteins enriched for interactions with submodule constituent proteins
utilizing a background network of protein-protein interactions in yeast.
Kinase/phosphatase shared interactors are candidate submodule regulators, while other
proteins may represent members of complexes.

Input: 4 plain csv text files

Submodule_constituents.csv (user created)
Csv format:
Submodule, ORF
Induced_....sP._hog1_induced_defective,YMR031C

Background_Network.csv   (provided on Github)
Csv format:
Protein1,Protein2,Interaction,Directed
YCR066W,YBR088C,kinase_substrate,1

Number_Interactions_Each_Protein.csv    (provided on Github)
Csv format:
Protein,Total
YPL031C,419

Annotation.csv            (provided on Github)
Csv format:
systematic_name_dash_removed,systematic_name
YKL135WA,YKL135W-A

Column order is unimportant, column names must match above.
User can modify the FDR cutoff.

Output:
-A table containing all shared interactors, without enrichment analysis, and their
interactions with submodule proteins.
-A table containing all enriched shared interactors identified for submodule proteins.
```
************************************************************************
### Classify Shared Interactors Inputs Outputs
```
Purpose: Classify shared interactors based on their interaction(s) directionality with
submodule constituent proteins. "Input" shared interacotrs have at least one directed
interaction (ie, kinase-substrate) toward a submodule constituent protein or one non-
directed interaction. Shared interactors that fail to meet these criteria are classified as
"outputs" and are unlikely to regulate submodules. Input shared interactor kinases and
phosphatases are likely candidate submodule regulators.

Input: plain csv text file

Csv format:
SI_submodule,Shared_Interactor,SI_name,Motif_Containing_Proteins,submodule_Name
,Interaction_Directionality

YLR164C_Repressed_..RR.s.No_Phenotype_Exists,YLR164C,Tpk1,YDR207C,
Repressed_..RR.s.No_Phenotype_Exists, kinase_substrate

Column order is unimportant, column names must match above.

Output: A table containing each SI-submodule pair and if the interaction
is an "Input" or "Output"
```
************************************************************************
### Create Fasta File
```
Purpose: Turn a list of modules and their phospho-peptide amino acid sequences into a
Fasta file that will be used later in the pipeline to generate position weight matrices
(PWMs).

Input: plain csv text file

Csv format:
Module,Name,Sequence
Induced_...sP.,YDR497C_S55,IQRAPASDDEDRI

Column order is unimportant, column names must match above.

Required parameters:
All peptide sequences should be the same length (13 amino acids). Module constituents
should be used here, not submodules. User must modify the directory where the output
Fasta Files for each module will be deposited.

Output: A directory containing a Fasta File for each module and its phospho-peptide
constituents.
```
************************************************************************
### Create PWMs From Module Fasta
```
Purpose: Generate PWMs for each module, using the module Fasta files. Module PWMs
can then be compared to PWMs for 63 known kinase recognition motifs (Mok et al.,
2010).

Input: A directory containing plain .txt files in Fasta format.

Required parameters: User must specify the directory where the
input Fasta Files for each module will be deposited. Results will print to screen and can
be moved to a .txt file manually.

Output: A PWM for each module that indicates the frequency of each amino acid at each
position.
```
************************************************************************
### Kullback Leibler Module toEachKinase
```
Purpose:  To quantify similarity between the Mok et al kinase PWMs and the module
PWMs, this script employs a previously described quantitative motif comparison method
called Kullback-Leibler divergence (KLD) (Thijs et al., 2002, Gupta et al., 2007). KLD
generates a similarity measure by comparing the Kullback-Leiber distance, or
information content, for each amino acid at each position between a query and
comparison PWM. The more alike two PWMs are, the closer to zero the score
approaches.

Input: A plain text .csv file that contains all module position weight matrices. Each
module PWM should have 20 rows, representing each of the 20 naturally occurring
amino acids. They are in a column called "AA" which stands for amino acid. There
should also be 13 columns, labeled 0-12 (representing the 13 amino acid sequence length
of the phospho-peptides used to build the position weight matrix) that contain the
frequency of each amino acid at each position.

Csv file format
Motif,AA,0,1,2,3,4,5,6,7,8,9,10,11,12
Induced_...sP.,P:,0.05,0.05,0.03, 0.05,0.05,0.03,0.05,0.05,0.03, 0.05,0.05,0.03

In addition, a directory that contains the Mok et al kinase PWMs. They have the identical
format as above. They have been pre-generated and are available for download on Github.
The repository is titled, "Mok_kinase_PWMs"

Output: A directory containing plain text .csv files named after each module (ie.
Induced_...sP..txt). Within the .csv files are 63 KLD scores representing how well the
63 Mok et al kinases match the module motif.
```
************************************************************************
### Kullback_Leibler Module toEachKinase Shuffled1000x
```
This is slow, for a faster command line version use: Shuffle-kullback-Leibler.py

Purpose:  To quantify similarity between the Mok et al kinase PWMs and the module
PWMs, this script employs a previously described quantitative motif comparison method
called Kullback-Leibler divergence (KLD) (Thijs et al., 2002, Gupta et al., 2007). KLD
generates a similarity measure by comparing the Kullback-Leiber distance, or
information content, for each amino acid at each position between a query and
comparison PWM. The more alike two PWMs are, the closer to zero the score
approaches. 1000 Shuffles of the Mok Kinase PWMs are performed by the script,
generating randomized PWMs that are compared against the Module PWMs, producing a
distribution of scores.

Input: A plain text .csv file that contains all module position weight matrices. Each
module PWM should have 20 rows, representing each of the 20 naturally occurring
amino acids. They are in a column called "AA" which stands for amino acid. There
should also be 13 columns, labeled 0-12 (representing the 13 amino acid sequence length
of the phospho-peptides used to build the matrix) that contain the frequency of each
amino acid at each position.

Csv file format
Motif,AA,0,1,2,3,4,5,6,7,8,9,10,11,12
Induced_...sP.,P:,0.05,0.05,0.03, 0.05,0.05,0.03,0.05,0.05,0.03, 0.05,0.05,0.03

In addition, a directory that contains the Mok et al kinase PWMs. They have the identical
format as above. They have been pre-generated and are available for download on Github.
The repository is titled, "Mok_kinase_PWMs"

Output: A directory containing plain text .csv files named after each module. Within
the .csv files are 63,000 KLD scores representing how well the 63 Mok et al kinases
match the module motif after 1000 permutations of each Mok kinase.
```
************************************************************************
### CalculateFDR Each Module to Each Kinase
```
Purpose: Identify FDR scores for each Mok et al kinase and each module by comparing
the non-shuffled scores to the distribution of shuffled scores. The user can then manually
define the FDR cutoff to call kinases "motif-match" or "non-match" for a given module.

Input: Two directories and a single plain text .csv file, called "Kinases_Not_In_Mok.csv"
that is provided on the Github page. The first directory contains plain text .csv files with
KLD scores for non-shuffled Mok et al kinases and Modules. The second directory
contains plain text .csv files containing KLD scores for shuffled Mok et al kinases and
Modules.

Csv format (For both Input Directories)

Scores,Kinase,Module,
13.25,cdc15,Induced.sP.

Output: A table that contains for each module, all yeast kinases, including those found in
the Mok et al dataset and those that were absent, and their FDR scores for each module.
Kinases not found in the Mok et al dataset are given an FDR score of 1.
```
************************************************************************
### Shuffle_kullback-Leibler.py

    Usage: Shuffle_kullback-Leibler.py -f <Position weight matrix file>  

    Shuffle kullback-Leibler results for use w/ FDR function.

    optional arguments:<br>
        -h, --help            show this help message and exit<br>
        -f , --file           Position Weight Matrix file, from Matt's script<br>
        -i ITER, --iterations Total number of iterations, this will be divided by number of processes.<br>
        -p PROC, --processes  Number of processes to run, be smart don't use more than you have!<br>                                                                      

************************************************************************
## phospho_subnet
The [phospho_subnet](phospho_subnet) subdirectory contains code, data, and results for the ILP
approach to infer the NaCl-dependent phosphoproteomic signaling subnetwork.

## pcsf
The [pcsf](pcsf) subdirectory contains code, data, and results for the
comparison with the prize-collecting Steiner forest network algorithm.


## References
Chasman D, Ho YH, Berry DB, et al. (2014) Pathway connectivity and signaling
coordination in the yeast stress-activated signaling network. Mol Syst Biol 10: 759.

Gupta S, Stamatoyannopoulos JA, Bailey TL & Noble WS (2007) Quantifying similarity
between motifs. Genome Biol 8: R24.

Mok J, Kim PM, Lam HY, et al. (2010) Deciphering protein kinase specificity through
large-scale analysis of yeast phosphorylation site motifs. Sci Signal 3: ra12.

Thijs G, Marchal K, Lescot M, Rombauts S, De Moor B, Rouze P & Moreau Y (2002) A Gibbs sampling method to detect overrepresented motifs in the upstream regions of
coexpressed genes. J Comput Biol 9: 447-464.
