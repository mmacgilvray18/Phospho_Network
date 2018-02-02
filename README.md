# Computational phospho-proteomic network inference pipeline

The scripts enclosed were used to analyze wild type and mutant phosphoproteomic data as submitted in
MacGilvray et al.  The pipeline takes a list of phospho-peptides from S. cerevisiae and defines groups of likely co-regulated peptides that share the same phosphorylation motif (see manuscript Methods for details).  For each of these ‘modules’ of peptides, the pipeline then identifies ‘shared interactors’, defined as proteins from a background network of protein interactions that show more physical interactions with module constituent proteins then expected by chance.  Shared interactor-module pairs serve as inputs for a previously developed Integer Programming (IP) approach that connects the sources to their downstream target submodules (Chasman et al., 2014).

The pipeline titled **“Phospho_Network_Inference.ipynb”** consists of the scripts used to analyze data from MacGilvray et al.  The output is a **Network.sif file** compatible with Cytoscape for visualization and analysis.

**Most users** with wild-type S. cerevisiae phospho-proteomic data will be interested in a second pipeline
titled,**“Shared_interactors.ipynb”**.  This pipeline takes as input user-defined groups of phospho-peptides
(e.g. those with increased phosphorylation in response stimulus or those with decreased phosphorylation
in response to stimulus).  The first script partitions each group into ‘modules’ of phospho-peptides that
share the same sequence motif around the phosphorylation site.  The method then identifies Shared Interactors
as described above for each identified module.


Please see our _bioRxiv_ preprint for additional information:

[Network inference reveals novel connections in pathways regulating growth and defense in the yeast salt response](https://doi.org/10.1101/176230).
Matthew E. MacGilvray<sup>+</sup>, Evgenia Shishkova<sup>+</sup>, Deborah Chasman, Michael Place, Anthony Gitter, Joshua J. Coon, Audrey P. Gasch.
_bioRxiv_ 2017. doi:10.1101/176230

<sup>+</sup>These authors contributed equally to the work.

## Prerequisites
The user should define differentially changing phospho-peptides in the "WT" or "Parent" strain using their own criteria (eg; fold-change, p-value, etc.), followed by grouping/clustering phospho-peptides based on similar directionality of abundance change.

# Code Prerequisites:
 * This code has been run successfully on Mac OS X, Ubuntu 14, but should work on any
 unix like system.
 * Python 3 (version 3.4 although any version 3.3 and above should work)
   Anaconda is an easy way to install python, https://www.anaconda.com/download/
 * python libraries required: 

   Code has successfully executed using these versions:
  > python                    3.4.4   
  > biopython                 1.68               
  > jupyter                   1.0.0               
  > jupyter_core              4.1.0              
  > numpy                     1.11.0            
  > pandas                    0.19.2           
  > beautifulsoup4            4.4.1  (required for motfix.py)  
  > statsmodels               0.6.1  
  > rpy2                      2.7.8              

   once anaconda has been downloaded the required libraries may be installed using:

    conda install -c anaconda biopython=1.68
    conda install panda=0.19.2
    conda install numpy=1.11.0
    conda install jupyter=1.0.0
    conda install -c anaconda beautifulsoup4 
    conda install -c conda-forge statsmodels


## To run pipeline requires that you clone the git repository
  Assuming you are running linux, Open a terminal and enter:  

    git clone https://github.com/mmacgilvray18/Phospho_Network.git

  This will copy the git repository into a folder called Phospho_Network.

1. Use the provided ipython notebook, Phospho_Network_inference.ipynb (jupyter notebook).  

   cd to the repository, the notebook is setup to access files and directories from here.
   Enter the directory name and full path in the first code cell once you start the notebook.

   For a Jupyter Notebook tutorial see:
    https://www.datacamp.com/community/tutorials/tutorial-jupyter-notebook
    or
    https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/index.html

2. Using the terminal run 'jupyter notebook' in that directory. If you installed anaconda python
   you may have to run ~/anaconda/bin/jupyter notebook.
   This will open a web browser with a page showing the directory contents.  
   Click on Phospho_Network_inference.ipynb which opens another page showing the script. 

   It's Recommended way to run the cells in order,
   **DO NOT CHOOSE RUN ALL CELLS UNLESS YOU KNOW WHAT YOU ARE DOING.""
   The results of the first step **"Identify Mofifs"**, requires user intervention to create the
   input file for the next step.  The next 4 steps can be run automatically, these steps are:
   Identify Modules and Submodules -> Identify Shared Interactors 
   -> Preparation for Kullback-Leibler script -> Run Kullback-Leibler Module

   The next step is to run Kullback-Leibler 1000x with random shuffling. This step will be very 
   slow with a single process.  Use the -i 1000 (feel free to change, but should be greater than 1)
   -p 4 (sets the number of processes to use). Be smart don't try to use more processes(cores)
   than you have.  I recommend you use a process number that will evenly divide the number of 
   iterations.  

   Once that is complete you can auto run the rest of the cells.

3. Follow the script and change any input and output names as desired. 
   All output directories and files will be created in your current directory.
1. Use the provided ipython notebook, Phospho_Network.ipynb or Shared_interactors.ipynb depending
   on what you want to run.  These notebooks were generated using jupyter notebook.  

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
   
   >fasta            -- files for 'Create PWMs From Module Fasta' step<br>
   >Kullback-Leibler -- Kullback-Leibler Module to Each Kinase results<br>
   >KL-Shuffle       -- random shuffle of  Kullback-Leibler Module used
                        to generate a FDR (False Discovery Rate)<br>

   Primary File outputs:

    >Motifx results             -- Single text file and 3 directories
    >Submodules.txt             -- All identified submodules<br>
    >Shared_interactors.csv     -- All shared interactors 
    >position_weight_matrix.txt -- position weight matrix file<br>
    >Shared_interactors_Kinase_FDR.csv -- FDR scores for each Module <br>
    >Network.sif               -- SIF file for cytoscape<br>

## Scripts
************************************************************************
### MotifxPreAlign.py
```
Purpose: Pre-align peptides to an orf fasta file for use with Motifx website.
         This is called by Motifx.py, user doesn't need to access this script.
         But it must be available to Motifx.py.

 ```
************************************************************************
### Motifx.py
```
Purpose: Automate submitting jobs to Motif-x website.

Input: Plain text file listing files to process, one file name per line.

file format:
Input: Plain text file listing excel files to process, one excel file name per line.


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
Purpose: Identifies co-regulated groups of phospho-peptides 

Output: Submodules.txt  -- submodule constituent file

```
************************************************************************
### Identify Shared Interactors
```
Purpose: Identify proteins enriched for interactions with submodule constituent proteins
utilizing a background network of protein-protein interactions in yeast.
Kinase/phosphatase shared interactors are candidate submodule regulators, while other
proteins may represent members of complexes.

Output: Shared_interactors.csv - all shared interactors identified for submodule proteins
```
************************************************************************
### Prep for Kullback_Leiber Script 
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
kullback-Leibler.py -f 'position_weight_matrix.txt' -m '/PathTo/Mok_kinase_PWMs/

Purpose:  To quantify similarity between the Mok et al kinase PWMs and the module
PWMs, this script employs a previously described quantitative motif comparison method
called Kullback-Leibler divergence (KLD) (Thijs et al., 2002, Gupta et al., 2007). KLD
generates a similarity measure by comparing the Kullback-Leiber distance, or
information content, for each amino acid at each position between a query and
comparison PWM. The more alike two PWMs are, the closer to zero the score
approaches.

Output: A directory containing plain text .csv files named after each module (ie.
Induced_...sP..txt). Within the .csv files are 63 KLD scores representing how well the
63 Mok et al kinases match the module motif.
```
************************************************************************
### Kullback_Leibler Module toEachKinase Shuffled1000x
```
kullback-Leibler.py -f 'position_weight_matrix.txt' -m '/PathTo/Mok_kinase_PWMs/ -i 1000 -p 4 -o '/pathTo/KL-shuffle/'

Input: A plain text .csv file that contains all module position weight matrices. Each
module PWM should have 20 rows, representing each of the 20 naturally occurring
amino acids. They are in a column called "AA" which stands for amino acid. There
should also be 13 columns, labeled 0-12 (representing the 13 amino acid sequence length
of the phospho-peptides used to build the matrix) that contain the frequency of each
amino acid at each position.

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

Output: Shared_interactors_Kinase_FDR.csv - A table that contains for each module, 
all yeast kinases, including those found in the Mok et al dataset and those that were 
absent, and their FDR scores for each module. 
```
************************************************************************
### kullback-Leibler.py  
```
usage: Shuffle_kullback-Leibler.py -f <Position weight matrix file>  

Shuffle kullback-Leibler results for use w/ FDR function.

optional arguments:
  -h, --help          show this help message and exit
  -f , --file         position_weight_matrix.txt file,
  -m , --mokdir       Full path to Mok_kinase_PWMs directory
  -i , --iterations   Total number of iterations, this will be divided by
                      number of processes.
  -o , --out          Output directory
  -p , --processes    Number of processes to run, be smart don't use 
                      more than you have!

Output: A table that contains for each module, all yeast kinases, including those found in
the Mok et al dataset and those that were absent, and their FDR scores for each module.
Kinases not found in the Mok et al dataset are given an FDR score of 1.
```

************************************************************************
## phospho_subnet folder  
The [phospho_subnet](phospho_subnet) subdirectory contains code, data, and results for the ILP  
approach to infer the NaCl-dependent phosphoproteomic signaling subnetwork.  

## path_linker folder
The [path_linker](path_linker) subdirectory contains code, data, and results for
the comparison with the PathLinker network algorithm.  

## pcsf folder  
The [pcsf](pcsf) subdirectory contains code, data, and results for the  
comparison with the prize-collecting Steiner forest network algorithm.  

## python folder  
Initial scripts used to create the ipython notebook pipeline.  

## reference folder
orf_trans_all.20150113.fasta file for running motifx w/ a more recent reference.
bgNtWk.csv   
Background_Network.csv  

## required folder  
Annotated kinase and phosphatase background files.

Kinase_Groups.csv
kinase_phosphatase_yeast.csv

## References
Chasman D, Ho YH, Berry DB, et al. (2014) Pathway connectivity and signaling
coordination in the yeast stress-activated signaling network. Mol Syst Biol 10: 759.

Gupta S, Stamatoyannopoulos JA, Bailey TL & Noble WS (2007) Quantifying similarity
between motifs. Genome Biol 8: R24.

Mok J, Kim PM, Lam HY, et al. (2010) Deciphering protein kinase specificity through
large-scale analysis of yeast phosphorylation site motifs. Sci Signal 3: ra12.

Thijs G, Marchal K, Lescot M, Rombauts S, De Moor B, Rouze P & Moreau Y (2002) A Gibbs sampling method to detect overrepresented motifs in the upstream regions of
coexpressed genes. J Comput Biol 9: 447-464.
