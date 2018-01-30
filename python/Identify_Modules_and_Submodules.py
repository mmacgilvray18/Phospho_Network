import pandas as pd 


'''This script identifies co-regulated groups of phospho-peptides using the following approach:

1) First, the script identifies 'modules', which are groups of phospho-peptides that exhibit the same directionality in stress-dependent abundance change (ie, increased 'Induced', or decreased 'Repressed') and the same motif.
The module nomenclature is as follows: Induced/Repressed- motif (ex: Induced..RK.s....).

2) Next, the script partitions modules into 'submodules' based on their phospho-peptide constituents dependency on a protein(s) for stress-dependent abundance changes (ie, phospho-peptides that exhibit increased 'amplified' or decreased 'defective' abundance in a deletion strain
compared to the 'WT' or  'Parental' type strain). These phenotypes are user defined. If two or more mutant phenotypes are recorded for a phospho-peptide then it's placed into two separate subModules (one for each mutant phenotype). If there was not a mutant phenotype at a user
defined threshold then the phenotype is 'No-Phenotype'
The submodule nomenclature is as follows: module name-mutant phenotype/No-Phenotype (ex: Induced..RK.s....Mutant_Defective).

Possible submodule phenotypes: Induced-Defective, Induced-Amplified, Repressed-Defective, Repressed-Amplified, Induced-No-Phenotype, Repressed-No-Phenotype 

This script is designed to define submodules based on phenotypes for three mutant strains that were interogated by phospho-proteomics.  

The input file, which is in .csv format, must use the following format:

Column headers 
Ppep, Cluster, Motif, Peptide, FirstMutantNamePhenotype, SecondMutantNamePhenotype, ThirdMutantNamePhenotype

Note: Ppep stands for phospho-peptide and should be the YORF followed by the phosphorylated residue(s). Examples - YLR113W_S115, YLR113W_S115_T179
Cluster - either 'Induced' or 'Repressed'
Motif - Identified by Motif-X
Peptide - 13 amino acid long phospho-peptide returned by Motif-X. The middle residue is the phosphorylated amino acid.
FirstMutantNamePhenotype - column name is 'hog1' in this script, which is the name of a gene for which we interogated a deletion strain by phospho-proteomics
SecondMutantNamePhenotype - column name is 'pde2' in this script, which is the name of a gene for which we interogated a deletion strain by phospho-proteomics
ThirdMutantNamePhenotype - column name is 'cdc14' in this script, which is the name of a gene for which we interogated a deletion strain by phospho-proteomics
see example input file: Identify_Modules_and_Submodules_InputFile.csv

'''
# Import Input file and create a dataframe. 
Data=pd.read_csv('/Users/mmacgilvray18/Desktop/WT_Induced_Repressed_2Clusters_Motifs_Matrix_Hog1_Pde2_Phenotypes_0.4_0.6_Cdc14_Phen_0.2_0.4_Matrix_Final.csv') # Define path to input file



def Slicedataframe():
    '''Define a function that slices the input dataframe into independent dataframes based on the Cluster names. Next, slice these dataframes based on the presence of the same motif, generating 'modules' '''  
    ClusterLST=Data['Cluster'].unique().tolist()                            # generate a list of unique Cluster names (ie, 'Induced' and 'Repressed')
    lst=[]                                                                 
    DF=Data.copy()                                                         
    for cluster in ClusterLST:                                              # Select the first 'cluster' on the list 
        DF2=DF.loc[DF['Cluster']== cluster]                                 # Create a new dataframe by selecting only those rows that contain the selected 'cluster' in the 'Cluster' column 
        MotifLST=DF2['Motif'].unique().tolist()                             # From the newly created dataframe, place each instance of a unique motif into a list
        cleanedMotifLST = [x for x in MotifLST if str(x) != 'nan']          # drop the string 'nan' from the list. 'nan' occurs for Ppeps that did not have an identified Motif from Motif-X. 
        for motif in cleanedMotifLST:                                       # Select a motif in the list
            DF3=DF2.loc[DF['Motif']== motif]                                # Filter the dataframe, selecting only those rows that contain 'motif' in the Motif column
            DF3['freq'] = DF3.groupby('Motif')['Motif'].transform('count')  # Produce a new column, called 'freq' that contains the number of rows, and thus phospho-peptides, that contain a given motif.
            lst.append(DF3) 
        
    return lst

SlicedDF_lst=Slicedataframe()



def ConcatenateDFs():
    ''' Define a function that appends the dataframes in the SlicedDF_list together. '''
    EmptyDF = pd.DataFrame()                                                # create an empty dataframe
    for df in SlicedDF_lst:                                                
        df=df.copy() 
        EmptyDF=EmptyDF.append(df)                                          # append to the empty DF the dataframe selected and overwrite the empty dataframe
    return EmptyDF

Final_DF=ConcatenateDFs()

#--------------------------------------------------------------------------------------------------------------------------
FinalDFV2=Final_DF.fillna(0)                                                # fill any NaN values with '0'

def Module_Motif_NoMutantPhenotypeExists(df):
    ''' Define a function that assigns no-phenotype submodules'''
    if (df['hog1']==0) & (df['pde2']==0) & (df['cdc14']==0):
        return 'No_Phenotype_Exists'
    
FinalDFV2['Phenotype']=FinalDFV2.apply(Module_Motif_NoMutantPhenotypeExists, axis=1) 
FinalDFV2=FinalDFV2.loc[FinalDFV2['Phenotype']=='No_Phenotype_Exists']        # Select all rows for which "No_Phenotype_Exists" in the 'Phenotype' column.
FinalDFV2['subModule']=FinalDFV2.Cluster.map(str) + "_" +
FinalDFV2.Motif + "_" + FinalDFV2.Phenotype                                   # create a new column, called submodule, that contains the concatenated strings in the 'Cluster', 'Motif', and 'Phenotype' columns.

 
FinalDF=Final_DF.dropna(subset = ['hog1', 'pde2', 'cdc14'], how='all')        # Remove rows that have NaN in all 3 columns representing mutant phenotpes. This steps removes theNo-phenotype submodules which were creat                                                                               # ed above. 
FinalDF=FinalDF.fillna(0)                                                     # fill any NaN that remain with '0'
lstCols=['hog1', 'pde2', 'cdc14']                                             # make a list that contains the column headers for the 3 mutants. 




def DefineMutantContribution(row):
    ''' Define a function that identifies for each phospho-peptide if it has a phenotype in more than one mutant strain'''
    dictData={} 
    for colname in lstCols:    
        if not row[colname]==0:                                                # if value is not equal to zero, there is a mutant phenotype (ex; Induced_defective)
            dictData[colname]=row[colname]  
    if len(dictData.keys())==0: return 0  
    else:
        return ":".join(dictData.keys())
    
FinalDF['Contribution']=FinalDF.apply(lambda x: DefineMutantContribution(x), axis=1) 



def DefinePhenotypeFromMutants(row):
    ''' Define a function that captures the mutant phenotype for Ppeps with multiple phenotypes and places it within a column'''
    dictData={}  
    for colname in lstCols:   
        if not row[colname]==0:
            dictData[colname]=row[colname] 
    if len(dictData.keys())==0: return 0 
    else:
        return ":".join(dictData.values()) 
    
FinalDF['Phenotype']=FinalDF.apply(lambda x: DefinePhenotypeFromMutants(x), axis=1)

#--------------------------------------------------------------------------------------------------------------------------------------------------
''' Determine all Ppeps that have 2 or more mutant phenotypes (Pde2/Hog1/Cdc14 phenotypes), then the script produces a new column with individual subModule names for 
the hog1 phenotype''' 

FinalDF_multiplePhenotypes=FinalDF[FinalDF['Contribution'].str.contains(":")]     # Select 'contribution column rows that contain ":", which means the Ppep has two mutant phenotypes since this is a separator between gene names
FinalDF_multiplePhenotypes_hog1=FinalDF_multiplePhenotypes[FinalDF_multiplePhenotypes['Contribution'].str.contains("hog1")] 
FinalDF_multiplePhenotypes_hog1['Hog1']='hog1' 
FinalDF_multiplePhenotypes_hog1['subModule']=FinalDF_multiplePhenotypes_hog1.Cluster.map(str) + "_" + FinalDF_multiplePhenotypes_hog1.Motif + "_" + FinalDF_multiplePhenotypes_hog1.Hog1 + "_" + FinalDF_multiplePhenotypes_hog1.hog1


#--------------------------------------------------------------------------------------------------------------------------------------------------
''' Define all Ppeps that have 2 or more mutant phenotypes (Pde2/Hog1/Cdc14 phenotypes), then the script produces a new column with individual subModule names for 
the pde2 phenotype''' 
FinalDF_multiplePhenotypes_pde2=FinalDF_multiplePhenotypes[FinalDF_multiplePhenotypes['Contribution'].str.contains("pde2")]
FinalDF_multiplePhenotypes_pde2['Pde2']='pde2'
FinalDF_multiplePhenotypes_pde2['subModule']=FinalDF_multiplePhenotypes_pde2.Cluster.map(str) + "_" + FinalDF_multiplePhenotypes_pde2.Motif + "_" + FinalDF_multiplePhenotypes_pde2.Pde2 + "_" + FinalDF_multiplePhenotypes_pde2.pde2


#--------------------------------------------------------------------------------------------------------------------------------------------------
''' Define all Ppeps that have 2 or more mutant phenotypes (Pde2/Hog1/Cdc14 phenotypes), then the script produces a new column with individual subModule names for 
the cdc14 phenotype'''

FinalDF_multiplePhenotypes_cdc14=FinalDF_multiplePhenotypes[FinalDF_multiplePhenotypes['Contribution'].str.contains("cdc14")]
FinalDF_multiplePhenotypes_cdc14['Cdc14']='cdc14'
FinalDF_multiplePhenotypes_cdc14['subModule']=FinalDF_multiplePhenotypes_cdc14.Cluster.map(str) + "_" + FinalDF_multiplePhenotypes_cdc14.Motif + "_" + FinalDF_multiplePhenotypes_cdc14.Cdc14 + "_" + FinalDF_multiplePhenotypes_cdc14.cdc14


#--------------------------------------------------------------------------------------------------------------------------------------------------

'''This section of code appends the above mutant dataframes together (ie, FinalDF_multiplePhenotypes_cdc14, etc.) (contained ":"). The result is Ppeps with phenotypes in more than one strain are listed on multiple lines rather than a single line'''

FinalDF_mutants=FinalDF_multiplePhenotypes_cdc14.append(FinalDF_multiplePhenotypes_hog1) 
FinalDF_mutants_Final=FinalDF_mutants.append(FinalDF_multiplePhenotypes_pde2)

FinalDF_mutants_Final=FinalDF_mutants_Final[['Ppep','Cluster','Motif','Peptide','hog1','pde2','cdc14','freq','Contribution','Phenotype','subModule']] # Only retain these columns 


#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Drop from the original dataframe rows containing Ppeps with multiple mutant phenotypes. 
FinalDF_minus_multiPhenotypePpeps=FinalDF[FinalDF.Contribution.str.contains(":")==False] # Removing all rows that contain ":", and thus are phospho-peptides with multiple mutant phenotypes


# Generate the final submodule names
FinalDF_minus_multiPhenotypePpeps['subModule']=FinalDF_minus_multiPhenotypePpeps.Cluster.map(str) + "_" + FinalDF_minus_multiPhenotypePpeps.Motif + "_" + FinalDF_minus_multiPhenotypePpeps.Contribution + "_" + FinalDF_minus_multiPhenotypePpeps.Phenotype


#-----------------------------------------------------------------------------------------------------------------------------------------------------
Ppeps_with_PhenotypesDF=FinalDF_minus_multiPhenotypePpeps.append(FinalDF_mutants_Final)  # Appending together the dataframes that originally had single mutant phenotypes, and the dataframe that started with multiple mutant Phentoypes, but now contains single listings for each Ppep-mutant phenotype



# Remove any submodule that only has a single Ppep constituent, since by default a submodule must contain 2 Ppeps. 
Ppeps_with_PhenotypesDF_subModules=Ppeps_with_PhenotypesDF[Ppeps_with_PhenotypesDF.duplicated(['subModule'], take_last=True) | Ppeps_with_PhenotypesDF.duplicated(['subModule'])]  # only retain duplicates, get rid of single entries 



# Append to the dataframe with phenotype subModules, all No-Phenotype submodules
Ppeps_with_Phenotypes_subModules_and_noPhenotypes_DF=Ppeps_with_PhenotypesDF_subModules.append(FinalDFV2) # append to dataframe
Ppeps_with_Phenotypes_subModules_and_noPhenotypes_DF=Ppeps_with_Phenotypes_subModules_and_noPhenotypes_DF[['Ppep', 'Cluster', 'Motif', 'Peptide', 'hog1', 'pde2', 'cdc14', 'freq', 'Contribution', 'Phenotype', 'subModule']]

#-----------------------------------------------------------------------------------------------------------------------------------------------------
''' Create a column with the 'Module' name '''
Ppeps_with_Phenotypes_subModules_and_noPhenotypes_DF['Module']=Ppeps_with_Phenotypes_subModules_and_noPhenotypes_DF.Cluster.map(str) + "_" + Ppeps_with_Phenotypes_subModules_and_noPhenotypes_DF.Motif 

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# define a function that will write out a dataframe as a tab separated file
def Dataframe_to_Tsv (dataframe, NewFileName):
    path ="/Users/mmacgilvray18/Desktop/"
    dataframe.to_csv (path+NewFileName,sep='\t')

Dataframe_to_Tsv(Ppeps_with_Phenotypes_subModules_and_noPhenotypes_DF, 'subModules_with_and_without_Phenotypes_hog1_pde2_0.4_0.6_cdc14_0.2_0.4_Cutoffs.csv') 
# The above file contains all modules and subModules with and without mutant phenotypes. 

