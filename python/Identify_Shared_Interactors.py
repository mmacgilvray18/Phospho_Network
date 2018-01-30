import pandas as pd
import numpy as np
from scipy.stats import hypergeom


''' This script identifies proteins enriched for interactions with Submodule constituent proteins, based on known interactions in the background network. We call these proteins 'Shared Interactors'. The background network is a protein
interaction network curated in yeast under mostly nutrient replete conditions that contains 4638 proteins and ~ 25,000 interactions, including directed (ex; kinase-substrate), and 
non-directed. 

Proteins enriched for interactions with Submodule proteins at a 5% FDR, determined by a hypergeometric test and BH correction, are considered shared interactors.
Shared Interactors represent numerous functional classes, including kinases and phosphatases. Kinase and phosphatase shared interactors represent potential Submodule regulators.
 
HyperG function:
distrib=hypergeom(N,M,n)
distrib.pmf(m)

N - population size  (4638 unique proteins in the Background network file - phospho_v4_bgnet_siflike_withdirections_Matt_Modified.csv)
M - total number of successes  (# of interactions for a given protein. ie. Protein A has 200 known interactions in the background network).
n - the number of trials (also called sample size) -  ie. (Number of proteins that reside within a submdoule)
m - the number of successes - for example: Protein A, a shared interactor, has 35 interactions with proteins in Submodule B. 

'''

Submodule_DF = pd.read_csv('/Users/mmacgilvray18/Desktop/Submodule_constituents.csv')                                                                       # File that contains Submodule names and their protein constituents
#print (Submodule_DF.head(5))
BgNet= pd.read_csv('/Users/mmacgilvray18/Desktop/Background_Network.csv')                                                                                   # Background network of protein interactions
#print (BgNet.head(5))
Num_Prot_Inter = pd.read_csv ('/Users/mmacgilvray18/Desktop/Number_Interactions_Each_Protein.csv')                                              # Number of protein interactions for each protein in the background network
Annotation_DF = pd.read_csv('/Users/mmacgilvray18/Desktop/Annotation_dashes_removed_for_SI_renaming.csv')                                                   # Yeast protein annotation file

  
Submodule_List=Submodule_DF['Submodule'].unique().tolist()                                                                                                  # Send the Submodules to a list, but filter out duplicates, which there will be many, since the Submodules will have been found in many proteins.


dicOrfs={}
for Submodule in Submodule_List:                                                                                                                            # Key (Submodule), Value (Yeast ORFs that are Submodule constituents). Filter ORFs found twice to single occurence (important for enrichment analysis)
    dicOrfs[Submodule]=(Submodule_DF.loc[Submodule_DF['Submodule'] == Submodule])['ORF'].unique().tolist()
        

dicOrfsCounts={}  
for k,v in dicOrfs.items():  
    if k not in dicOrfsCounts:  
        value=len(v)            
        dicOrfsCounts[k]=value


        
df_Submodule_Size=pd.DataFrame(list(dicOrfsCounts.items()),                                                                                                  # convert dict to dataframe.
                      columns=['Submodule','n'])


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

def DF_to_TSV(dataframe, NewFileName): 
    ''' Define a function that writes out a dataframe as TSV'''
    path ='/Users/mmacgilvray18/Desktop/' 
    dataframe.to_csv (path+NewFileName,sep='\t')

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
def SliceDataframe():
    ''' For each Submodule identify all proteins that interact with the Submodule proteins in the backgroudn network '''
    lst = []
    for key in dicOrfs.keys():                                                                                                                             #Select the key, which is a Submodule, from the dict
        CurrentDF=BgNet.copy() 
        x=CurrentDF[CurrentDF['Protein1'].isin(dicOrfs[key])].rename(columns={'Protein1':'Submodule_Containing_Proteins', 'Protein2':'Possible_Shared_Interactors'})                              #Create a new dataframe that is a slice of the salt background network, and only contains proteins that were passed in "dicOrfs[key]". At the same time, rename the columns                                
        x['Submodule']=key 
        lst.append(x)
        
    return lst

Sliced_dataframe_list= SliceDataframe()
    
      
def Add_n():    
    ''' Function adds 'n', the number of proteins in the Submodule, to each dataframe'''
    lst= []
    for df in Sliced_dataframe_list:
        NewDF=df.merge(df_Submodule_Size)
        lst.append(NewDF)
        
    return lst

Sliced_dataframe_list= Add_n()

            
def Identify_Shared_Interactors():
    ''' Function identifies proteins that interact with at least 2 protein constituents of each submodule'''
    
    lst=[] 
    for df in Sliced_dataframe_list: 
        NewDF=df.copy()
        NewDF2=NewDF[NewDF.duplicated(['Possible_Shared_Interactors'], take_last=True)| NewDF.duplicated(['Possible_Shared_Interactors'])]                  # Only retain proteins that interact with at least 2 submodule protein constituents
        x=NewDF2.sort('Possible_Shared_Interactors', ascending=True) 
        lst.append(x)
       
    return lst

Shared_Interactors_lst=Identify_Shared_Interactors()


def AppendDFs_that_Contain_AllSharedInteractors_and_their_targets():
    ''' Function appends all submodules and their shared interactors together into a single file'''
    EmptyDF = pd.DataFrame() 
    for df in Shared_Interactors_lst:  
        df=df.copy() 
        EmptyDF=EmptyDF.append(df)
    return EmptyDF

SI_andTargets=AppendDFs_that_Contain_AllSharedInteractors_and_their_targets()


SI_andTargets_FINAL=pd.merge(left=SI_andTargets, right=Annotation_DF, how='left',
                              left_on='Possible_Shared_Interactors', right_on='systematic_name_dash_removed')                                               # complete a merge so I can get the dashes back in the names, which are not included in the background network
del SI_andTargets_FINAL['Possible_Shared_Interactors']                                                                                                      # drop because  lacks the dashes which are needed for the correct naming convention
del SI_andTargets_FINAL['systematic_name_dash_removed']                                                                                                     # drop because carried over from the merge
del SI_andTargets_FINAL['Directed']

SI_andTargets_FINAL.columns = ['Submodule_Containing_Proteins', 'Interaction', 'Submodule', 'n','Possible_Shared_Interactors']                        # rename columns



DF_to_CSV(SI_andTargets_FINAL, 'SI_Identification_NaCl_SubmoduleS_Possible_SIs_and_Targets_Dashes_Removed_4638_proteins_01_18_17_0.1_FDR.csv')              # All interactions between SIs and their submodule constituent proteins. No enrichment at this step.


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
''' Preparing dataframe for Hypergeometric test'''

def Add_N_and_m():
    ''' Function adds 'N' and calculates 'm' values, which are inputs for the hypergeometric test, to the datframe'''
    lst=[]
    for df in Shared_Interactors_lst:
        NewDF=df.copy()
        NewDF['N'] = 4638                                                                                                                                   # of proteins in the background network
        NewDF['m'] = NewDF.groupby('Possible_Shared_Interactors')['Possible_Shared_Interactors'].transform('count')
        lst.append(NewDF)
    
    return lst

Dataframes_list_with_n_N_m=Add_N_and_m()


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------


def Drop_dups():
    ''' For each dataframe, which contains a single submodule, it's protein constituents, and shared interactors, drop duplicate entries for identified SI proteins
    . This leaves a single entry for each shared interactor protein. '''
    lst=[]
    for df in Dataframes_list_with_n_N_m:
        NewDF=df.copy()
        Final_DF=NewDF.drop_duplicates('Possible_Shared_Interactors')
        Final_DF=Final_DF.rename(columns={'Possible_Shared_Interactors':'Shared_Interactor'})
        lst.append(Final_DF)
        
    return lst

Drop_Dups_lst=Drop_dups()


def Return_M():
    ''' Function identifies 'M' (the total number of interactions for each Shared Interactor protein in the background network) and adds that number
    to the dataframe'''
    lst=[]
    for df in Drop_Dups_lst:
        NewDF=df.copy()
        NewDF2=df.copy()
        NewDF_lst=NewDF['Shared_Interactor'].tolist()                                                                                                            # place all proteins in the 'Shared_Interactor' column in a list 
        Shared_Interactors=Num_Prot_Inter[Num_Prot_Inter['Protein'].isin(NewDF_lst)].rename(columns={'Protein':'Shared_Interactor', 'Total':'M'})
        Shared_Interactor_merge=Shared_Interactors.merge(NewDF2, on='Shared_Interactor')
        Shared_Interactor_merge=Shared_Interactor_merge.sort('Shared_Interactor', ascending=True)
        lst.append(Shared_Interactor_merge)
        
    return lst

Return_M_lst=Return_M()

#-----------------------------------------------------------------------------------------------------------------------------------------
def hyper(N,M,n,m): 
    ''' Function defines the parameters for a hypergeometric test that returns a p-value representing the chances of identifying >= x, where x is the number of successes '''  
    frozendist=hypergeom(N,M,n)
    ms=np.arange(m, min(n+1, M+1))
    rv=0;
    for single_m in ms: rv=rv+frozendist.pmf(single_m)
    return rv

def run_hyper():
    ''' Function calls the hypergeometric function above  on each shared interactor for each submodule'''
    lst=[]
    for df in Return_M_lst:
        if not df.empty:
            NewDF=df.copy()
            NewDF['p-value'] = NewDF.apply(lambda row: hyper(row['N'], row['M'], row['n'], row['m']), axis=1)
            lst.append(NewDF)
        
    return lst 

run_hyper_lst=run_hyper()

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def AppendDFs():
    ''' Append DFs for each submodule and it's SIs together into a single DF'''   
    EmptyDF = pd.DataFrame() #
    for df in run_hyper_lst: 
        df=df.copy() 
        EmptyDF=EmptyDF.append(df)
    return EmptyDF

Final=AppendDFs()

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
''' Prepping for Benjamini Hochberg procedure. Below code is ranking p-values from 1 to n based on lowest to highest p-value score'''

Final=Final.sort(['p-value'],ascending=[True])                                                                                              # Sort p-values from lowest to highest
Final_resetIndex=Final.reset_index()                                                                                                        # Reset the index after the sort
Final_resetIndex.index +=1                                                                                                                  # start numbering at 1 for index

       
NewDF=Final_resetIndex
NewDF_Allp_values=Final_resetIndex
NewDF=NewDF[['p-value']]                                                                                                                    # select only the p-value column of the dataframe 
NewDF_dropdups=NewDF.drop_duplicates('p-value')                                                                                             # drop duplicate p-values
NewDF_dropdups=NewDF_dropdups.reset_index()                                                                                                 # reset the index
NewDF_dropdups.index +=1                                                                                                                    # start numbering at 1 for index
NewDF_dropdups['Rank(i)'] = NewDF_dropdups.index                                                                                            # #Add a rank column that will be filled with index values. 
NewDF_dropdups=NewDF_dropdups.drop('index', 1)                                                                                              # Drop the additional column 'index' that is not sorted.
NewDF_merge=NewDF_Allp_values.merge(NewDF_dropdups, on='p-value')                                                                           # create a new dataframe that is a merge of the dataframe with all p-values, and the dataframe with unique p-values and their ranks. 
NewDF_merge=NewDF_merge.drop('index',1)                                                                                                     # drop the index that was added from the merge. This leaves all p-values ordered from lowest to highest with their ranking.
     



'''Add parameters necessary for completing Benjamini-Hochberg procedure '''

NewDF=NewDF_merge
NewDF['m_(number_of_tests)']=(len(NewDF))                                                                                                   # Add 'm (number of tests)' column 
NewDF['Q_(FDR)']=0.05                                                                                                                       # Add Q (FDR) column. This can be changed manually.
NewDF['(i/m)Q']=((NewDF['Rank(i)']/NewDF['m_(number_of_tests)'])*NewDF['Q_(FDR)'])                                                          # add the (i/m)Q column 
NewDF['BH_significant']=NewDF.apply(lambda x: 1 if x['p-value']<x['(i/m)Q'] else 0, axis=1)                                                 # Identify which proteins are  significant. 
NewDF=pd.merge(left=NewDF, right=Annotation_DF, how='left', left_on='Shared_Interactor', right_on='systematic_name_dash_removed')           # complete a merge to recover dashed version of YORFs
del NewDF['Shared_Interactor'] 
del NewDF['systematic_name_dash_removed']
del NewDF['Directed']
NewDF.columns = ['M','Submodule_Containing_Proteins', 'Interaction', 'Submodule', 'n','N','m','p-value','Rank(i)', 'm_(number_of_tests)', 'Q_(FDR)','(i/m)Q','BH_significant', 'Shared_Interactor'] # rename columns


DF_to_CSV(NewDF, 'SIs_NaCl_Network_SubmoduleS_01_18_17_FDR_0.05_BH_done_At_Once_sig_and_not_Dashes_Removed_4638_Nodes_background_Network.csv') # Write out final file with enriched shared interactors for each submodule








