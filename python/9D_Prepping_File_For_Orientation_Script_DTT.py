import pandas as pd 

''' 

The purpose of this script is to identify for each enriched SI (according to the hypergeometric and BH correction) what it's interacting proteins, and potential "targets"  are in the background network. Next, the script identifies the correct orientation for interactions between a Shared Interactor and it's partners, which are reversed in an earlier script. Thus, the script ensures all interactions are in the correct orientation. 

The output file from this script can then be used to determine if a Shared Interactor is acting as an input (all interactions going towards a subModule, or is an output (interactions go away), etc. An input is a potential regulator of a submodule whereas an output is unlikely to regulate a submodoule.

'''
##################################################################################################################################

Enriched_SIs_subModules_Only_Sig=pd.read_csv('/Users/mmacgilvray18/Desktop/DTT_T120_SIs_All_Sept2017_Enriched.csv')  # Import ONLY enriched SIs (drop non-enriched SIs)

Enriched_SIs_subModules_Only_Sig["SI_subModule"] = Enriched_SIs_subModules_Only_Sig["Shared_Interactor"].map(str) + "_" + Enriched_SIs_subModules_Only_Sig["subModule"] # creating a new column by merging two columns together.



All_enriched_and_not_SIs_andTargets=pd.read_csv('/Users/mmacgilvray18/Desktop/SI_Identification_T120_DTT_Sept2017_T120_Possible_SIs_and_Targets_Dashes_Removed_4638_proteins_BOTH_Reps_Ppeps_NORMALIZED.csv') # Import all identified Shared Interactors (both enriched and not) and their interacting partners from submodules.

All_enriched_and_not_SIs_andTargets["SI_subModule"] = All_enriched_and_not_SIs_andTargets["Possible_Shared_Interactors"].map(str) + "_" + All_enriched_and_not_SIs_andTargets["subModule"]

Merged_left = pd.merge(left=Enriched_SIs_subModules_Only_Sig,right=All_enriched_and_not_SIs_andTargets, how='left', left_on='SI_subModule', right_on='SI_subModule') # 

Merged_left_FINAL=Merged_left[['SI_subModule', 'Shared_Interactor', 'Motif_Containing_Proteins_y', 'subModule_y']] # retain these columns only
Merged_left_FINAL.columns=['SI_Module','Shared_Interactor', 'Motif_Containing_Proteins', 'subModule_Name'] # rename columns



def DF_to_CSV(dataframe, NewFileName): 
    path ='/Users/mmacgilvray18/Desktop/' 
    dataframe.to_csv (path+NewFileName,sep='\t')



''' Add the SI common name to the file '''
Annotation_File_DF=pd.read_csv('/Users/mmacgilvray18/Desktop/Annotation file_dashes_remain_No_duplicate_Common_names.csv') 

Merged_left_Again = pd.merge(left=Merged_left_FINAL,right=Annotation_File_DF, how='left', left_on='Shared_Interactor', right_on='Protein_Name') # complete the merge toget the common names
Merged_left_Again=Merged_left_Again[['SI_Module', 'Shared_Interactor', 'Common_Name', 'Motif_Containing_Proteins', 'subModule_Name']] # retain only these columns
Merged_left_Again.columns=['SI_Module', 'Shared_Interactor', 'SI_Name', 'Motif_Containing_Proteins', 'subModule_Name'] # rename columns 

Merged_left_Again['Protein1:Protein2']=Merged_left_Again['Shared_Interactor'].map(str) + ":" + Merged_left_Again['Motif_Containing_Proteins'] # adding column for merge section below.

Merged_left_Again['Protein1:Protein2'] = Merged_left_Again['Protein1:Protein2'].str.replace('-', '') # Removing the dashes from the names because if they remain in this column, the merge below will fail, because the background network lacks dashes in the gene annotations. 

###########################################################################################


Background_network_correct_Orientation=pd.read_csv('/Users/mmacgilvray18/Desktop/phospho_v4_bgnet_siflike_withdirections_fix_Matt_Modifications_ForPipeline.csv') # import the salt background network with the correct orientations (this file lacks dashes in gene annotations!)

Merged_left_Again_Get_Correct_Protein_Orientiations=pd.merge(left=Merged_left_Again, right=Background_network_correct_Orientation, how='left', left_on='Protein1:Protein2', right_on='Protein1:Protein2') # merge based on the columns to the left

Merged_left_Again_Get_Correct_Protein_Orientiations=Merged_left_Again_Get_Correct_Protein_Orientiations[['SI_Module', 'Shared_Interactor', 'SI_Name', 'Motif_Containing_Proteins', 'subModule_Name','Interaction']] # retain only these columns


Merged_left_Again_Get_Correct_Protein_Orientiations.columns=['SI_Module', 'Shared_Interactor', 'SI_name', 'Motif_Containing_Proteins', 'subModule_Name','Interaction_Directionality'] # rename columns 
print (Merged_left_Again_Get_Correct_Protein_Orientiations.head(2))

DF_to_CSV(Merged_left_Again_Get_Correct_Protein_Orientiations, 'DTT_T120_Prep_for_Orientation_Script_Sept2017.csv')
