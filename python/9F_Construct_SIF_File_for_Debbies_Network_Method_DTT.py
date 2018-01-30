import pandas as pd

''' The function of this script is to produce a SIF file that Debbie can use to Infer a signaling network using an ILP progamming method. The output SIF file can also be opended with Cytoscape to view a signaling network that has NOT been inferred.  



***Overview/Important Notes***
-The final SIF file indicates directionality between interactions or information flow between entitites. For example, the protein in column "A" acts upon the protein/subModule in column "B". The submodule in column "A" contains the proteins listed in column "B"

-There are 6 Interaction types in this file, listed below:

A) Motif-matched: This is a Kinase-SI that recognizes the phosphorylation motif for a subModule at a predetermined FDR cutoff. This interaction type only exists for 
kinases that are Input for a subModule (and thus can potentially regulate them). If a kinase is a match to a subModule, but is an output, it is unlikely to regulate that subModule, since
it is downstream of the subModule, and would be considered an Output (see below).

B) Unknown-recognition-motif: A Kinase or Phoshphatase SI for which we have no information about the phosphorylation motif it recognizes (ie, not in the Mok et al Dataset)
This interaction type is input only (Kinase/Phosphatase SI -> subModule)
All of our SI-Phosphatase inputs fall into this group, since we do not know their recognized phosphorylation motif.
 
C) Motif-unmatched: A Kinase SI which did NOT meet our FDR cutoff for a subModule. This is also only for Input Kinases.

D) Output: subModule -> Kinase/Phosphatase SI. Key here is that all subModule-SI interactions face TOWARDS the SI, and thus the SI is an output and unlikely to regulate the submodule.

E) Constituent: A subModule to it's protein constituents (proteins that are part of the subModule). For this group, many of the constituent proteins will NOT be SIs.
   
E) Shared_Interaction: Non-Kinase/Phosphatase SI connected to it's subModule. Can be either SI -> subModule or subModule -> SI (so an Input/Output based on interaction directionality)


Files that will always be used to create the SIF file:

-Annotation_kinases.csv = this file contains all kinases, including which kinases are in the Mok Dataset.
-kinase_phosphatase_yeast.csv = this file contains all kinases and phosphatases in yeast (annotated as kinase/phosphatase catalytic-from Mike Tyers Kinome project)
-List of enriched SIs and their subModules
-List of subModules and their protein constituents 
-KL scoring system for kinases to their subModules (which is actually based on comparing Kinases to their Modules).

Example input files:

All_SI_subModule_relationships_DF= pd.read_csv('/Users/mmacgilvray18/Desktop/Input/SIs_for_Making_SIF_File.csv')
subModule_Constituent_Proteins_DF= pd.read_csv('/Users/mmacgilvray18/Desktop/Input/subModules_with_and_without_Phenotypes_hog1_pde2_0.4_0.6_cdc14_0.2_0.4_Cutoffs.csv')
Annotation_Kinases_DF=pd.read_csv('/Users/mmacgilvray18/Desktop/Input/Annotation_kinases.csv') 
Kinases_Phosphatases_yeast_DF=pd.read_csv('/Users/mmacgilvray18/Desktop/Input/kinase_phosphatase_yeast.csv')
KL_Matching_Kinases_Modules_DF=pd.read_csv('/Users/mmacgilvray18/Desktop/Input/All_Mok_Kinases_Shuffled1000x_Newest_Method_Compared_17Modules_FDR_Scores_For_Making_SIF_file.csv')

'''
All_SI_subModule_relationships_DF= pd.read_csv('SIs_subModule_Relationships_Defined_DTT_T120_Network_Input_for_making_SIF_Sept2017.csv') 


subModule_Constituent_Proteins_DF= pd.read_csv('/DTT_submodule_constituents_Sept2017.csv')

#DF contains all kinases, including the 3 Pho85-Co-activator varieties, and includes whether a kinase is in the Mok dataset or not.
Annotation_Kinases_DF=pd.read_csv('Annotation_kinases_Updated_Correct.csv') 

#DF contains all protein annotated as kinase catalytic or phosphatase catalytic from Mike Tyers Kinome project_base
Kinases_Phosphatases_yeast_DF=pd.read_csv('kinase_phosphatase_yeast.csv')


KL_Matching_Kinases_Modules_DF=pd.read_csv('MotifMatch_Scores_for_FASTA.csv')


#-------------------------------------------------------------------------------------------------------------------------------
def DF_to_CSV(dataframe, NewFileName): 
    dataframe.to_csv (NewFileName,sep='\t')  

#-------------------------------------------------------------------------------------------------------------------------------
#Add a column to the SI_File
All_SI_subModule_relationships_DF['SI']='Yes'

#-------------------------------------------------------------------------------------------------------------------------------
#Generating Interaction type: Output (subModule -> subModule constituent proteins)

subModule_Constituent_Proteins_DF=subModule_Constituent_Proteins_DF[['Protein', 'subModule']]  # Only retain the listed columns in the dataframe
subModule_Constituent_Proteins_DF['Protein_Constituent_subModule']=subModule_Constituent_Proteins_DF.Protein.map(str) + "_" + subModule_Constituent_Proteins_DF.subModule  # create a new column that merges the protein name and subModule name together
subModule_Constituent_Proteins_DF_dupes_removed=subModule_Constituent_Proteins_DF.drop_duplicates(subset='Protein_Constituent_subModule') # drop duplicates, so only a single occurrence is listed for each SI-subModule
subModule_Constituent_Proteins_DF_dupes_removed_reorganized=subModule_Constituent_Proteins_DF_dupes_removed.rename(columns={'subModule':'Interactor_A', 'Protein':'Interactor_B'}) # Change column names
subModule_Constituent_Proteins_DF_dupes_removed_reorganized['Edge_Type']='Constituent'  # Add the Edge type (Interaction type-column). This line has been modified since the original script was created.
subModule_Constituent_Proteins_DF_dupes_removed_reorganized=subModule_Constituent_Proteins_DF_dupes_removed_reorganized[['Interactor_A','Edge_Type','Interactor_B','Protein_Constituent_subModule']] # reorganize columns
subModule_Constituent_Proteins_merged_left_DF=pd.merge(left=subModule_Constituent_Proteins_DF_dupes_removed_reorganized,right=Kinases_Phosphatases_yeast_DF, how='left', left_on='Interactor_B', right_on='ORF')  # do a merge, so I get information about which proteins are kinases/phosphatases
subModule_Constituent_Proteins_merged_left_DF=subModule_Constituent_Proteins_merged_left_DF[['Interactor_A', 'Edge_Type', 'Interactor_B', 'Annotation', 'Protein_Constituent_subModule']]  # drop the columns I don't want 
subModule_Constituent_Proteins_merged_left_DF_Add_SI_Information=pd.merge(left=subModule_Constituent_Proteins_merged_left_DF, right=All_SI_subModule_relationships_DF, how='left', left_on='Protein_Constituent_subModule', right_on='SI_Module')
subModule_Constituent_Proteins_merged_left_DF_Add_SI_Information=subModule_Constituent_Proteins_merged_left_DF_Add_SI_Information[['Interactor_A', 'Edge_Type', 'Interactor_B', 'Annotation', 'SI']]

#-------------------------------------------------------------------------------------------------------------------------------
# Split the SI_subModule_relationships_DF by Input/Output

#Splitting the Dataframe based on if the SI is acting as an Input or Output
#If a SI_subModule relationship is an output, then it should be listed that Interactor_A is the subModule and Interactor_B is the SI. 

SIs_subModules_Output=All_SI_subModule_relationships_DF.loc[All_SI_subModule_relationships_DF['Shared_Interactor_subModule_Relationship']=='Output'] # All interactions that go from subModule -> SI
SIs_subModules_Input=All_SI_subModule_relationships_DF.loc[All_SI_subModule_relationships_DF['Shared_Interactor_subModule_Relationship']=='Input']

#---------------------------------------------------------------------------------------------------------------------------------
# Generating Interaction type: Output (subModule -> SIs-Kinase (that have all Interactions facing from subModule to SI)    

SIs_subModules_Output_renamed=SIs_subModules_Output.rename(columns={'subModule_Name':'Interactor_A','Shared_Interactor':'Interactor_B'}) #rename columns appropriately. 
#print (SIs_subModules_Output_renamed.head(5))
SIs_subModules_Output_renamed['Edge_Type']='Output'  # Add the edge_type, which is output in this case 

SIs_subModules_Output_renamed=SIs_subModules_Output_renamed[['Interactor_A', 'Edge_Type', 'Interactor_B', 'SI_Module']] # drop unwanted columns
SIs_subModules_Output_renamed_merge_left=pd.merge(left=SIs_subModules_Output_renamed,right=Kinases_Phosphatases_yeast_DF, how='left', left_on='Interactor_B', right_on='ORF')  # do a merge, so I get information about which proteins are kinases/phosphatases
SIs_subModules_Output_renamed_merge_left=SIs_subModules_Output_renamed_merge_left[['Interactor_A', 'Edge_Type', 'Interactor_B', 'Annotation', 'SI_Module']] #Drop unwanted columns from the above merge 

SIs_subModules_Output_renamed_merge_left_left_again=pd.merge(left=SIs_subModules_Output_renamed_merge_left, right=All_SI_subModule_relationships_DF, how='left', left_on='SI_Module', right_on='SI_Module') # merge so we get SI information 
SIs_subModules_Output_renamed_merge_left_left_again=SIs_subModules_Output_renamed_merge_left_left_again[['Interactor_A', 'Edge_Type', 'Interactor_B', 'Annotation', 'SI']] # drop columns we don't want
SIs_subModules_Output_renamed_merge_left_left_again_Kinase_SIs_Only=SIs_subModules_Output_renamed_merge_left_left_again.loc[SIs_subModules_Output_renamed_merge_left_left_again['Annotation']=='Kinase']
SIs_subModules_Output_renamed_merge_left_left_again_Kinase_SIs_Only=SIs_subModules_Output_renamed_merge_left_left_again_Kinase_SIs_Only.loc[SIs_subModules_Output_renamed_merge_left_left_again_Kinase_SIs_Only['Annotation']=='Kinase'] # only keep annotations that are Kinases!


#---------------------------------------------------------------------------------------------------------------------------------
# Generating Interaction type: Output (subModule -> SIs_Phosphatase (that have all interactions facing from submodule to SI)
SIs_subModules_Output_renamed_merge_left_left_Twice=pd.merge(left=SIs_subModules_Output_renamed_merge_left, right=All_SI_subModule_relationships_DF, how='left', left_on='SI_Module', right_on='SI_Module')
SIs_subModules_Output_renamed_merge_left_left_Twice=SIs_subModules_Output_renamed_merge_left_left_Twice[['Interactor_A', 'Edge_Type', 'Interactor_B', 'Annotation', 'SI']] # drop columns we don't want
SIs_subModules_Output_renamed_merge_left_left_Twice_Phosphatase_SIs_Only = SIs_subModules_Output_renamed_merge_left_left_Twice.loc[SIs_subModules_Output_renamed_merge_left_left_Twice['Annotation']=='Phosphatase']


#---------------------------------------------------------------------------------------------------------------------------------
# Generating Interaction type: Shared_Interaction: Non-Kinase/Phosphatase Shared Interactors and their Module associations. can be directed as follows: subModule -> non-kinase/phosphatase SI  or non-kinase/phosphatase SI -> subModule
# Note: These are All Shared Interactors, so we can simply add a SI column to this file. 

# This section of code is making the subModule -> SI direction 
SIs_subModules_Output_renamed_merge_left=SIs_subModules_Output_renamed_merge_left[['Interactor_A', 'Edge_Type', 'Interactor_B', 'Annotation', 'SI_Module']] 
SIs_subModules_Output_renamed_merge_left_non_Kin_Phos=SIs_subModules_Output_renamed_merge_left.loc[SIs_subModules_Output_renamed_merge_left['Annotation']!='Kinase']
SIs_subModules_Output_renamed_merge_left_non_Kin_Phos=SIs_subModules_Output_renamed_merge_left_non_Kin_Phos.loc[SIs_subModules_Output_renamed_merge_left_non_Kin_Phos['Annotation']!='Phosphatase'] # Return all proteins with annotations that are NOT kinases!
SIs_subModules_Output_renamed_merge_left_non_Kin_Phos['SI']='SI' # Since all of these proteins are SIs, add a column that indicates they are SIs
SIs_subModules_Output_renamed_merge_left_non_Kin_Phos=SIs_subModules_Output_renamed_merge_left_non_Kin_Phos[['Interactor_A', 'Edge_Type', 'Interactor_B', 'Annotation', 'SI']]
SIs_subModules_Output_renamed_merge_left_non_Kin_Phos['Edge_Type']='Shared_Interaction'
#print (SIs_subModules_Output_renamed_merge_left_non_Kin_Phos)

#This section of code is making the SI -> subModule direction 

#print (len(SIs_subModules_Input))
SIs_subModules_Input_merge_left_left_again=pd.merge(left=SIs_subModules_Input,right=Kinases_Phosphatases_yeast_DF, how='left', left_on='Shared_Interactor', right_on='ORF') # Merge so we get Kinae/phosphatase information for SIs
SIs_subModules_Input_merge_left_left_again_non_Kin_Phos=SIs_subModules_Input_merge_left_left_again.loc[SIs_subModules_Input_merge_left_left_again['Annotation']!='Kinase']  # Drop kinases
SIs_subModules_Input_merge_left_left_again_non_Kin_Phos=SIs_subModules_Input_merge_left_left_again_non_Kin_Phos.loc[SIs_subModules_Input_merge_left_left_again_non_Kin_Phos['Annotation']!='Phosphatase'] # drop phosphatases 
SIs_subModules_Input_merge_left_left_again_non_Kin_Phos['Edge_Type']='Shared_Interaction' # add the edge type 

SIs_subModules_Input_merge_left_left_again_non_Kin_Phos=SIs_subModules_Input_merge_left_left_again_non_Kin_Phos[['Shared_Interactor', 'Edge_Type', 'subModule_Name', 'Annotation', 'SI']] # drop unwanted columns
SIs_subModules_Input_merge_left_left_again_non_Kin_Phos_renamed=SIs_subModules_Input_merge_left_left_again_non_Kin_Phos.rename(columns={'Shared_Interactor':'Interactor_A', 'subModule_Name':'Interactor_B'}) # rename columns so I can merge all dataframes later on.


#---------------------------------------------------------------------------------------------------------------------------------
# Generating Interaction Type: Motif-Matched  (SI-Kinase -> subModule)

SIs_subModules_Input_merge_left_left_again_Kinase_Only=SIs_subModules_Input_merge_left_left_again.loc[SIs_subModules_Input_merge_left_left_again['Annotation']== 'Kinase'] # subset the dataframe so we are only working with Kinases.['a'] = df['a'].apply(lambda x: x.split('-')[0])
SIs_subModules_Input_merge_left_left_again_Kinase_Only['SI_V2']=SIs_subModules_Input_merge_left_left_again_Kinase_Only['SI_Module'].apply(lambda x: x.split('_')[0])  # Get term before the first "_"
SIs_subModules_Input_merge_left_left_again_Kinase_Only['Cluster']=SIs_subModules_Input_merge_left_left_again_Kinase_Only['SI_Module'].apply(lambda x: x.split('_')[1]) # Get term after the first "_"
SIs_subModules_Input_merge_left_left_again_Kinase_Only['motif']=SIs_subModules_Input_merge_left_left_again_Kinase_Only['SI_Module'].apply(lambda x: x.split('_')[2]) # Get term after the second "_"
SIs_subModules_Input_merge_left_left_again_Kinase_Only['SI_MODULE']=SIs_subModules_Input_merge_left_left_again_Kinase_Only.SI_V2.map(str) + "_" + SIs_subModules_Input_merge_left_left_again_Kinase_Only.Cluster + "_" + SIs_subModules_Input_merge_left_left_again_Kinase_Only.motif
DF_to_CSV(SIs_subModules_Input_merge_left_left_again_Kinase_Only, "Test1.csv")
SIs_subModules_Input_merge_left_left_again_Kinase_Only_merge_left_KL=pd.merge(left=SIs_subModules_Input_merge_left_left_again_Kinase_Only, right=KL_Matching_Kinases_Modules_DF, how='left', left_on='SI_MODULE', right_on='Kinase_Module') # Perform a merge so I get information about matching Kinases to subModule motifs (KL script output).
#print (SIs_subModules_Input_merge_left_left_again_Kinase_Only_merge_left_KL)
#DF_to_CSV(SIs_subModules_Input_merge_left_left_again_Kinase_Only_merge_left_KL, 'Test2.csv')
SIs_subModules_Input_merge_left_left_again_Kinase_Only_merge_left_KL=SIs_subModules_Input_merge_left_left_again_Kinase_Only_merge_left_KL[['Shared_Interactor', 'motif_match', 'subModule_Name', 'Annotation', 'SI', 'FDR']] # drop columns I don't want in the final version
SIs_subModules_Input_merge_left_left_again_Kinase_Only_merge_left_KL_renamed=SIs_subModules_Input_merge_left_left_again_Kinase_Only_merge_left_KL.rename(columns={'Shared_Interactor':'Interactor_A', 'motif_match':'Edge_Type', 'subModule_Name':'Interactor_B'}) # rename columns so I can merge all dataframes in the future.
SIs_subModules_Input_merge_left_left_again_Kinase_Only_merge_left_KL_renamed['Edge_Type']=SIs_subModules_Input_merge_left_left_again_Kinase_Only_merge_left_KL_renamed['Edge_Type'].map({'yes':'motif_match', 'no':'motif_unmatched', 'unknown_recognition_motif':'unknown_recognition_motif'})
#print (SIs_subModules_Input_merge_left_left_again_Kinase_Only_merge_left_KL_renamed)
#DF_to_CSV(SIs_subModules_Input_merge_left_left_again_Kinase_Only_merge_left_KL_renamed, 'New_Script_file_for_merge.csv')
#print (KL_Matching_Kinases_Modules_DF.head(5))

#---------------------------------------------------------------------------------------------------------------------------------
# Generating Interaction Type: unknown-recognition motif  (SI-Phosphatase -> subModule)

SIs_subModules_Input_merge_left_left_again_Phosphatase_Only_unknown_recogntion_motif=SIs_subModules_Input_merge_left_left_again.loc[SIs_subModules_Input_merge_left_left_again['Annotation']== 'Phosphatase'] # subset the dataframe so we are only working with Phosphatases.
SIs_subModules_Input_merge_left_left_again_Phosphatase_Only_unknown_recogntion_motif['Edge_Type']='unknown_recognition_motif'  # add the edge type 
SIs_subModules_Input_merge_left_left_again_Phosphatase_Only_unknown_recogntion_motif=SIs_subModules_Input_merge_left_left_again_Phosphatase_Only_unknown_recogntion_motif[['Shared_Interactor', 'Edge_Type','subModule_Name', 'Annotation', 'SI']]
SIs_subModules_Input_merge_left_left_again_Phosphatase_Only_unknown_recogntion_motif_renamed=SIs_subModules_Input_merge_left_left_again_Phosphatase_Only_unknown_recogntion_motif.rename(columns={'Shared_Interactor':'Interactor_A', 'subModule_Name':'Interactor_B'})
#print (SIs_subModules_Input_merge_left_left_again_Phosphatase_Only_unknown_recogntion_motif_renamed)

#---------------------------------------------------------------------------------------------------------------------------------
#Adding empty columns to a few of the dataframes so that the final append will work (requires that all files have the same column names,otherwise you'll end up with Column_X, Column_Y)

subModule_Constituent_Proteins_merged_left_DF_Add_SI_Information['FDR_Score']=""
subModule_Constituent_Proteins_merged_left_DF_Add_SI_Information['Match_FDR']=""
#print (subModule_Constituent_Proteins_merged_left_DF_Add_SI_Information.head(5))


SIs_subModules_Output_renamed_merge_left_left_again_Kinase_SIs_Only['FDR_Score']=""
#print (SIs_subModules_Output_renamed_merge_left_left_again_Kinase_SIs_Only.head(10))
SIs_subModules_Output_renamed_merge_left_left_again_Kinase_SIs_Only['Match_FDR']=""
#print (SIs_subModules_Output_renamed_merge_left_left_again_Kinase_SIs_Only.head(2))

SIs_subModules_Output_renamed_merge_left_non_Kin_Phos['FDR_Score']=""
SIs_subModules_Output_renamed_merge_left_non_Kin_Phos['Match_FDR']=""
#print (SIs_subModules_Output_renamed_merge_left_non_Kin_Phos.head(2))

SIs_subModules_Input_merge_left_left_again_non_Kin_Phos_renamed['FDR_Score']=""
SIs_subModules_Input_merge_left_left_again_non_Kin_Phos_renamed['Match_FDR']=""
#print (SIs_subModules_Input_merge_left_left_again_non_Kin_Phos_renamed.head(4))


SIs_subModules_Input_merge_left_left_again_Kinase_Only_merge_left_KL_renamed=SIs_subModules_Input_merge_left_left_again_Kinase_Only_merge_left_KL_renamed.rename(columns={'FDR':'FDR_Score'})
SIs_subModules_Input_merge_left_left_again_Kinase_Only_merge_left_KL_renamed['Match_FDR']=""
#DF_to_CSV(SIs_subModules_Input_merge_left_left_again_Kinase_Only_merge_left_KL_renamed, 'Test101.csv')
#print (SIs_subModules_Input_merge_left_left_again_Kinase_Only_merge_left_KL_renamed.head(4))

SIs_subModules_Input_merge_left_left_again_Phosphatase_Only_unknown_recogntion_motif_renamed['FDR_Score']=""
SIs_subModules_Input_merge_left_left_again_Phosphatase_Only_unknown_recogntion_motif_renamed['Match_FDR']=""
#print (SIs_subModules_Input_merge_left_left_again_Phosphatase_Only_unknown_recogntion_motif_renamed.head(4))

SIs_subModules_Output_renamed_merge_left_left_Twice_Phosphatase_SIs_Only['FDR_Score']=""
SIs_subModules_Output_renamed_merge_left_left_Twice_Phosphatase_SIs_Only['Match_FDR']=""

FinalDF=subModule_Constituent_Proteins_merged_left_DF_Add_SI_Information.append(SIs_subModules_Output_renamed_merge_left_left_again_Kinase_SIs_Only)
#DF_to_CSV(FinalDF, "Test.csv")
FinalDF_2=FinalDF.append(SIs_subModules_Output_renamed_merge_left_non_Kin_Phos)
FinalDF_3=FinalDF_2.append(SIs_subModules_Input_merge_left_left_again_non_Kin_Phos_renamed)
FinalDF_4=FinalDF_3.append(SIs_subModules_Input_merge_left_left_again_Kinase_Only_merge_left_KL_renamed)
FinalDF_5=FinalDF_4.append(SIs_subModules_Input_merge_left_left_again_Phosphatase_Only_unknown_recogntion_motif_renamed)
FinalDF_6=FinalDF_5.append(SIs_subModules_Output_renamed_merge_left_left_Twice_Phosphatase_SIs_Only)
DF_to_CSV(FinalDF_6, 'DTT_T120_SIF_FINAL_Sept2017.csv')
