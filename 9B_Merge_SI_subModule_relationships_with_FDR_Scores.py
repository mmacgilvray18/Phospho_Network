import pandas as pd 


''' This script takes the output of the KLD FDR script and adds these FDR scores to the identified SI-submodule pairs. 
'''

FDR_Scores_DF=pd.read_csv('/Users/mmacgilvray18/Desktop/FDR_Input_T120_ForMerge_withSIs.csv')  # All FDR scores for Mok kinase and module PWM comparison
SIs_DF=pd.read_csv('/Users/mmacgilvray18/Desktop/DTT_T120_SIs_All_Sept2017.csv')  # SIs - All, enriched and not enriched
Kinase_Names_DF=pd.read_csv('/Users/mmacgilvray18/Desktop/Kinases_Mok_andNOT_In_Mok.csv') # Contains all kinases in Mok and not in Mok. Contains common names and YORFs, also contains modifications for Pho85 naming (ex: Pho85-Pcl is now Pho85 in one column)

def DF_to_CSV(dataframe, NewFileName): 
    path ='/Users/mmacgilvray18/Desktop/' 
    dataframe.to_csv (path+NewFileName,sep='\t') 
####################################################################################################################################
''' Merging the FDR_Scores_DF with the Kinase_Names_DF '''
# This step is completed so that the correct Pho85 nomenclature is used for a subsquent merge with the SI dataframe. This is necessary because there are 3 pho85-cofactor variants in the Mok et al, dataset.

FDR_Scores_DF_merged_left = pd.merge(left=FDR_Scores_DF,right=Kinase_Names_DF, how='left', left_on='Kinase', right_on='Kinase') # completing a merge so that all the kinase nomenclature from the Kinase_Names_DF

####################################################################################################################################
''' Merge the Kinase name (common name,ex. Hog1) with the Module name to create a new column called "Candidate_Kinase_Regulators" '''

FDR_Scores_DF_merged_left['Candidate_Kinase_Regulators'] = FDR_Scores_DF_merged_left.Module.map(str) + "_" + FDR_Scores_DF_merged_left.Kinase_Pho85_renamed # creating a new column that is the result of a merge between Kinase_Pho85_renamed, and Module
#print (FDR_Scores_DF_merged_left.head(1))
DF_to_CSV(FDR_Scores_DF_merged_left, 'yay3.csv')
####################################################################################################################################
''' Filtering out non-enriched SIs'''

SIs_Filtered=SIs_DF.loc[SIs_DF['BH_significant'] == 1] # Filter out non-significant shared interactors (anything with a "0") 

####################################################################################################################################
''' Adding to the SIs file, the "common" name for the proteins, rather than just using the YORF designation'''
SIs_Filtered_merged_left= pd.merge(left=SIs_Filtered,right=Kinase_Names_DF, how='left', left_on='Shared_Interactor', right_on='Kinase_YORF')


####################################################################################################################################
''' Splitting 'Motif' column and producing a new column, called "Module" that only lists the Induced/Repressed WT phenotype and the motif '''

def Split_After_2nd_Occurence_In_A_String_Retaining_Beginning():
    lst=[] # create an empty list 
    for string in SIs_Filtered_merged_left['Motif']: # Select the string from the "Motif" column 
        strip_character ="_"  # define character where strip will occur
        lst.append(strip_character.join(string.split(strip_character)[:2])) # append to the list the text before the second occurence of the character "_"
    Series_Object = pd.Series(lst) # put the list into a series 
    SIs_Filtered_merged_left['Module'] = Series_Object.values # append the series values to the already existing DF in a new column
    return SIs_Filtered_merged_left
        
SIs_Filtered_merged_left_String_Split=Split_After_2nd_Occurence_In_A_String_Retaining_Beginning()
#print (SIs_Filtered_merged_left_String_Split)
####################################################################################################################################    
####################################################################################################################################    
'''Creating New Columns that can be used for a merge'''
SIs_Filtered_merged_left_String_Split['Kinase_subModules'] = SIs_Filtered_merged_left_String_Split.Motif.map(str) + "_" + SIs_Filtered_merged_left_String_Split.Kinase_Pho85_renamed


SIs_Filtered_merged_left_String_Split['Kinase_Modules'] = SIs_Filtered_merged_left_String_Split.Module.map(str) + "_" + SIs_Filtered_merged_left_String_Split.Kinase_Pho85_renamed

#################################################################################################################################### 
''' Perform a merge where of the FDR Scores Dataframe with the SIs_Filtered_merged_left_String_Split DF. 
This will reveal if a kinase-Module relationship, from the FDR Score Dataframe, which contains all possible Kinase-Module relationships, exist in the users 
kinase-subModule file (so the SIs file)'''
    
merged_left= pd.merge(left=FDR_Scores_DF_merged_left,right=SIs_Filtered_merged_left_String_Split, how='left', left_on='Candidate_Kinase_Regulators', right_on='Kinase_Modules')


#################################################################################################################################### 
''' Drop columns that are not needed or redundant '''

merged_left=merged_left[['Scores', 'Kinase_x', 'Module_x', 'Candidate_Kinase_Regulators', 'Counts_Less_Than', 'Number_of_Scores', 'FDR', 'Kinase_subModules']]


#################################################################################################################################### 
''' Drop NaN values '''
merged_left=merged_left.dropna(subset=['Kinase_subModules']) # Drop the NaN values, so that the dataframe only contains Kinases that were connected to subModules.


####################################################################################################################################
''' Drop any duplicates that occur in TWO columns - this is done only because of Pho85 being listed 3 times (because of it's co-factor interactions) and that affects the merge'''
merged_left=merged_left.drop_duplicates(subset=['Kinase_x', 'Kinase_subModules']) # only drop duplicates that are found in BOTH columns

DF_to_CSV(merged_left, 'All_DTT_T120_Kinase_Module_FDR_Scores_and_their_Kinase_SI_subModules_Sept2017.csv')
