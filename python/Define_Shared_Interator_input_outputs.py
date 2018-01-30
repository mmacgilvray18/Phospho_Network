''' 
The function of this script is as follows: For each SI and it's interactions a submodules constituents, the script determine if the 
SI is a likely submodule regulator (that is, the Shared Interactor has at least 1 directional interaction, or ppi interaction, with a subModule protein aimed from the SI to the submodule), or if the subModule proteins are act upon the SI (that is, all interactions between the SI and subModule proteins have the 'Reverse' designation', indicating that the subModule proteins act upon the SI).  

-If all of the interactions are reversed, then the script will define the relationship between the SI and the subModule as "Output", indicating that the SI is likely downstream of the submodule and is not likely regulating the submodule.
-If there is at least one interaction that is NOT reverse (ie, kinase-substrate) or is and non-directed ppi, the relationship between the SI and the subMOdule is defined as "Input", suggesting there is a possibility that the SI can regulate the submodule protein phosphorylation state.


This script takes an input file that contains the following:
- All enriched Shared Interactors (SIs) (according to HyperG) and their connections to subModules.
- All known protein interactions for each SI (ppi, kinase-substrate, etc)
- Many of these interactions are directed (kinase-substrate, metabolic pathway, etc). PPI are not a directed interaction.

'''
import pandas as pd

Input_df=pd.read_csv('/Users/mmacgilvray18/Desktop/DTT_T120_Prep_for_Orientation_Script_Sept2017.csv')

# Split the input DF into independent DFs based on the term in the SI_Module column (this columns contains the SI and it's connection to each subModule). 
def Split_based_on_SI_Module_Column():
    DF_lst =[]
    for SI_Module in Input_df['SI_Module'].unique():
        DF=Input_df.loc[Input_df['SI_Module']==SI_Module]
        DF_lst.append(DF)
    return DF_lst

DF_lst=Split_based_on_SI_Module_Column()

# This function counts, for each DF, how many of the interactions are reversed. It also counts the length of the dataframe, and then
# subtracts the the length of the dataframe from the counts. If the resultant value is 0, then all of the interactions were reversed.
def Count_Instances_of_Reverse_Interaction():
    DF_Counts_lst=[]
    for df in DF_lst:
        df=df.copy()
        df['Counts']=df.Interaction_Directionality.str.contains('Reversed').sum()  # Count the number of interactions that are "Reversed"
        x=len(df)
        df['Length']=x
        df['Counts_Length']=df['Counts']-df['Length']
        
        DF_Counts_lst.append(df)
    return DF_Counts_lst

DF_Counts_lst=Count_Instances_of_Reverse_Interaction()
print (DF_Counts_lst)

# This function assigns 'Input' and 'Output' classifications based on the 'Counts_Length' column in the dataframe. 
def Only_Reverse_Interactions_Move_to_Outgoing_Columns():
    df_Modified_Outgoing_lst=[]
    for df in DF_Counts_lst:
        #print (df.dtypes)
        for value in df['Counts_Length'].unique():
            #print (value)
            if value == 0:
                df['Shared_Interactor_subModule_Relationship']= 'Output'
                df_Modified_Outgoing_lst.append(df)
            else:
                df['Shared_Interactor_subModule_Relationship']= 'Input'
                df_Modified_Outgoing_lst.append(df)
                
    return df_Modified_Outgoing_lst
    
df_Modified_Outgoing_lst=Only_Reverse_Interactions_Move_to_Outgoing_Columns()

# This function concatenates the dataframes back together, leaving a single DF. 
def ConcatenateDFs():   
    EmptyDF = pd.DataFrame() # create an empty dataframe
    for df in df_Modified_Outgoing_lst:  # select a dataframe in the list 
        df=df.copy() # make a copy of that dataframe 
        EmptyDF=EmptyDF.append(df) # append to the empty DF the dataframe selected and overwrite the empty dataframe
    return EmptyDF

Final=ConcatenateDFs()

print (Final)

#Function writes out a Dataframe to a CSV file. 
def DF_to_CSV(dataframe, NewFileName): 
    path ='/Users/mmacgilvray18/Desktop/' 
    dataframe.to_csv (path+NewFileName,sep='\t') 
    

Final_Keep_Columns_Needed_For_SIF=Final[['SI_Module', 'Shared_Interactor', 'subModule_Name', 'Shared_Interactor_subModule_Relationship']] 
Final_Keep_Columns_Needed_For_SIF=Final_Keep_Columns_Needed_For_SIF.drop_duplicates('SI_Module')


DF_to_CSV(Final_Keep_Columns_Needed_For_SIF, 'SIs_subModule_Relationships_Defined_DTT_T120_Network_Input_for_making_SIF_Sept2017.csv')
 
