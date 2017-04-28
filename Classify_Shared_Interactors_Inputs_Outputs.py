import pandas as pd


''' 
The function of this script is as follows: For each SI and it's connections with submodule protein constituents, determine if the 
SI acts upon the submodule (that is, the Shared Interactor has at least 1 directional interaction, or ppi interaction, with a submodule protein), or if the 
submodule acts upon the SI (that is, all interactions between the SI and submodule proteins have the 'Reverse' designation', indicating that the submodule proteins act
upon the SI).  

-If all of the interactions are reversed, then the script will define the relationship between the SI and the submodule as "Output"
-If there is at least one interaction that is directed from SI towards submodule, or is a ppi, the relationship between the SI and the submodule is defined as "Input"


This script takes an input file that contains the following:
- All enriched Shared Interactors (SIs) (according to HyperG) and their connections to submodules.
- All known protein interactions for each SI (ppi, kinase-substrate, etc)
- Many of these interactions are directed (kinase-substrate, metabolic pathway, etc). PPI are not a directed interaction.

'''
Input_df=pd.read_csv('/Users/mmacgilvray18/Desktop/ClassA_Prep_for_Orientation_Script.csv')

def Split_based_on_SI_submodule_Column():
    ''' Function splits the input DF into independent DFs based on the SI-submodule column pairs. Thus, each SI and it's submodule protein interactions are
    in independent dataframes '''
    DF_lst =[]
    for SI_submodule in Input_df['SI_submodule'].unique():
        DF=Input_df.loc[Input_df['SI_submodule']==SI_submodule]
        DF_lst.append(DF)
    return DF_lst

DF_lst=Split_based_on_SI_submodule_Column()


def Count_Instances_of_Reverse_Interaction():
    ''' Function counts, for each DF, and thus each SI-submodule pair, how many of the interactions are 'reversed', or facing from submodule TOWARDS SI. 
        It also counts the length of the dataframe, and then subtracts the the length of the dataframe from the counts. If the resultant value is 0, then all of the interactions 
        were reversed '''
    DF_Counts_lst=[]
    for df in DF_lst:
        df=df.copy()
        df['Counts']=df.Interaction_Directionality.str.contains('Reversed').sum()                                                               # Count the number of interactions that are "Reversed"
        x=len(df)
        df['Length']=x
        df['Counts_Length']=df['Counts']-df['Length']
        
        DF_Counts_lst.append(df)
    return DF_Counts_lst

DF_Counts_lst=Count_Instances_of_Reverse_Interaction()

    
 
def Only_Reverse_Interactions_Move_to_Outgoing_Columns():
    '''Function assigns 'Input' and 'Output' classifications based on the 'Counts_Length' column in the dataframe. '0' values are 'outputs', all other's are 'inputs' '''
    df_Modified_Outgoing_lst=[]
    for df in DF_Counts_lst:
        for value in df['Counts_Length'].unique():
           
            if value == 0:
                df['Shared_Interactor_submodule_Relationship']= 'Output'
                df_Modified_Outgoing_lst.append(df)
            else:
                df['Shared_Interactor_submodule_Relationship']= 'Input'
                df_Modified_Outgoing_lst.append(df)
           
    return df_Modified_Outgoing_lst
    
        
df_Modified_Outgoing_lst=Only_Reverse_Interactions_Move_to_Outgoing_Columns()




def AppendDFs(): 
    '''Function appends all dataframes back together '''
    EmptyDF = pd.DataFrame()
    for df in df_Modified_Outgoing_lst: 
        df=df.copy() 
        EmptyDF=EmptyDF.append(df) 
    return EmptyDF

Final=AppendDFs()




def DF_to_TSV(dataframe, NewFileName): 
    '''Function writes out a Dataframe to a tsv file'''
    path ='/Users/mmacgilvray18/Desktop/' 
    dataframe.to_csv (path+NewFileName,sep='\t') 
    


Final_Keep_Columns_Needed_For_SIF=Final[['SI_submodule', 'Shared_Interactor', 'submodule_Name', 'Shared_Interactor_submodule_Relationship']]  
Final_Keep_Columns_Needed_For_SIF=Final_Keep_Columns_Needed_For_SIF.drop_duplicates('SI_submodule')                                                                     # Dropping duplicates entries, which are created because for each SI-submodule interaction there are numerous interactions with protein constituent. Only want a single interaction, input or output, for each SI and it's submodule. 
DF_to_TSV(Final_Keep_Columns_Needed_For_SIF, 'SIs_submodule_Relationships_Define_ClassA_Network_Input_for_making_SIF_8_30_16.csv')
 
