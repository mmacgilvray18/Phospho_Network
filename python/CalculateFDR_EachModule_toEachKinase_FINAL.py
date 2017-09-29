import pandas as pd 
import glob
import itertools



''' This script is calculates an FDR for each Mok Kinase to each Module. 
The script takes the 63,000 shuffled Mok et al kinase-module scores and determines for each kinase where that unshuffled kinase-module score
falls in the shuffled distribution.For example,if a Kinase has an shuffled score of 14.7 to a module , and only 63 unshuffled kinase-module scores
are below that value, then this kinase has has an FDR of 0.1% (63/63,000 scores). We can then use the FDR values for all kinases to a module
to determine an FDR cutoff for that module. Thus, we can say only these x kinases are a good match to the module. Calling an FDR threshold is 
done manually by the user.
'''



#Read in input files (which are the non-shuffled scores for all Mok kinases compared to each Module) 

path=r"/Users/mmacgilvray18/Desktop/ClassA_NoShuffle_KL_Unwanted_Characters_Removed/"
filenames = glob.glob(path + "*.csv")


dfs_input_lst = []
for filename in filenames:
    dfs_input_lst.append(pd.read_csv(filename, sep=","))


def Concatenate_Input_DFs(): 
    ''' Append the dataframes to one another '''
    EmptyDF = pd.DataFrame() 
    for df in dfs_input_lst:  
        df=df.copy() 
        EmptyDF=EmptyDF.append(df) 

Input=Concatenate_Input_DFs()
Input=Input[['Scores','Kinase','Module']]





def SplitInput_df_byModule():
    ''' Function splits the input dataframe, by Module, into independent dataframes'''
    DF_Input_lst =[]
    for Module in Input['Module'].unique():
        DF=Input.loc[Input['Module']==Module]
        DF_Input_lst.append(DF)
    return DF_Input_lst

DF_Input_lst=SplitInput_df_byModule()



# Importing the Shuffled Mok kinase -Module Score csv files individually and creating dataframes
path=r"/Users/mmacgilvray18/Desktop/Module_PWMs_KLD_Output_ClassA_Unwanted_Characters_Removed/"
filenames = glob.glob(path + "*.csv")
#print (filenames)

dfs_lst = []
for filename in filenames:
    dfs_lst.append(pd.read_csv(filename, sep=","))

dfs_lst2=[]
for df in dfs_lst:
    df=df[['Scores','Kinase', 'Module']]
    dfs_lst2.append(df)




def CountScores_Below():
    '''Function is calculating the number of scores in the shuffled distribution below a non-shuffled Kinase_Module score '''
    Final_DFs_lst=[]   #
    for DF_shuffled in dfs_lst2:  
        DF_shuffled=DF_shuffled.copy()
       
        for df_noShuffle in DF_Input_lst:  
            df_noShuffle=df_noShuffle.copy()
       
            if df_noShuffle['Module'].unique().all() == DF_shuffled['Module'].unique().all():      # if all of the values match in the module column for each df, then and only then, perform the below steps
                Input_Score_lst=df_noShuffle['Scores'].tolist()  
                
                Scores_lst=[]  
                for score in Input_Score_lst: 
                    Scores_lst_individual=[]
                    num_smaller_items = (DF_shuffled['Scores']<score).sum()                        # create a variable that is the sum of all scores below the score in the Shuffled_Scores dataframe
                 
                    Scores_lst_individual.append(num_smaller_items)                                # append the number of scores below a given kinase-module score to the individual list.
                    Scores_lst.append(Scores_lst_individual)                                       # append the individual scores to a list.
                    merged = list(itertools.chain(*Scores_lst))                                    # Flatten the list of lists. 
                    
                df_noShuffle['Counts_Less_Than']=merged
                df_noShuffle['Number_of_Scores']=len(DF_shuffled)                                  # take all of the summed scores, one per kinase from the kinase-module no-shuffle dataframe, and create a new column.
                Final_DFs_lst.append(df_noShuffle)  
    return Final_DFs_lst

DFs_with_CountsBelow_lst=CountScores_Below()

            
   
    
def Calculate_FDR():
    ''' Function calculates an FDR value by dividing the number of shuffled scores for a kinase-module
    that are smaller than the non-shuffled kinase-module score by all shuffled scores (63,000)  '''
    DFs_with_CountsBelow_lst2=[]
    for DF in DFs_with_CountsBelow_lst:
        DF['FDR']=DF['Counts_Less_Than']/DF['Number_of_Scores']
        DF=DF.sort(['FDR'], ascending=[True])
        DFs_with_CountsBelow_lst2.append(DF)
    return DFs_with_CountsBelow_lst2
 
 
DFs_with_CountsBelow_lst2=Calculate_FDR()



# Import the kinases not in the Mok et al dataset.
Kinases_Not_In_Mok_DF=pd.read_csv('/Users/mmacgilvray18/Desktop/Kinases_Not_In_Mok.csv')


def ConcatenateDFs_with_Kinases_Not_In_Mok():
    '''Function adds the kinases not in the Mok et al dataset to the dataframes for each module'''
    DFs_with_CountsBelow_lst3=[]
    for DF in DFs_with_CountsBelow_lst2:
        DF=DF.copy()
        FinalDF=DF.append(Kinases_Not_In_Mok_DF)
        DFs_with_CountsBelow_lst3.append(FinalDF)
    return DFs_with_CountsBelow_lst3

DFs_with_CountsBelow_lst3=ConcatenateDFs_with_Kinases_Not_In_Mok() 
print (DFs_with_CountsBelow_lst3)   



def ConcatenateDFs(): 
    '''Function appends all of the dataframes, for each module, together into one dataframe''' 
    EmptyDF = pd.DataFrame() 
    for df in DFs_with_CountsBelow_lst3: 
        df=df.copy() 
        EmptyDF=EmptyDF.append(df) 
    return EmptyDF

Final=ConcatenateDFs()





def DF_to_CSV(dataframe, NewFileName): 
    ''' Write out dataframe as a tab separated file.'''
    path ='/Users/mmacgilvray18/Desktop/' 
    dataframe.to_csv (path+NewFileName,sep='\t') 
    
DF_to_CSV(Final, 'All_Mok_Kinases_Shuffled1000x_New_Method_ClassA_Modules_FDR_Scores_8_30_16.csv')
