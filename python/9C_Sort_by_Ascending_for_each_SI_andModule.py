import pandas as pd 

''' This is a quick script that cleans up the Output of the 3_Kullback_Leibler Shuffle Script. It removes unwanted names that trail
the subModule name-these are leftovers from a previous script. It then also sorts each subModule by ascending for the SI scores

example input file: All_Mok_Kinases_Shuffled1000x_Newest_Method_Compared_17Modules_FDR_Scores_SIs_andScores.csv
example output file :All_Mok_Kinases_Shuffled1000x_Newest_Method_Compared_17Modules_FDR_Scores_SIs_andScores_Sorted_Ascending.csv '''

Input=pd.read_csv('/Users/mmacgilvray18/Desktop/All_DTT_T120_Kinase_Module_FDR_Scores_and_their_Kinase_SI_subModules_Sept2017.csv')
print (Input.dtypes)

# Remove the last occurrence of a character, and the text that follows
def Remove_Text_After_Last_Occurence_of_Character():
    Value_lst=[]
    for value in Input['Kinase_subModules']:
        value="_".join(value.split("_")[:-1]) # return everything minus the last occurrence of the "_" and what trailed
        Value_lst.append(value)
    Input['Kinase_subModules']=Value_lst
        #sep = '_'
        #value = value.split(sep, 5)[-1]
    return (Input)
        
        
Input=Remove_Text_After_Last_Occurence_of_Character()
print (Input)

#Split the dataframe into separate dataframes by the subModule name
def SplitInput_df_by_subModule():
    DF_Input_lst =[]
    for subModule in Input['Kinase_subModules'].unique():
        DF=Input.loc[Input['Kinase_subModules']==subModule]
        DF_Input_lst.append(DF)
    return DF_Input_lst

DF_Input_lst=SplitInput_df_by_subModule()
print (DF_Input_lst)


#Sort each dataframe within the list of dataframes by ascending for the FDR column
def Sort_by_Ascending():
    DF_Input_lst2=[]
    for DF in DF_Input_lst:
        DF=DF.copy()
        DF=DF.sort(['FDR'], ascending=[True])
        DF_Input_lst2.append(DF)
    return DF_Input_lst2

DF_Input_lst2=Sort_by_Ascending()



#Concatenate the dataframes back together into one so they can be printed out as a single dataframe.
def ConcatenateDFs():    #Concatenate the DFs together 
    EmptyDF = pd.DataFrame() # create an empty dataframe
    for df in DF_Input_lst2:  # select a dataframe in the list 
        df=df.copy() # make a copy of that dataframe 
        EmptyDF=EmptyDF.append(df) # append to the empty DF the dataframe selected and overwrite the empty dataframe
    return EmptyDF

Final=ConcatenateDFs()
#print (Final)


def DF_to_CSV(dataframe, NewFileName): 
    path ='/Users/mmacgilvray18/Desktop/' 
    dataframe.to_csv (path+NewFileName,sep='\t')
    
DF_to_CSV(Final, 'All_DTT_T120_Kinase_Module_FDR_Scores_and_their_Kinase_SI_subModules_Sorted_by_Ascending_Sept2017.csv')

