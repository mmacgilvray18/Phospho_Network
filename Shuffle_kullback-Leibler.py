#!/home/mplace/anaconda3/bin/python
"""
Purpose:  


@author: mplace
"""
import multiprocessing
import argparse	                # handle command line args
import glob
import math
import numpy as np
import os
import pandas as pd
import random
import re
import sys

def SplitCompareTOMotifs_df(Input):
    ''' Function splits the Compare_To DF by Motif, which is listed in the "Motif" column, and puts the new dataframes into a list'''
    DF_CompareTo_lst =[]
    for Motif in Input['Motif'].unique():
        DF=Input.loc[Input['Motif']==Motif]
        DF_CompareTo_lst.append(DF)
    return DF_CompareTo_lst

def Shuffle(df):
    '''Shuffle each column by row, after first creating independent dataframes for each column (position within the PWM)''' 

    df_Index = df[['Motif','AA']]  
    df_data_0 = df['0']#,'1','2','3','4','5','7','8','9','10','11','12']]
    df_data_1 = df['1']#'2','3','4','5','7','8','9','10','11','12']]
    df_data_2 = df['2']
    df_data_3 = df['3']
    df_data_4 = df['4']
    df_data_5 = df['5']
    df_data_6 = df['6']
    df_data_7 = df['7']
    df_data_8 = df['8']
    df_data_9 = df['9']
    df_data_10 = df['10']
    df_data_11 = df['11']
    df_data_12 = df['12']
    Shuffled_Input_0=df_data_0.iloc[np.random.permutation(len(df_data_0))] 
    Shuffled_Input_1=df_data_1.iloc[np.random.permutation(len(df_data_1))]                                  # Shuffle the data by row
    Shuffled_Input_2=df_data_2.iloc[np.random.permutation(len(df_data_2))] 
    Shuffled_Input_3=df_data_3.iloc[np.random.permutation(len(df_data_3))] 
    Shuffled_Input_4=df_data_4.iloc[np.random.permutation(len(df_data_4))] 
    Shuffled_Input_5=df_data_5.iloc[np.random.permutation(len(df_data_5))] 
    Shuffled_Input_6=df_data_6.iloc[np.random.permutation(len(df_data_6))] 
    Shuffled_Input_7=df_data_7.iloc[np.random.permutation(len(df_data_7))] 
    Shuffled_Input_8=df_data_8.iloc[np.random.permutation(len(df_data_8))]
    Shuffled_Input_9=df_data_9.iloc[np.random.permutation(len(df_data_9))] 
    Shuffled_Input_10=df_data_10.iloc[np.random.permutation(len(df_data_10))] 
    Shuffled_Input_11=df_data_11.iloc[np.random.permutation(len(df_data_11))] 
    Shuffled_Input_12=df_data_12.iloc[np.random.permutation(len(df_data_12))] 

    Shuffled_Input_0_reset_Index= Shuffled_Input_0.reset_index(drop=True)                                   # reset the index for a later merge
    Shuffled_Input_1_reset_Index= Shuffled_Input_1.reset_index(drop=True)
    Shuffled_Input_2_reset_Index= Shuffled_Input_2.reset_index(drop=True) 
    Shuffled_Input_3_reset_Index= Shuffled_Input_3.reset_index(drop=True)
    Shuffled_Input_4_reset_Index= Shuffled_Input_4.reset_index(drop=True) 
    Shuffled_Input_5_reset_Index= Shuffled_Input_5.reset_index(drop=True)
    Shuffled_Input_6_reset_Index= Shuffled_Input_6.reset_index(drop=True) 
    Shuffled_Input_7_reset_Index= Shuffled_Input_7.reset_index(drop=True)
    Shuffled_Input_8_reset_Index= Shuffled_Input_8.reset_index(drop=True) 
    Shuffled_Input_9_reset_Index= Shuffled_Input_9.reset_index(drop=True)
    Shuffled_Input_10_reset_Index= Shuffled_Input_10.reset_index(drop=True)
    Shuffled_Input_11_reset_Index= Shuffled_Input_11.reset_index(drop=True) 
    Shuffled_Input_12_reset_Index= Shuffled_Input_12.reset_index(drop=True)
    
    result = pd.concat([df_Index, Shuffled_Input_0_reset_Index, Shuffled_Input_1_reset_Index, Shuffled_Input_2_reset_Index, Shuffled_Input_3_reset_Index,
                       Shuffled_Input_4_reset_Index, Shuffled_Input_5_reset_Index, Shuffled_Input_6_reset_Index, Shuffled_Input_7_reset_Index,
                       Shuffled_Input_8_reset_Index, Shuffled_Input_9_reset_Index, Shuffled_Input_10_reset_Index, Shuffled_Input_11_reset_Index, Shuffled_Input_12_reset_Index], axis=1) # concatenate the dataframes - they will sit side by side since have the same index (numbering)
    result_reordered=result[[ 'Motif', 'AA','0','1', '2','3','4','5','6','7','8','9','10','11','12']]       # reorder the columns so in the correct PWM order.
    
    result_reordered_Index=result_reordered[['Motif','AA']]
    result_reordered_Frame=result_reordered[['0','1', '2','3','4','5','6','7','8','9','10','11','12']]
    cols = result_reordered_Frame.columns.tolist()                                                          # send column headers to a list
   
    random.shuffle(cols)                                                                                    # shuffle the columns, which are in a list by name, and return a different order
 
    
    
    FinalDF=result_reordered_Frame[cols]                                                                    # make a new dataframe with randomly shuffled columns
    FinalDF.columns = ['0','1', '2','3','4','5','6','7','8','9','10','11','12']                             # reset the column names, so that they have the original names
    FinalDataframe= pd.concat([result_reordered_Index, FinalDF],axis=1)
    return FinalDataframe

def mergeInputMotifFile_withDF_CompareTo(df_Input,df_CompareTo):
    ''' Merge the query and comparison PWMs so that KLD can be calculated by comparing column values'''
    df_merged=df_Input.merge(df_CompareTo, on='AA')
    #df_lst.append(df_merged)
    return df_merged
   
def Calculate_log_x_y(df):
        ''' Function takes the log2 value of the amino acid frequency at each position of the query/comparison motifs'''
        df['0_log(x/y)'] = df.apply(lambda x: math.log(x['0_x'],2) - math.log(x['0_y'],2), axis=1)
        df['1_log(x/y)'] = df.apply(lambda x: math.log(x['1_x'],2) - math.log(x['1_y'],2), axis=1)
        df['2_log(x/y)'] = df.apply(lambda x: math.log(x['2_x'],2) - math.log(x['2_y'],2), axis=1)
        df['3_log(x/y)'] = df.apply(lambda x: math.log(x['3_x'],2) - math.log(x['3_y'],2), axis=1)
        df['4_log(x/y)'] = df.apply(lambda x: math.log(x['4_x'],2) - math.log(x['4_y'],2), axis=1)
        df['5_log(x/y)'] = df.apply(lambda x: math.log(x['5_x'],2) - math.log(x['5_y'],2), axis=1)
        df['6_log(x/y)'] = df.apply(lambda x: math.log(x['6_x'],2) - math.log(x['6_y'],2), axis=1)
        df['7_log(x/y)'] = df.apply(lambda x: math.log(x['7_x'],2) - math.log(x['7_y'],2), axis=1)
        df['8_log(x/y)'] = df.apply(lambda x: math.log(x['8_x'],2) - math.log(x['8_y'],2), axis=1)
        df['9_log(x/y)'] = df.apply(lambda x: math.log(x['9_x'],2) - math.log(x['9_y'],2), axis=1)
        df['10_log(x/y)'] = df.apply(lambda x: math.log(x['10_x'],2) - math.log(x['10_y'],2), axis=1)
        df['11_log(x/y)'] = df.apply(lambda x: math.log(x['11_x'],2) - math.log(x['11_y'],2), axis=1)
        df['12_log(x/y)'] = df.apply(lambda x: math.log(x['12_x'],2) - math.log(x['12_y'],2), axis=1)
        return df

def Calculate_log_y_x(df):
        ''' Function takes the log2 value of the amino acid frequency at each position of the comparison/query motifs'''
        df['0_log(y/x)'] = df.apply(lambda x: math.log(x['0_y'],2) - math.log(x['0_x'],2), axis=1)
        df['1_log(y/x)'] = df.apply(lambda x: math.log(x['1_y'],2) - math.log(x['1_x'],2), axis=1)
        df['2_log(y/x)'] = df.apply(lambda x: math.log(x['2_y'],2) - math.log(x['2_x'],2), axis=1)
        df['3_log(y/x)'] = df.apply(lambda x: math.log(x['3_y'],2) - math.log(x['3_x'],2), axis=1)
        df['4_log(y/x)'] = df.apply(lambda x: math.log(x['4_y'],2) - math.log(x['4_x'],2), axis=1)
        df['5_log(y/x)'] = df.apply(lambda x: math.log(x['5_y'],2) - math.log(x['5_x'],2), axis=1)
        df['6_log(y/x)'] = df.apply(lambda x: math.log(x['6_y'],2) - math.log(x['6_x'],2), axis=1)
        df['7_log(y/x)'] = df.apply(lambda x: math.log(x['7_y'],2) - math.log(x['7_x'],2), axis=1)
        df['8_log(y/x)'] = df.apply(lambda x: math.log(x['8_y'],2) - math.log(x['8_x'],2), axis=1)
        df['9_log(y/x)'] = df.apply(lambda x: math.log(x['9_y'],2) - math.log(x['9_x'],2), axis=1)
        df['10_log(y/x)'] = df.apply(lambda x: math.log(x['10_y'],2) - math.log(x['10_x'],2), axis=1)
        df['11_log(y/x)'] = df.apply(lambda x: math.log(x['11_y'],2) - math.log(x['11_x'],2), axis=1)
        df['12_log(y/x)'] = df.apply(lambda x: math.log(x['12_y'],2) - math.log(x['12_x'],2), axis=1)
        return df

def Calculate_Faax_times_log_x_y(df):
    ''' Function multiplies the frequency of an amino acid (Faax) "Xa" at a specific position in the query motif against the log(Xa/Ya) for that amino acid
     It is calculating this part of the function  "Xalog(Xa/Ya)" '''
    df['F(aax)*0_log(x/y)']=df['0_x']*df['0_log(x/y)']
    df['F(aax)*1_log(x/y)']=df['1_x']*df['1_log(x/y)']
    df['F(aax)*2_log(x/y)']=df['2_x']*df['2_log(x/y)']
    df['F(aax)*3_log(x/y)']=df['3_x']*df['3_log(x/y)']
    df['F(aax)*4_log(x/y)']=df['4_x']*df['4_log(x/y)']
    df['F(aax)*5_log(x/y)']=df['5_x']*df['5_log(x/y)']
    df['F(aax)*6_log(x/y)']=df['6_x']*df['6_log(x/y)']
    df['F(aax)*7_log(x/y)']=df['7_x']*df['7_log(x/y)']
    df['F(aax)*8_log(x/y)']=df['8_x']*df['8_log(x/y)']
    df['F(aax)*9_log(x/y)']=df['9_x']*df['9_log(x/y)']
    df['F(aax)*10_log(x/y)']=df['10_x']*df['10_log(x/y)']
    df['F(aax)*11_log(x/y)']=df['11_x']*df['11_log(x/y)']
    df['F(aax)*12_log(x/y)']=df['12_x']*df['12_log(x/y)']
    return df

def Calculate_Faay_times_log_y_x(df):
    ''' Function multiplies the frequency of an amino acid (Faay) "Ya" at a specific position in the query motif against the log(Ya/Xa) for that amino acid
     It is calculating this part of the function  "Yalog(Ya/Xa)" '''
    df['F(aay)*0_log(y/x)']=df['0_y']*df['0_log(y/x)']
    df['F(aay)*1_log(y/x)']=df['1_y']*df['1_log(y/x)']
    df['F(aay)*2_log(y/x)']=df['2_y']*df['2_log(y/x)']
    df['F(aay)*3_log(y/x)']=df['3_y']*df['3_log(y/x)']
    df['F(aay)*4_log(y/x)']=df['4_y']*df['4_log(y/x)']
    df['F(aay)*5_log(y/x)']=df['5_y']*df['5_log(y/x)']
    df['F(aay)*6_log(y/x)']=df['6_y']*df['6_log(y/x)']
    df['F(aay)*7_log(y/x)']=df['7_y']*df['7_log(y/x)']
    df['F(aay)*8_log(y/x)']=df['8_y']*df['8_log(y/x)']
    df['F(aay)*9_log(y/x)']=df['9_y']*df['9_log(y/x)']
    df['F(aay)*10_log(y/x)']=df['10_y']*df['10_log(y/x)']
    df['F(aay)*11_log(y/x)']=df['11_y']*df['11_log(y/x)']
    df['F(aay)*12_log(y/x)']=df['12_y']*df['12_log(y/x)']
    return df

def Column_SUM(df):
    ''' Function sums the values calculated by the previous two functions for each position, or column, in the PWMs  '''
    df['sum_0']=sum(df['F(aax)*0_log(x/y)'])+sum(df['F(aay)*0_log(y/x)'])
    df['sum_1']=sum(df['F(aax)*1_log(x/y)'])+sum(df['F(aay)*1_log(y/x)'])
    df['sum_2']=sum(df['F(aax)*2_log(x/y)'])+sum(df['F(aay)*2_log(y/x)'])
    df['sum_3']=sum(df['F(aax)*3_log(x/y)'])+sum(df['F(aay)*3_log(y/x)'])
    df['sum_4']=sum(df['F(aax)*4_log(x/y)'])+sum(df['F(aay)*4_log(y/x)'])
    df['sum_5']=sum(df['F(aax)*5_log(x/y)'])+sum(df['F(aay)*5_log(y/x)'])
    df['sum_6']=sum(df['F(aax)*6_log(x/y)'])+sum(df['F(aay)*6_log(y/x)'])
    df['sum_7']=sum(df['F(aax)*7_log(x/y)'])+sum(df['F(aay)*7_log(y/x)'])
    df['sum_8']=sum(df['F(aax)*8_log(x/y)'])+sum(df['F(aay)*8_log(y/x)'])
    df['sum_9']=sum(df['F(aax)*9_log(x/y)'])+sum(df['F(aay)*9_log(y/x)'])
    df['sum_10']=sum(df['F(aax)*10_log(x/y)'])+sum(df['F(aay)*10_log(y/x)'])
    df['sum_11']=sum(df['F(aax)*11_log(x/y)'])+sum(df['F(aay)*11_log(y/x)'])
    df['sum_12']=(sum(df['F(aax)*12_log(x/y)'])+sum(df['F(aay)*12_log(y/x)']))
    return df

def TotalScore(df):
    ''' Function calculates the total score by summing the summed values for each position in the PWM (13 positions)'''
    df['FinalScore']=df['sum_0']+df['sum_1']+df['sum_2']+df['sum_3']+df['sum_4']+df['sum_5']+df['sum_6']+df['sum_7']+df['sum_8']+df['sum_9']+df['sum_10']+df['sum_11']+df['sum_12']
    Lst=df['FinalScore'].unique()
    n=Lst[0]
    return n

def runShuffle(iterations, DF_CompareTo_lst, dfs_lst, n ):
    """
    """
        # CHANGE THE NUMBER OF ITERATIONS HERE IF DESIRED
    ITER_NUM=iterations                               # interations of this function
    dict_Final={}
    for df2 in dfs_lst:                                                                 # this is the dataframe that has Mok Kinase PWMs
        
        subModule_name=df2['Motif'].unique()
        for df in DF_CompareTo_lst:                                                     # select one of the compare to dataframes (Modules)
            Kinase_name=[]
            Kinase_name=df['Motif'].unique()
       
            for iteration in range (ITER_NUM):                                          # for iteration x 
        
                Shuffled=Shuffle(df2)                                                   # shuffle the dataframe by row and column
                df_merged=mergeInputMotifFile_withDF_CompareTo(Shuffled, df)   
                DF_1=Calculate_log_x_y(df_merged)
                DF_2=Calculate_log_y_x(DF_1)
                DF_3=Calculate_Faax_times_log_x_y(DF_2)
                DF_4=Calculate_Faay_times_log_y_x(DF_3)
                DF_5=Column_SUM(DF_4)
            
                n=TotalScore(DF_5)
                #print (n)
                test_tup = (n, subModule_name[0])
                if Kinase_name[0] in dict_Final:
                    dict_Final[Kinase_name[0]].append(test_tup)
                else:
                    lst=[]
                    dict_Final[Kinase_name[0]] = lst
                    dict_Final[Kinase_name[0]].append(test_tup)
    
    path=r"Shuffle_TEST/"  # path to output directory
                  
    for k, v in dict_Final.items():  
        newFile=path+ k + '.csv'       
        with open(newFile, 'a') as output:  
            output.write(k)
            output.write("\n")
            for x in v:
                #for subModule_name in filenames:
                output.write(str(x))
                output.write("\n")


def main():
    """
    Process command line arguments and run program.
    """    
    cmdparser = argparse.ArgumentParser(description="Shuffle kullback-Leibler results for use w/ FDR function.",
                                        usage='%(prog)s -f <Position weight matrix file>  ', prog='Shuffle_kullback-Leibler.py'  )                                  
    cmdparser.add_argument('-f', '--file', action='store', dest='FILE', 
                           help='Position Weight Matrix file, from Matt\'s script', metavar='')
    cmdResults = vars(cmdparser.parse_args())
    
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)
                
    # Get the inputfile, which contains a list of file names to process
    if cmdResults['FILE']:
        Input = pd.read_csv( cmdResults['FILE'])
    else:
        print('')
        cmdparser.print_help()
        sys.exit(1)    
        
    DF_CompareTo_lst=SplitCompareTOMotifs_df(Input)
    
    # Import the Mok Kinases PWM .csv files individually and create dataframes
    filenames = glob.glob('/home/mplace/projects/forMatt/Phospho_Network/Mok_kinase_PWMs/' + '*.csv') 
    
    dfs_lst = []
    for filename in filenames:
        dfs_lst.append(pd.read_csv(filename, sep=","))
    
    # write out the final dictionary to a folder where each key and value pair is an independent csv file. 
    if not os.path.exists('Shuffle_TEST'):
        os.mkdir('Shuffle_TEST')

    procs = 3
    iterations = 5    
    jobs = []
    for i in range(0,procs):
        process = multiprocessing.Process(target=runShuffle, args=(iterations,DF_CompareTo_lst, dfs_lst, str(i)))
        jobs.append(process)
    
    for j in jobs:
        j.start()
    
    for j in jobs:
        j.join()
        
    

    # run shuffle
    #runShuffle(5, DF_CompareTo_lst, dfs_lst) 


if __name__ == "__main__":
    main()