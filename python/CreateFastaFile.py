import pandas as pd

'''Script written by Matt MacGilvray 10/26/15 

As input, script takes a file containing Modules, their constituent Ppep names, and their Ppep amino acid sequences.
The script then splits the file by Module name into individual dataframes.
The script then creates individual Fasta files from each dataframe (module) and saves them into a user defined directory. 
The output Fasta format files are named with their module designation.


'''
Input_df=pd.read_csv('/Users/mmacgilvray18/Desktop/DTT_Modules_for_making_FASTAfiles.csv')



def Split_Into_SeparateDFs():
    ''' Function splits the input dataframe, based on the module name, into independent dataframes for each module'''
    df_lst=[]
    for Module in Input_df['Module'].unique():
        DF=Input_df.loc[Input_df['Module']==Module]
        df_lst.append(DF)
        
    return df_lst

df_lst=Split_Into_SeparateDFs()
print (df_lst)


def CreateIndividualFastaFiles():
    '''Function creates individual fasta files for each module nd writes them out to a user defined directory'''
    for df in df_lst:                
        Module_lst=df["Module"].tolist()    
        for name in Module_lst:          #
            ofile= open("/Users/mmacgilvray18/Desktop/DTT_Single_Rep_FastaFiles_Modules/"+name+".txt", "w")             # open a new file that contains the module name. USER can Change directory here.
        
            df_lstName=df['Name'].tolist()                                                                              # send the module names to a list
            df_lstSeq=df['Sequence'].tolist()                                                                           # send the peptide sequences to a list 
            
            for i in range(len(df_lstSeq)):                    
                
                ofile.write(">" + df_lstName[i] + "\n" + df_lstSeq[i] + "\n")                                            # create a fasta file where the peptide name will be followed by the peptide sequence, on a new line
           
        df_lstName=[]                                                                                                    # empty each of the lists for the next iteration
        df_lstSeq=[]
        Module_lst=[]
        ofile.close      
       
        
    return 
      
   
CreateIndividualFastaFiles()

