#!/home/mplace/anaconda3/bin/python
"""
Program: motifx.py

Purpose: Automate submitting jobs to Motif-x website. http://motif-x.med.harvard.edu/motif-x.html
         Each job will be run with 3 central characters S, T, Y.  This means each input file will
         have 3 different central characters used.  
         
         NOTE: The default SGD proteome file used by motifx is different (older) 
         than the current SGD file.  This results in a few motifs which cannot be 
         mapped back to identify the gene.  As a work around the input peptides
         are pre-aligned to the R64-2-1 version of SGD proteome (current as of 12/9/2016).
         
         Duplicate peptides are discarded to avoid biasing the motifx.
         for example:
             YOR051C_S33 
             YOR051C_S33_T34
             
             will produce 3 peptides, 2 of which are identical, one of these
             is discarded.

         located here: /home/GLBRCORG/mplace/scripts/motifx/orf_trans_all.20150113.fasta  
         
Input : One or more peptide excel files, listed in a text file, column order is 
        unimportant but column name must match, each file looks like

        Ppep            Group      Localized_Sequence        Motif_X_Input_Peptide
        YGL076C_T8_S11  Induced    AAEKILtPEsQLKK            AAEKILT*PES*QLKK
        YNR047W_T428    Induced    AASEPNGLQLASATSPtSSSAR    AASEPNGLQLASATSPT*SSSAR
    
Output : A final text table named after the input file containing all the motifs matched to a gene.
        
         Given an input file named,  motifx_sample.xlsx the final results file
         will be: motifx_sample-Motifx-results.txt
         
         The other results are put in a directory.  For instance if your input file is called
         motifx_sample.xlsx, 3 directories will be created one for each central character.
         
         motifx_sample_T
         motifx_sample_S
         motifx_sample_Y
         
         These contain the LOGO pngs and the original html results page.
         
NOTE: if the output file exists, it will not be overwritten and the program will exit.

Mofif-x input parameters: 
            upload pre-aligned
            central character
            width       13
            occurances  10
            significance 0.000001
            upload background  SGD R64-2-1 proteome fasta
            background fasta
            background central character
            
Dependencies: Python 3 
              Python modules: argparse, BeautifulSoup4, requests, xlrd, unicodecsv
              MotifxPreAlign.py (Pre-aligns peptide sequence via exact match)

author: Mike Place
Date:   12/12/2016
"""
import argparse	                # handle command line args
from bs4 import BeautifulSoup       # html parser
from collections import OrderedDict 
from collections import defaultdict
import MotifxPreAlign as mpa        # written by Mike Place to pre-align & extend peptides for input to motifx
import os
import re                           # regex
import requests                     # HTTP library
import shutil                       # 
import sys
import time                         # sleep
from tqdm import tqdm               # progress bar
import unicodecsv                   # handle unicode chars from excel file
import xlrd                         # handle excel files

class Motifx ( object ):
    """
    Methods and data structures for Motifx class 
    """
    def __init__(self, file, occ, sig, prot, wid ):
        """
        Set up Motifx object
        text        : text version of input excel file
        pep         : list of peptides to submit to motifx website
        occurrence  : occurrences, default is 10
        sig         : significance threshold, default = 0.000001
        proteome    : Yeast ORF fasta file to use, default is the motifx website default
        width       : width of motif, default = 13       
        """
        self.occurance   = occ
        self.sig         = sig
        self.proteome    = prot
        self.width       = wid
        self.dir         = os.getcwd()
        #self.text        = self.excelToText(file)    # get text version of input file
        self.text        = file                       # get input file
        self.centralRes  = ['Y','T','S']              # central character on motif-x web form required by motifx site
        self.fileName    = re.sub(r'.csv', '', file) # filename 
        # pre-align peptides
        self.prealign    = mpa.MotifxPreAlign(self.proteome,   self.text, self.width)
        self.prealign.matchSeq()
        self.prealign.cleanResult()
        # map peptide to gene name in dictionary
        self.mapPep      = self.pepTideToGene()       # looks like 'ALSRSPSNQQYLL': 'YIL135C_S303'
        self.group       = self.getGroup()
        self.logo        = []                         # list of locations for logo images
        self.peptideFile = 'pepFile.txt'              # temp peptide file name        
        self.result      = OrderedDict()
        
    def pepTideToGene(self):
        """ Create a dict w/ the key = the peptide and the value = gene name
            Used to output final results table. 
        """
        genePep = dict()                        # key = peptide, value = gene 
        for gene, pep in self.prealign.pepInfo.items():
            gninfo = gene.split('_')         # split gene name & position  YIL135C_S303 is broken up into a list
            name = gninfo.pop(0)             # capture gene name
            for p,g in zip(pep['extended'], gninfo):    # recombine name and appropriate position
                gname = name + '_' + g
                genePep[str(p)] = gname                
        return genePep               
        
    def getGroup(self):
        """
        Get Group type information induced/repressed 
        """
        with open(self.text, 'r') as f:
            for line in f:
                if line.startswith('Ppep'):   # skip header if present
                    line = f.readline()
                    break
        row = line.split(',')
        return row[1]
       
    def pepFile(self):
        """
        Write peptide to file to use when submitting job to Motifx website.
        This is just a temporary file, will be overwritten if multiple files are processed.
        The last processed file's peptides will be in this file.
        """
        with open('pepFile.txt', 'w') as pf:
            for key in self.mapPep.keys():
                pf.write("%s\n" %(key))
        pf.close()

    def excelToText( self, file ):
        """
        Converts an Excel file to a CSV file.
        If the excel file has multiple worksheets, only the first worksheet is converted.
        Uses unicodecsv, so it will handle Unicode characters.
        Uses a recent version of xlrd, so it should handle old .xls and new .xlsx equally well.
        """        
        wb = xlrd.open_workbook(file)
        sh = wb.sheet_by_index(0)

        fh = open('Temp_file.txt',"wb")
        csv_out = unicodecsv.writer(fh, encoding='utf-8')
        
        # go through the file line by line
        for row_number in range (sh.nrows):
            csv_out.writerow(sh.row_values(row_number))

        fh.close()
        return 'Temp_file.txt'
        
    def submitMotifX( self, char ):
        """
        Submit protein list to the Motif-x website.
                
        url:  http://motif-x.med.harvard.edu/cgi-bin/multimotif-x.pl
        
        result url looks similar to:
        http://motif-x.med.harvard.edu/  plus  /cgi-bin/jobres.pl?jobid=20160830-8155-04402686 (example only)
        which is found in tag:  <A HREF="/cgi-bin/jobres.pl?jobid=20160830-8155-04402686" TARGET="_blank">Check results</A>
        
        Inputs to the web form:     
            fgfile        -- upload file        # prealigned peptides
            fgtype        -- prealigned
            fgcentralres  -- central char  
            width         -- width
            occurances    -- occurrences
            bgdb          -- uploaded 
            bgtype        -- type fasta in this case
            bgcentralres  -- central Character
            bgfile        -- proteome fasta            
            
        return the url of the result page
        """
        url       =  'http://motif-x.med.harvard.edu/cgi-bin/multimotif-x.pl'           # submission url
        
        # set values required for the form submission  'submit':'get          
        form_data = {  'fgtype' : 'prealigned', 'fgcentralres' : char, 'width' : self.width,
                     'occurrences' : self.occurance, 'significance': self.sig, 'bgdb' : 'uploaded', 
                     'bgtype' : 'fasta', 'bgcentralres' : char}  
            # Now we have 2 files to submit to motifx
        file     = { 'fgfile': open( self.dir + '/' + self.peptideFile, 'rb'),           # peptide upload file, one peptide per line
                         'bgfile': open( self.proteome, 'rb')   }           # user provided proteome fasta
                
        # send value to website, post will time out after 6 min seconds with no response 
        try:
            response  = requests.post( url, data = form_data, files = file, allow_redirects = True, timeout = 360 )
            soup      = BeautifulSoup( response.text, 'html.parser')     # store resulting webpage as text, need to parse to get result page
            tag       = str(soup.a)
            jobID     = re.findall('"([^"]*)"', tag )       # get the results page location jobID[0]
        except requests.exceptions.RequestException as e:
            print(e)
            sys.exit(1)
        
        time.sleep(240)                                  # give motifx time to run the job
  
        # return the results location
        return 'http://motif-x.med.harvard.edu/' + jobID[0] 
                  
    def parseResults( self, resultURL, char ):
        """
        Parse the results page from a Motif-x job, resultURL is the results webpage.
        Write the html results page to a file.
        Download all the logo images to current directory.
        """
        centralChar = char
        # make directory for each central character's results
        resultsDir = self.fileName + '_' + centralChar        
        os.mkdir( resultsDir )
        
        # get the results page
        try:
            response = requests.get(resultURL, timeout = 60)          # html results 
        except requests.exceptions.RequestException as e:             
            print(e)
            sys.exit(1)
            
            
        # write the html page to file, to provide a record of query.
        with open( resultsDir + '/' + self.fileName + '-' + centralChar + '.html', 'w') as html:
            html.write(resultURL)
            html.write(response.text)
       
        # parse the returned result html
        tree = BeautifulSoup( response.text, "lxml" )
        
        # check if there are results to report
        check = tree.getText('BODY').split('\n')
        if not 'Motifs Found: None'in check:
            # get the body of the page
            data = tree.body.find_all('font')[3]
            
            # find motif asscociated with peptide
            for m in data.findNextSiblings('a'):
                if m.text not in self.result:          # add key to dictionary, key is the motif ex: ".....T....S"
                    self.result[m.text] = []
                for i in m.children:
                    self.result[m.text].append(str(i.next))         # append a string will all the motifs combined
          
            # Get logo png files
            for logo in tree.find_all('a'):
                for img in logo.find_all('img'):
                    path = 'http://motif-x.med.harvard.edu/' + img['src']      # get web address for logo image
                    req = requests.get(path, stream = True)                    # download image
                    if req.status_code == 200:       
                        imageFile = re.sub(r'/logos', '', img['src']) 
                        # good to go, now download
                        with open( resultsDir + imageFile, 'wb') as png:
                            req.raw.decode_content = True
                            shutil.copyfileobj(req.raw, png)
                    else:                                                      # if unable to download logo image record path
                        with open("logo-image.log", 'a') as log:
                            log.write("Unable to download image: %s" %(path))
        else:
            cntlChar = [ i for i in check if i.startswith('fgcentralres')]
            
            with open('LOG.txt', 'a') as log:
                log.write('No results returned for Motifx using Central Character: %s' %(cntlChar[0]) )
            log.close()
        
                
    def writeResults( self ):
        """
        Write peptide table, 3 columns, comma separated.
        peptide,Group,motif
        """
        if bool(self.result):     # check to see if we have any results to write out
            with open(self.fileName + '-Motifx-results.txt', 'w') as out:
                for k,v in self.result.items():
                    motifs = v[0].split('\n')
                    [ motifs.pop(0) for x in range(2)]
                    [ motifs.pop() for x in range(3)]
                    for row in motifs:
                        line = "".join(row)                           # this joins the peptides into a single string
                        gene = self.mapPep[line]
                        line = gene + "," + line + "," + self.group + "," + k
                        out.write("%s\n" %(line))
      
            
                
def main():
    """
    Process command line arguments and run program.
    """
    cmdparser = argparse.ArgumentParser(description="Automate submitting jobs to Motif-x website.",
                                        usage='%(prog)s -f <File listing excel files>  ', prog='Motifx.py'  )                                  
    cmdparser.add_argument('-f', '--file',      action='store', dest='FILE',     help='file listing excel files to process expected', metavar='')
    cmdparser.add_argument('-o', '--occurrence',action='store', dest='OCC',      help='occurrences, default is 10',  metavar='')
    cmdparser.add_argument('-s', '--sig',       action='store', dest='SIG',      help='significance, default = 0.000001', metavar='', type=float)
    cmdparser.add_argument('-u', '--upload',    action='store', dest='UPLOAD',   help='Upload Yeast ORF fasta file to use', metavar='')
    cmdparser.add_argument('-w', '--width',     action='store', dest='WIDTH',    help='width, default = 13', metavar='', type=int)      
    cmdparser.add_argument('-i', '--info',      action='store_true', dest='INFO',help='Detailed description of program.')
    cmdResults = vars(cmdparser.parse_args())
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)
             
    if cmdResults['INFO']:
        print("\n  Motifx.py ")
        print("\n  Automate submitting jobs to Motif-x website.")
        print("\n  Input  : Plain text file listing excel files to process, one excel file name per line.")
        print("\n\tExcel file format:")
        print("\tPpep	Group	Localized_Sequence	Motif_X_Input_Peptide")
        print("\tYGL076C_T8_S11	Induced	AAEKILtPEsQLKK	AAEKILT*PES*QLKK")
        print()
        print("\tColumn order is unimportant, column names must match above.")       
        print()
        print(" usage: Motifx.py -f inputfiles ")
        print()
        print(" Optional Arguments:")
        print("\t-o Minimum number of times each of your extracted motifs to occur in the data set (10)")
        print("\t-s P-value threshold for the binomial probability (.000001)")
        print("\t-u upload your own version of SGD proteome (orf_trans.fasta).")
        print("\t-w Number of total characters in motif, (13)")
        print("\n  Output :")
        print("\tA final text table named after the input file containing all the motifs matched to a gene.")
        print("\tGiven an input file named,  motifx_sample.xlsx the final results file")
        print("\twill be: motifx_sample-Motifx-results.txt ")
        print("\n\tThe other results are put in a directory." )
        print("\tFor instance if your input file is called motifx_sample.xlsx\n")
        print("\t3 directories will be created one for each central character:")
        print("\n\t\tmotifx_sample_T    motifx_sample_S    motifx_sample_Y")
        print("\n\tThese contain the LOGO pngs and the original html results page.")
        print("\n\tTo see Python Docs and get a better explaination of the program:")
        print("\n\tOpen python console and enter")
        print("\timport sys")
        print("\tsys.path.append('/full/path/to/script')")
        print("\timport Motifx")
        print("\thelp(Motifx)")
        print("\n\tSee Mike Place for help or suggestions.\n")
        sys.exit(1)
    
    # Get the inputfile, which contains a list of file names to process
    if cmdResults['FILE']:
        userFile = cmdResults['FILE']
    else:
        print('')
        cmdparser.print_help()
        sys.exit(1)
        
    # check for the occurance parameter 
    if cmdResults['OCC']:
        occurance = cmdResults['OCC']
        if not occurance.isdigit():
            print('\n\t-o occurance is a number, default = 10 \n')
            print('\tThe occurrence threshold refers to the minimum number of times'
                  '\n\tyou wish each of your extracted motifs to occur in the data set.\n')
            cmdparser.print_help()
            sys.exit(1)
    else:
        occurance = 10
            
    # check for significance parameter
    if cmdResults['SIG']:
        sig = float(cmdResults['SIG'])
    else:
        sig = '{:f}'.format(0.000001)
        
    # check for an alternate proteome file, a value of default will use the motifx default SGD proteome
    if cmdResults['UPLOAD']:
        proteome = cmdResults['UPLOAD']
        proteome = os.getcwd() + '/' + proteome    # user provides a proteome fasta
        if not os.path.exists(proteome):
            print('\n\tAlternate proteome file does not exist %s.\n' %(proteome))
            cmdparser.print_help()
            sys.exit(1)
    else:
        proteome = os.getcwd() + '/reference/orf_trans_all.20150113.fasta'  # use SGD R64-2-1 proteome fasta
    
        
    # check of an alternate window width
    if cmdResults['WIDTH']:
        width = cmdResults['WIDTH']
    else:
        width = 13
        
    # open input file & process each file listed
    with open(userFile, 'r') as f:
        for line in f:
            line = line.rstrip()
            print("Processing file: %s" %(line))                # prints the name of the file currently being processed
            d = Motifx(line, occurance, sig, proteome, width )  # create a Mofifx object 
            d.pepFile()                                         # write peptide to temp file for use in job submission.            
            # submit job to motifx website
            for char in tqdm(d.centralRes):                     # loop through all potential central characters
                resultsPage = d.submitMotifX(char)
                d.parseResults(resultsPage, char)

            d.writeResults()                                    # write out results to file
            
    
if __name__ == "__main__":
    main()
