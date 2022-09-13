# load packages
import pandas as pd
import numpy as np
import os
import shutil
import sys
from datetime import datetime

# Grab current time before running code
startTime = datetime.now()

# Parse command line arguments
Inputfile = sys.argv[1]
Output_directory = sys.argv[2]

# Create new working directory
os.makedirs(Output_directory, exist_ok = True)  # Create directory
shutil.copy(Inputfile, Output_directory)    # Move Input file to new dir
os.chdir(Output_directory)  # Make new dir working dir

## Open input file & log number of queries to be performed
Input_df = pd.read_csv(Inputfile)
Runs = len(Input_df.index.tolist())   # Temp to make it easier to check
RowList = list(range(Runs))
fastaRowList = RowList[:]   # Copies list so it can be manipulated without affecting RowList

## create running Pandas dataframe
# Includes dummy 1st row that will be deleted later
d1, d2, d3 = {'Entry': [0], 'Db Species Name': ['Template'], 'A': [-1]}, dict.fromkeys('BCDEFGHIJKLM',[0]), {'Vs': [0], 'Similarity': [0], 'DB': 'Silva', 'Accession No.': [0], 'Date Accessed': [0], '16S Sequence': ['AAAAAGGGG']}
Output_df = pd.DataFrame(dict(d1, **d2, **d3))
Output_df.to_csv("rdf.csv", sep=',', index=False)

## convert Input sheet into fasta format for vsearch
Fastabase = Inputfile[:]    # Copy Inputfile
base = os.path.splitext(Fastabase)[0]   # remove .csv
Fastaname = base + '.fa'    # add .fa
def querytofasta():
    # parse csv file to create required variables
    Input_df = pd.read_csv(Inputfile)
    Row = fastaRowList[0]   # Take the 1st value from a list which sets the program to take the 1st row
    Query = str(Input_df['Input Index'].values[Row])
    Name = str(Input_df['DB Name'].values[Row])
    Sequence = str(Input_df['16S Sequence'].values[Row])
    # write to text file
    with open(Fastaname, 'a') as f:
        f.write(f">Q{Query} {Name}\n{Sequence}\n")
    # delete 1st entry from RowList so the next Row can be used
    del fastaRowList[0]
for i in range(Runs):
  querytofasta()

## Run vsearch using a system call
cmd = 'vsearch --allpairs_global ' + Fastaname + ' --blast6out vsearch.csv --acceptall'
os.system(cmd)

## convert vsearchoutput into format for this program
# Load 1st 3 columns of vsearch output
vsearch_df = pd.read_csv("vsearch.csv", usecols=[0, 1, 2], sep='\t', names=["Query Seq", "Target Seq", "Sim Target"])
# Add column to identify each row uniquely - based on their order within the vsearch output
firstorder = pd.Series(np.arange(1,46,1))
vsearch_df.insert(0, "Vsearch Entry", firstorder, True)
# Load vsearch output a 2nd time with a different name
vsearch2_df = pd.read_csv("vsearch.csv", usecols=[0, 1, 2], sep='\t', names=["Target Seq", "Query Seq", "Sim Target"])
# Switch Query and Target columns
vsearch2_df = vsearch2_df[vsearch2_df.columns[[1,0,2]]]
# Add column to identify each row uniquely - based on their order within the vsearch output
firstorder = pd.Series(np.arange(46,91,1))
vsearch2_df.insert(0, "Vsearch Entry", firstorder, True)
# Concatenate the 2 df
vsearch_df = pd.concat([vsearch_df, vsearch2_df], ignore_index=True, axis=0)
# Add dummy 1st row to allow next bit to work
vsearch_df.loc[-1] = [0, 'Q1', 'Q0', 0]   # adds new row
vsearch_df.index = vsearch_df.index + 1  # shifting index
vsearch_df = vsearch_df.sort_index()  # sorting by index
# Convert Query columns to allow comparisons - e.g. change Q1 to 1
vsearch_df ['Query Seq'] = vsearch_df ['Query Seq'].str.replace(r'\D', '')   # Strip all non-numeric values
vsearch_df ['Target Seq'] = vsearch_df ['Target Seq'].str.replace(r'\D', '')
vsearch_df ['Query Seq'] = vsearch_df ['Query Seq'].astype(int)   # convert strings to integers
vsearch_df ['Target Seq'] = vsearch_df ['Target Seq'].astype(int)

## Add Columns finding closest LOWER sequence & similarity%
# Create some parameters for function
vRuns = len(vsearch_df.index.tolist())   # Number of rows in input file
vRowList = list(range(vRuns))
simlist = [1000]   # Create dummy 1st entry in list
seqlist = [1000]
# Create function to generate list of similarities
def closestsequence():
    vRow = vRowList[0]
    vQuery = vsearch_df['Query Seq'].values[vRow]
    close1_df = vsearch_df[vsearch_df['Query Seq'] == vQuery]  # Find all entries for that Q
    close2_df = close1_df[close1_df['Query Seq'] > close1_df['Target Seq']]  # Confine to lower Queries
    max_sim = close2_df.loc[close2_df['Sim Target'].idxmax()]
    simclo = max_sim['Sim Target']
    closeq = max_sim['Target Seq']
    simlist.append(simclo)
    seqlist.append(closeq)
    del vRowList[0]
for i in range(vRuns):
    closestsequence()
simlist.pop(0)  # remove dummy 1st entry
seqlist.pop(0)
# add columns from list
vsearch_df ['Closest Seq'] = seqlist
vsearch_df ['Sim Close'] = simlist
vsearch_df = vsearch_df.iloc[1:]
# Re-order df to be more readable (group queries and order ascendingly)
vsearch_df = vsearch_df.sort_values(['Query Seq', 'Target Seq'])
# Save df to csv
vsearch_df.to_csv("AllvsAll.csv", sep=',', index=False)

## Generate new names for each query
def namequery():
    # Import RunningDataframe
    Output_df = pd.read_csv("rdf.csv")
    # Open file - create variables from Row that is being used as Input
    Input_df = pd.read_csv(Inputfile)
    Row = RowList[0]   # Take the 1st value from a list which sets the program to take the 1st row
    Input = Input_df['Input Index'].values[Row]
    Species_name = Input_df['DB Name'].values[Row]
    Db_Found = Input_df['DB Found'].values[Row]
    Accession = Input_df['DB Accession Number'].values[Row]
    Date_Accessed = Input_df['Date Accessed'].values[Row]
    Sequence = Input_df['16S Sequence'].values[Row]
    # Read from vsearch list to find closest sequence & its similarity
    close_df = pd.read_csv("AllvsAll.csv")
    Query = int(Input)
    ComparisonRow_df = close_df.loc[close_df['Query Seq'] == Query]
    ComparisonRow = ComparisonRow_df.values[0]
    Versus = int(ComparisonRow[4])
    similarity = float(ComparisonRow[5])
    ComparisonSequence = Output_df['Db Species Name'].values[Versus]

    # Assign labels for the individual values & convert to integers
    Ai, Bi, Ci, Di, Ei, Fi, Gi, Hi, Ii, Ji, Ki, Li, Mi = int(Output_df['A'].values[Versus]), int(Output_df['B'].values[Versus]), int(Output_df['C'].values[Versus]), int(Output_df['D'].values[Versus]), int(Output_df['E'].values[Versus]), int(Output_df['F'].values[Versus]), int(Output_df['G'].values[Versus]), int(Output_df['H'].values[Versus]), int(Output_df['I'].values[Versus]), int(Output_df['J'].values[Versus]), int(Output_df['K'].values[Versus]), int(Output_df['L'].values[Versus]), int(Output_df['M'].values[Versus])
    # Set the thresholds - 99.94% is the highest threshold required for a single nt difference across a full 16S sequence
    A, B, C, D, E, F, G, H, I, J, K, L, M = 60, 70, 80, 85, 90, 95, 98, 99, 99.5, 99.6, 99.7, 99.8, 99.9

    # Algorithm to create new name
    if similarity < A:  # if the query fails to meet the threshold then
        Ao = int(max(Output_df['A'].values)) + 1  # creates a new bin +1 to the highest currently existing in the db, within this category
        Bo = Co = Do = Eo = Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = Oo = Po = Qo = 0  # all others are 0 as this is a new bin
    else:
        Ao = Ai   # if query passes threshold then it lies within the same bin as the reference sequence
        if similarity < B:   # as with failure within cat A
            df_b = Output_df[Output_df['A'] == Ai]
            Bo = int(max(df_b['B'].values)) + 1
            Co = Do = Eo = Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = Oo = Po = Qo = 0
        else:
            Bo = Bi
            if similarity < C:
                df_c = Output_df[(Output_df['A'] == Ai) & (Output_df['B'] == Bi)]
                Co = int(max(df_c['C'].values)) + 1
                Do = Eo = Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = Oo = Po = Qo = 0
            else:
                Co = Ci
                if similarity < D:
                    df_d = Output_df[(Output_df['A'] == Ai) & (Output_df['B'] == Bi) & (Output_df['C'] == Ci)]
                    Do = int(max(df_d['D'].values)) + 1
                    Eo = Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = Oo = Po = Qo = 0
                else:
                    Do = Di
                    if similarity < E:
                        df_e = Output_df[(Output_df['A'] == Ai) & (Output_df['B'] == Bi) & (Output_df['C'] == Ci) & (Output_df['D'] == Di)]
                        Eo = int(max(df_e['E'].values)) + 1
                        Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = Oo = Po = Qo = 0
                    else:
                        Eo = Ei
                        if similarity < F:
                            df_f = Output_df[(Output_df['A'] == Ai) & (Output_df['B'] == Bi) & (Output_df['C'] == Ci) & (Output_df['D'] == Di) & (Output_df['E'] == Ei)]
                            Fo = int(max(df_f['F'].values)) + 1
                            Go = Ho = Io = Jo = Ko = Lo = Mo = No = Oo = Po = Qo = 0
                        else:
                            Fo = Fi
                            if similarity < G:
                                df_g = Output_df[(Output_df['A'] == Ai) & (Output_df['B'] == Bi) & (Output_df['C'] == Ci) & (Output_df['D'] == Di) & (Output_df['E'] == Ei) & (Output_df['F'] == Fi)]
                                Go = int(max(df_g['G'].values)) + 1
                                Ho = Io = Jo = Ko = Lo = Mo = No = Oo = Po = Qo = 0
                            else:
                                Go = Gi
                                if similarity < H:
                                    df_h = Output_df[(Output_df['A'] == Ai) & (Output_df['B'] == Bi) & (Output_df['C'] == Ci) & (Output_df['D'] == Di) & (Output_df['E'] == Ei) & (Output_df['F'] == Fi) & (Output_df['G'] == Gi)]
                                    Ho = int(max(df_h['H'].values)) + 1
                                    Io = Jo = Ko = Lo = Mo = No = Oo = Po = Qo = 0
                                else:
                                    Ho = Hi   # The rest of the letters need changin as above - left for time reasons
                                    if similarity < I:
                                        df_i = Output_df[(Output_df['A'] == Ai) & (Output_df['B'] == Bi) & (Output_df['C'] == Ci) & (Output_df['D'] == Di) & (Output_df['E'] == Ei) & (Output_df['F'] == Fi) & (Output_df['G'] == Gi) & (Output_df['H'] == Hi)]
                                        Io = int(max(df_i['I'].values)) + 1
                                        Jo = Ko = Lo = Mo = No = Oo = Po = Qo = 0
                                    else:
                                        Io = Ii
                                        if similarity < J:
                                            df_j = Output_df[(Output_df['A'] == Ai) & (Output_df['B'] == Bi) & (Output_df['C'] == Ci) & (Output_df['D'] == Di) & (Output_df['E'] == Ei) & (Output_df['F'] == Fi) & (Output_df['G'] == Gi) & (Output_df['H'] == Hi) & (Output_df['I'] == Ii)]
                                            Jo = int(max(df_j['J'].values)) + 1
                                            Ko = Lo = Mo = No = Oo = Po = Qo = 0
                                        else:
                                            Jo = Ji
                                            if similarity < K:
                                                df_k = Output_df[(Output_df['A'] == Ai) & (Output_df['B'] == Bi) & (Output_df['C'] == Ci) & (Output_df['D'] == Di) & (Output_df['E'] == Ei) & (Output_df['F'] == Fi) & (Output_df['G'] == Gi) & (Output_df['H'] == Hi) & (Output_df['I'] == Ii) & (Output_df['J'] == Ji)]
                                                Ko = int(max(df_k['K'].values)) + 1
                                                Lo = Mo = No = Oo = Po = Qo = 0
                                            else:
                                                Ko = Ki
                                                if similarity < L:
                                                    df_l = Output_df[(Output_df['A'] == Ai) & (Output_df['B'] == Bi) & (Output_df['C'] == Ci) & (Output_df['D'] == Di) & (Output_df['E'] == Ei) & (Output_df['F'] == Fi) & (Output_df['G'] == Gi) & (Output_df['H'] == Hi) & (Output_df['I'] == Ii) & (Output_df['J'] == Ji) & (Output_df['K'] == Ki)]
                                                    Lo = int(max(df_l['L'].values)) + 1
                                                    Mo = No = Oo = Po = Qo = 0
                                                else:
                                                    Lo = Li
                                                    if similarity < M:
                                                        df_m = Output_df[(Output_df['A'] == Ai) & (Output_df['B'] == Bi) & (Output_df['C'] == Ci) & (Output_df['D'] == Di) & (Output_df['E'] == Ei) & (Output_df['F'] == Fi) & (Output_df['G'] == Gi) & (Output_df['H'] == Hi) & (Output_df['I'] == Ii) & (Output_df['J'] == Ji) & (Output_df['K'] == Ki) & (Output_df['L'] == Li)]
                                                        Mo = int(max(df_m['M'].values)) + 1
                                                        No = Oo = Po = Qo = 0
                                                    else:
                                                        Mo = Mi

    # Append new line to pandas Output_df
    # First create a temporary df for the new line
    MaxEntry = max(Output_df['Entry'].values)   # Find the highest entry in the df
    NewEntry = MaxEntry + 1
    Updata = {'Entry': [NewEntry], 'Db Species Name': [Species_name], 'A': [Ao], 'B': [Bo], 'C': [Co], 'D': [Do], 'E': [Eo], 'F': [Fo], 'G': [Go], 'H': [Ho], 'I': [Io], 'J': [Jo], 'K': [Ko], 'L': [Lo], 'M': [Mo], 'Vs': [Versus], 'Similarity': [similarity], 'DB': [Db_Found], 'Accession No.': [Accession], 'Date Accessed': [Date_Accessed], '16S Sequence': [Sequence]}
    Temp_df = pd.DataFrame(Updata)  # create temporary df
    # Second, concatenate the temporary df to the df
    Output_df = pd.concat([Output_df, Temp_df], ignore_index=True, axis=0)
    #with pd.option_context('display.max_columns', None): print(Output_df)
    Output_df.to_csv("rdf.csv", sep=',', index=False)
    del RowList[0]      # delete 1st entry from RowList so the next Row can be used
for i in range(Runs):       # Run the program as many times as there are rows (query inputs) in the Input csv
  namequery()

# Delete the Template row
Output_df = pd.read_csv("rdf.csv")
Output_df = Output_df.drop(0)
# Edit df and save to Summary csv file
Output_df['Entry'] = 'Q' + Output_df['Entry'].astype(str)
Output_df['Vs'] = 'Q' + Output_df['Vs'].astype(str)
Output_df.to_csv("Summary.csv", sep=',', index=False)
# create Output directory and move files there
os.mkdir('data')
shutil.move('rdf.csv', 'data')
shutil.move(Fastaname, 'data')

# Grab time and print runtime of program
programlength = datetime.now() - startTime

# Summary report for Bash
print(f'{sys.argv[0]} has completed in {programlength}. {sys.argv[1]} was processed and the Summary.csv file can be found in the new {sys.argv[2]} directory.')