# load packages
import random as ra
import numpy as np
import pandas as pd
from numpy.random import seed
from numpy.random import randint

# Find how many times the program has to run
Input_df = pd.read_csv("Input10.csv")
num_of_rows = len(Input_df.index.tolist())
RowList = list(range(num_of_rows))

#create function of program
def jpprogram():
    # Import RunningDataframe
    Output_df = pd.read_csv("RunningDataframe.csv")

    # Open file
    Input_df = pd.read_csv("Input10.csv")
    Row = RowList[0]   # Take the 1st value from a list which sets the program to take the 1st row
    Input = Input_df['Input Index'].values[Row]
    Species_name = Input_df['DB Name'].values[Row]
    Accession = Input_df['DB Accession Number'].values[Row]
    Date_Accessed = Input_df['Date Accessed'].values[Row]
    Sequence = Input_df['16S Sequence'].values[Row]

    # Temporary section - here I am going to randomise the reference sequence that the query is being built against
    Current = len(Output_df.index.tolist())    # Find how many rows currently in df
    valuesfromCurrent = randint(0, Current, 1000)   # Produce list of random integers
    ra.shuffle(valuesfromCurrent)   # Shuffle list
    RandomRow = int(valuesfromCurrent[0])   # select 1st number from list & convert into integer

    # Assign labels for the individual values & convert to integers
    Ai, Bi, Ci, Di, Ei, Fi, Gi, Hi, Ii, Ji, Ki, Li, Mi, Ni, Oi, Pi, Qi = int(Output_df['A'].values[RandomRow]), int(Output_df['B'].values[RandomRow]), int(Output_df['C'].values[RandomRow]), int(Output_df['D'].values[RandomRow]), int(Output_df['E'].values[RandomRow]), int(Output_df['F'].values[RandomRow]), int(Output_df['G'].values[RandomRow]), int(Output_df['H'].values[RandomRow]), int(Output_df['I'].values[RandomRow]), int(Output_df['J'].values[RandomRow]), int(Output_df['K'].values[RandomRow]), int(Output_df['L'].values[RandomRow]), int(Output_df['M'].values[RandomRow]), int(Output_df['N'].values[RandomRow]), int(Output_df['O'].values[RandomRow]), int(Output_df['P'].values[RandomRow]), int(Output_df['Q'].values[RandomRow])
    Am, Bm, Cm, Dm, Em, Fm, Gm, Hm, Im, Jm, Km, Lm, Mm, Nm, Om, Pm, Qm = int(max(Output_df['A'].values)), int(max(Output_df['B'].values)), int(max(Output_df['C'].values)), int(max(Output_df['D'].values)), int(max(Output_df['E'].values)), int(max(Output_df['F'].values)), int(max(Output_df['G'].values)), int(max(Output_df['H'].values)), int(max(Output_df['I'].values)), int(max(Output_df['J'].values)), int(max(Output_df['K'].values)), int(max(Output_df['L'].values)), int(max(Output_df['M'].values)), int(max(Output_df['N'].values)), int(max(Output_df['O'].values)), int(max(Output_df['P'].values)), int(max(Output_df['Q'].values))

    # Process similarity variable
    # Temporarily this will be a random number
    # seed random number generator
    seed(1)
    # generate some integers
    values = randint(50, 99, 1000)
    ra.shuffle (values)
    FirstNumber_Sim = values[0]
    similarity = FirstNumber_Sim

    # Set the thresholds - 99.94% is the highest threshold required for a single nt difference across a full 16S sequence
    A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q = 60, 70, 80, 85, 90, 95, 98, 99, 99.5, 99.6, 99.7, 99.8, 99.9, 99.91, 99.92, 99.93, 99.94

    # Algorithm to create new name
    if similarity < A:   # if the query fails to meet the threshold then
        Ao = Am + 1   # creates a new bin +1 to the highest currently existing in the db, within this category
        Bo = Co = Do = Eo = Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = Oo = Po = Qo = 0   # all others are 0 as this is a new bin
    else:
        Ao = Ai   # if query passes threshold then it lies within the same bin as the reference sequence
        if similarity < B:   # as with failure within cat A
            Bo = Bm + 1
            Co = Do = Eo = Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = Oo = Po = Qo = 0
        else:
            Bo = Bi
            if similarity < C:
                Co = Cm + 1
                Do = Eo = Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = Oo = Po = Qo = 0
            else:
                Co = Ci
                if similarity < D:
                    Do = Dm + 1
                    Eo = Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = Oo = Po = Qo = 0
                else:
                    Do = Di
                    if similarity < E:
                        Eo = Em + 1
                        Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = Oo = Po = Qo = 0
                    else:
                        Eo = Ei
                        if similarity < F:
                            Fo = Fm + 1
                            Go = Ho = Io = Jo = Ko = Lo = Mo = No = Oo = Po = Qo = 0
                        else:
                            Fo = Fi
                            if similarity < G:
                                Go = Gm + 1
                                Ho = Io = Jo = Ko = Lo = Mo = No = Oo = Po = Qo = 0
                            else:
                                Go = Gi
                                if similarity < H:
                                    Ho = Hm + 1
                                    Io = Jo = Ko = Lo = Mo = No = Oo = Po = Qo = 0
                                else:
                                    Ho = Hi
                                    if similarity < I:
                                        Io = Im + 1
                                        Jo = Ko = Lo = Mo = No = Oo = Po = Qo = 0
                                    else:
                                        Io = Ii
                                        if similarity < J:
                                            Jo = Jm + 1
                                            Ko = Lo = Mo = No = Oo = Po = Qo = 0
                                        else:
                                            Jo = Ji
                                            if similarity < K:
                                                Ko = Km + 1
                                                Lo = Mo = No = Oo = Po = Qo = 0
                                            else:
                                                Ko = Ki
                                                if similarity < L:
                                                    Lo = Lm + 1
                                                    Mo = No = Oo = Po = Qo = 0
                                                else:
                                                    Lo = Li
                                                    if similarity < M:
                                                        Mo = Mm + 1
                                                        No = Oo = Po = Qo = 0
                                                    else:
                                                        Mo = Mi
                                                        if similarity < N:
                                                            No = Nm + 1
                                                            Oo = Po = Qo = 0
                                                        else:
                                                            No = Ni
                                                            if similarity < O:
                                                                Oo = Om + 1
                                                                Po = Qo = 0
                                                            else:
                                                                Oo = Oi
                                                                if similarity < P:
                                                                    Po = Pm + 1
                                                                    Qo = 0
                                                                else:
                                                                    Po = Pi
                                                                    if similarity < Q:
                                                                        Qo = Qm + 1
                                                                    else:
                                                                        Qo = Qi
    # Append new line to pandas Output_df
    # First create a temporary df for the new line
    MaxEntry = max(Output_df['Entry'].values)
    NewEntry = MaxEntry + 1

    # Match existing format & populate with correct data
    Updata = {'Entry': [NewEntry],
            'Db Species Name': [Species_name],
            'A': [Ao],
            'B': [Bo],
            'C': [Co],
            'D': [Do],
            'E': [Eo],
            'F': [Fo],
            'G': [Go],
            'H': [Ho],
            'I': [Io],
            'J': [Jo],
            'K': [Ko],
            'L': [Lo],
            'M': [Mo],
            'N': [No],
            'O': [Oo],
            'P': [Po],
            'Q': [Qo],
            'Vs': [RandomRow],
            'Similarity': [similarity],
            '16S Sequence': [Sequence]
            }

    Temp_df = pd.DataFrame(Updata)   # create temporary df
    # Second, concatenate the temporary df to the df
    Output_df = pd.concat([Output_df, Temp_df],ignore_index = True, axis = 0)

    # Print output - this is old code used to check that the output was working correctly
    #print(f"The thresholds used here are A={A}% B={B}% C={C}% D={D}% E={E}%")
    #print(f"Your query '{Species_name}' is {similarity}% similar to 'Template' - {Ai}A {Bi}B {Ci}C {Di}D {Ei}E {Fi}F {Gi}G {Hi}H {Ii}I {Ji}J {Ki}K {Li}L {Mi}M {Ni}N {Oi}O {Pi}P {Qi}Q")
    #print(f"{Ao}A {Bo}B {Co}C {Do}D {Eo}E {Fo}F {Go}G {Ho}H {Io}I {Jo}J {Ko}K {Lo}L {Mo}M {No}N {Oo}O {Po}P {Qo}Q is the name now assigned within the database to this query")

    #with pd.option_context('display.max_columns', None): print(Output_df)
    Output_df.to_csv("RunningDataframe.csv", sep=',', index=False)

    #delete 1st entry from RowList so the next Row can be used
    del RowList[0]

for i in range(num_of_rows):
  jpprogram()
