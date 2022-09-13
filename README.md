# 16S-name
Code for naming system for rumen 16S samples

To run this program Python 3 and vsearch need to be installed

To run the program enter the following:
python Prog.py Input.csv userdefined_dir

Input.csv is the name of the csv file that the program will be working on
userdefined_dir is the directory where the data will be stored and can be set by the user

N.B. The program currently generates a new database based on the inputted data, and currently does not support updating the database.

An additional program '16SGenerator' is included here that generates some randomish 16S sequences to test the main program 'Prog.py'.
To run this use 
python 16SGenerator.py number

The number here is the number of 16S sequences that the user wants the program to generate.
The output 'Input.csv' can then be fed into the main program 'Prog.py'
