# 16S-name
7 Dec 2022
I have added two programs today as an idea of where my code is at.

create_db.py is the program to take the Hungate 1000 sequences (or anything else) and create a database.
This is basically (might be a couple of changes that need to be made as I haven't tested it) a working program.

map_variable_regions.py is the program that takes a sheet of sequences from a variable region (shouldn't matter which I don't think, but has been built using V1-V3) and maps them to a sequence within the db (created by create_db.py). This one is still being coded as there are a couple of problems. The vsearch output is not as I am expecting it to be and it's taking me a while to work out why.

CuratedHungate1000.csv is the file that would go into create_db.py, and db_497.csv is the file that would come out. CuratedHungate1000.csv would also be used for map_variable_regions.py (with different sequences). 


Old
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
