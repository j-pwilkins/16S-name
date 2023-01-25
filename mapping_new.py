import pandas as pd
import numpy as np
import os
import shutil
import sys
from datetime import datetime
pd.options.display.max_columns = None

def main():
    start_time = datetime.now()
    input_sequences_list, input_database, vsearch_output, curated_vsearch_output, output_directory, output_file, region = read_inputs()
    organise_directories(output_directory, input_sequences_list, input_database)
    database_df, query_df = prepare_then_run_vsearch(input_sequences_list, input_database, vsearch_output, region)

#L1
def map_vsearch_output_to_database(curated_vsearch_output, database_df, output_file, query_df, vsearch_output):
    parse_vsearch_usearchglobal_blast6(vsearch_output, curated_vsearch_output)
    add_similarity_columns(curated_vsearch_output)
    create_output_file(query_df, database_df, curated_vsearch_output, output_file)

#L1
def read_inputs():
    # input_sequences_list = sys.argv[1]
    input_sequences_list = 'test_30.csv'
    # input_database = sys.argv[2]
    input_database = 'db_497.csv'
    vsearch_output = 'vsearch_blast6_output.csv'
    curated_vsearch_output = 'Query_vs_All.csv'
    text_insert = '_vs_'
    region = 'unknown'
    output_directory = input_sequences_list.rstrip('.csv') + text_insert + input_database.rstrip('.csv')
    output_file = 'Summary_' + output_directory + '.csv'
    return input_sequences_list, input_database, vsearch_output, curated_vsearch_output, output_directory, output_file, region

#L1
def organise_directories(output_directory, input_sequences_list, input_database):
    if os.path.exists(output_directory):
        shutil.rmtree(output_directory)
        print(f'The previous {output_directory} folder has been deleted')
    os.makedirs(output_directory, exist_ok=True)
    shutil.copy(input_sequences_list, output_directory)
    shutil.copy(input_database, output_directory)
    os.chdir(output_directory)

#L1 # Read in .csv files, create fasta, run vsearch
def prepare_then_run_vsearch(input_sequences_list, input_database, vsearch_output, region):
    query_df, database_df, length_of_query, length_of_database, query_fastaname, database_fastaname = import_and_prepare_files(input_sequences_list, input_database, region)
    for i in range(length_of_query):
        convert_query_to_fasta(query_df, query_fastaname, i)
    for i in range(length_of_database):
        convert_db_to_fasta(database_df, database_fastaname, i)
    run_vsearch(query_fastaname, database_fastaname, vsearch_output)
    return database_df, query_df

#L2
def import_and_prepare_files(input_sequences_list, input_database, region):
    query_df = read_csv(input_sequences_list)
    query_df = query_df.replace(np.nan, '', regex=True)
    database_df = read_csv(input_database)
    length_of_query = len(query_df.index.tolist())
    length_of_database = len(database_df.index.tolist())
    query_fastaname, database_fastaname = create_fasta_filenames(input_sequences_list, input_database)
    query_df.insert(0, 'Query#', np.arange(1, len(query_df) + 1))
    print(f'The {input_sequences_list} of {length_of_query} {region} region 16S rRNA sequences has been read in. It will be compared against the database {input_database} of {length_of_database} full length 16S rRNA gene sequences. ')
    return query_df, database_df, length_of_query, length_of_database, query_fastaname, database_fastaname

#L3
def create_fasta_filenames(input_sequences_list, input_database):
    query_base_name = os.path.splitext(input_sequences_list)[0]
    query_fastaname = query_base_name + '.fa'
    database_base_name = os.path.splitext(input_database)[0]
    database_fastaname = database_base_name + '.fa'
    return query_fastaname, database_fastaname

#L2
def convert_query_to_fasta(query_df, query_fastaname, i):
    firstline = '>Query_' + str(query_df['Query#'].values[i]) + ' ' + str(query_df['DB Accession Number'].values[i]) + ' ' + str(query_df['DB Name'].values[i])
    secondline = str(query_df['16S rRNA sequence'].values[i])
    with open(query_fastaname, 'a') as f:
        f.write(f"{firstline}\n{secondline}\n")

#L2
def convert_db_to_fasta(database_df, database_fastaname, i):
    firstline = '>DB#_' + str(database_df['DB#'].values[i]) + ' ' + str(database_df['DB Accession Number'].values[i]) + " " + str(database_df['DB Name'].values[i])
    secondline = str(database_df['16S rRNA sequence'].values[i])
    with open(database_fastaname, 'a') as f:
        f.write(f"{firstline}\n{secondline}\n")

#L2 # Run vsearch - might need further function for --minwordmatches decision
def run_vsearch(query_fastaname, database_fastaname, vsearch_output):
    # cmd = 'vsearch --usearch_global ' + query_fastaname + ' --db ' + database_fastaname + ' --blast6out ' + vsearch_output + ' --acceptall --id 0.0 --maxaccepts 0 --minwordmatches 3'
    # os.system(cmd)
    print('vsearch has finished')

## helper functions
def read_csv(filename):
    try:
        df = pd.read_csv(filename)
        return df
    except FileNotFoundError:
        print("The file '{}' could not be found. Make sure the file is in the correct location and try again.".format(filename))
        exit()
    except pd.errors.EmptyDataError:
        print("The file '{}' is empty. Make sure the file contains data and try again.".format(filename))
        exit()
    except:
        print("An unknown error occurred while trying to read the file '{}'.".format(filename))
        exit()

def write_csv(df, filename):
    df.to_csv(filename, sep=',', index=False)

def check_unique_values(df1, df2, column_header):
    unique_values_df1 = df1[column_header].unique().tolist()
    unique_values_df1.sort()
    unique_values_df2 = df2[column_header].unique().tolist()
    unique_values_df2.sort()
    extra_unique_values_df1 = list(set(unique_values_df1) - set(unique_values_df2))
    extra_unique_values_df2 = list(set(unique_values_df2) - set(unique_values_df1))
    unique_values_in_both_lists = list(set(unique_values_df1).intersection(unique_values_df2))
    number_of_shared_unique_values = len(unique_values_in_both_lists)
    if unique_values_df1 == unique_values_df2:
        print(f"{df1.name} and {df2.name} are compatible. They contain {number_of_shared_unique_values} unique values in the '{column_header}' column.")
    elif len(extra_unique_values_df1) > 0 and len(extra_unique_values_df2) > 0:
        print(f"{df1.name} and {df2.name} both contain extra unique values in the '{column_header}' column. They share {number_of_shared_unique_values} unique values. {df1.name} has {len(extra_unique_values_df1)} extra unique values. {df2.name} has {len(extra_unique_values_df2)} extra unique values.")
    elif len(extra_unique_values_df1) > 0 and len(extra_unique_values_df2) == 0:
        print(f"{df1.name} contains extra unique values in the '{column_header}' column. They share {number_of_shared_unique_values} unique values. {df1.name} has {len(extra_unique_values_df1)} extra unique values.")
    elif len(extra_unique_values_df1) == 0 and len(extra_unique_values_df2) > 0:
        print(f"{df2.name} contains extra unique values in the '{column_header}' column. They share {number_of_shared_unique_values} unique values. {df2.name} has {len(extra_unique_values_df2)} extra unique values.")

def check_df_column_unique_entries_only(df, column_header):
    if df[column_header].duplicated().any():
        duplicated_rows = df[column_header].duplicated().sum()
        print(f"{duplicated_rows} queries in the {df.name} '{column_header}' column are duplicated")
    else:
        print(f"All entries in the {df.name} '{column_header}' column are unique")

def concatenate_dfs(*dfs):
    column_headers = [list(df.columns) for df in dfs]
    if all(headers == column_headers[0] for headers in column_headers):
        concatenated_df = pd.concat(dfs)
        return concatenated_df
    else:
        df_names = ', '.join([str(df) for df in dfs])
        print(f"{df_names} could not be concatenated as the column headers do not match. Please rectify")


main()