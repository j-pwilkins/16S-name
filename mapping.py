import pandas as pd
import numpy as np
import os
import shutil
import sys
from datetime import datetime

def main():
    start_time = datetime.now()
    input_sequences_list = sys.argv[1]
    input_database = sys.argv[2]
    vsearch_output = 'vsearch_blast6_output.csv'
    curated_vsearch_output = 'Query_vs_All.csv'
    text_insert = '_vs_'
    output_directory = input_sequences_list.rstrip('.csv') + text_insert + input_database.rstrip('.csv')
    output_file = 'Summary_' + output_directory + '.csv'
    if os.path.exists(output_directory):
        shutil.rmtree(output_directory)
        print(f'The previous {output_directory} folder has been deleted')
    os.makedirs(output_directory, exist_ok=True)
    shutil.copy(input_sequences_list, output_directory)
    shutil.copy(input_database, output_directory)
    os.chdir(output_directory)
    import_csv(input_sequences_list, input_database, output_file)
    query_df, database_df, length_of_query, length_of_database, query_fastaname, database_fastaname = import_csv(
        input_sequences_list, input_database, output_file)
    for i in range(length_of_query):
        convert_query_to_fasta(query_df, query_fastaname, i)
    for i in range(length_of_database):
        convert_db_to_fasta(database_df, database_fastaname, i)
    run_vsearch (query_fastaname, database_fastaname, vsearch_output)
    parse_vsearch_usearchglobal_blast6(vsearch_output, curated_vsearch_output)
    add_similarity_columns(curated_vsearch_output)
    create_output_file(query_df, database_df, curated_vsearch_output, output_file)
    tidy_files(curated_vsearch_output, query_fastaname, database_fastaname, vsearch_output)
    summary_statement(output_directory, start_time)

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

def import_csv(input_sequences_list, input_database, output_file):
    query_df = read_csv(input_sequences_list)
    database_df = read_csv(input_database)
    query_df = query_df.replace(np.nan, '', regex=True)
    length_of_query = len(query_df.index.tolist())
    length_of_database = len(database_df.index.tolist())
    query_base_name = os.path.splitext(input_sequences_list)[0]
    query_fastaname = query_base_name + '.fa'
    database_base_name = os.path.splitext(input_database)[0]
    database_fastaname = database_base_name + '.fa'
    query_length_as_list = list(range(length_of_query))
    input_list = [x + 1 for x in query_length_as_list]
    query_df.insert(0, 'Query#', input_list)
    return query_df, database_df, length_of_query, length_of_database, query_fastaname, database_fastaname

def convert_query_to_fasta(query_df, query_fastaname, i):
    firstline = '>Query_' + str(query_df['Query#'].values[i]) + ' ' + str(query_df['DB Accession Number'].values[i]) + ' ' + str(query_df['DB Name'].values[i])
    secondline = str(query_df['16S rRNA sequence'].values[i])
    with open(query_fastaname, 'a') as f:
        f.write(f"{firstline}\n{secondline}\n")

def convert_db_to_fasta(database_df, database_fastaname, i):
    firstline = '>DB#_' + str(database_df['DB#'].values[i]) + ' ' + str(database_df['DB Accession Number'].values[i]) + " " + str(database_df['DB Name'].values[i])
    secondline = str(database_df['16S rRNA sequence'].values[i])
    with open(database_fastaname, 'a') as f:
        f.write(f"{firstline}\n{secondline}\n")

def run_vsearch(query_fastaname, database_fastaname, vsearch_output):
    cmd = 'vsearch --usearch_global ' + query_fastaname + ' --db ' + database_fastaname + ' --blast6out ' + vsearch_output + ' --acceptall --id 0.0 --maxaccepts 0 --minwordmatches 3'
    os.system(cmd)
    print('vsearch has finished')

def parse_vsearch_usearchglobal_blast6(vsearch_output, curated_vsearch_output):
    query_vs_targetdb_df = pd.read_csv(vsearch_output, usecols=[0, 1, 2], sep='\t',
                                            names=["Query", "DB#", "Similarity(%)"])
    df_length = len(query_vs_targetdb_df)
    vsearch_output_order = pd.Series(np.arange(1, df_length + 1, 1))
    query_vs_targetdb_df.insert(0, 'Vsearch Order', vsearch_output_order, True)
    query_vs_targetdb_df['Q'] = query_vs_targetdb_df['Query'].str[6:].astype(int)
    query_vs_targetdb_df['T'] = query_vs_targetdb_df['DB#'].str[4:].astype(int)
    query_vs_targetdb_df = query_vs_targetdb_df.sort_values(['Q', 'T'])
    query_vs_targetdb_df.to_csv(curated_vsearch_output, sep=',', index=False)

def add_similarity_columns(curated_vsearch_output):
    curated_df = pd.read_csv(curated_vsearch_output)

    def closest_match_to_query(input):
        selected_df = curated_df[:]
        selected_df = selected_df.loc[selected_df['Query'] == input]
        selected_df = selected_df.loc[selected_df['Q'] != selected_df['T']]
        match_row = selected_df.loc[selected_df['Similarity(%)'].idxmax()]
        closest_match = match_row['DB#']
        return closest_match

    def highest_similarity(input):
        selected_df = curated_df[:]
        selected_df = selected_df.loc[selected_df['Query'] == input]
        selected_df = selected_df.loc[selected_df['Q'] != selected_df['T']]
        match_row = selected_df.loc[selected_df['Similarity(%)'].idxmax()]
        highest_similarity = match_row['Similarity(%)']
        return highest_similarity

    def number_of_equal_highest_similarities(input):
        selected_df = curated_df[:]
        selected_df = selected_df.loc[selected_df['Query'] == input]
        selected_df = selected_df.loc[selected_df['Q'] != selected_df['T']]
        match_row = selected_df.loc[selected_df['Similarity(%)'].idxmax()]
        highest_similarity = match_row['Similarity(%)']
        number_of_equals = (selected_df['Similarity(%)'].values == highest_similarity).sum()
        return number_of_equals

    curated_df['Closest DB Match'] = curated_df['Query'].apply(closest_match_to_query)
    curated_df['Highest Similarity(%)'] = curated_df['Query'].apply(highest_similarity)
    curated_df['Alternative Matches'] = curated_df['Query'].apply(number_of_equal_highest_similarities) - 1
    curated_df['M'] = curated_df['Closest DB Match'].str[4:].astype(int)
    curated_df = curated_df.iloc[:, [0, 1, 2, 3, 6, 7, 8, 4, 5, 9]]
    curated_df.to_csv(curated_vsearch_output, sep=',', index=False)

def remove_duplicates(df):
    df = df.drop_duplicates(subset=['Q'])
    df = df.reset_index(drop=True)
    return df

def rearrange_vsearch_output(df):
    df['temp_index'] = df['Query'].str.replace("Query_", "").astype(int)
    df = df.set_index('temp_index')
    selected_columns = ['Q', 'M', 'Highest Similarity(%)', 'Alternative Matches']
    df = df[selected_columns]
    return df

def rearrange_query(df):
    df.insert(0, 'temp_index', df['Query#'])
    df = df.set_index('temp_index')
    return df

def combine_query_and_vsearch_dfs(df1, df2, df3):
    output_df = pd.concat([df1,df2], axis=1, join='outer')
    database_df = df3
    output_df['Vs Name'] = output_df['M'].map(database_df.set_index('DB#')['DB Name'])
    output_df['Vs ID'] = output_df['M'].map(database_df.set_index('DB#')['DB Accession Number'])
    output_df = output_df.replace(np.nan, 'N/A', regex=True)
    output_df = output_df.rename(columns={'M': 'Vs DB#', 'Highest Similarity(%)': 'Similarity(%)'})
    retained_columns = ['Query#', 'Hun#', 'Fas#', 'DB Name', 'Vs DB#', 'Vs Name', 'Similarity(%)', 'Alternative Matches',
                 'DB Accession Number', 'DB Found', '16S rRNA sequence']
    output_df = output_df[retained_columns]
    return output_df

def create_output_file(query_df, database_df, curated_vsearch_output, output_file):
    edited_vsearch_output_df = read_csv(curated_vsearch_output)
    edited_vsearch_output_df = remove_duplicates(edited_vsearch_output_df)
    edited_vsearch_output_df = rearrange_vsearch_output(edited_vsearch_output_df)
    query_df = rearrange_query(query_df)
    output_df = combine_query_and_vsearch_dfs(query_df, edited_vsearch_output_df, database_df)
    write_csv(output_df, output_file)

def tidy_files(curated_vsearch_output, query_fastaname, database_fastaname, vsearch_output):
    os.mkdir('data')
    shutil.move(curated_vsearch_output, 'data')
    shutil.move(vsearch_output, 'data')
    shutil.move(query_fastaname, 'data')
    shutil.move(database_fastaname, 'data')

def summary_statement(output_directory, start_time):
    program_length = datetime.now() - start_time
    print(f'{sys.argv[0]} has completed in {program_length}. {sys.argv[1]} was processed and the Summary.csv file can be '
          f'found in the new {output_directory} directory.')

main()
