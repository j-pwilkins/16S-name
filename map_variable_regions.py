import pandas as pd
import numpy as np
import os
import shutil
from datetime import datetime

def main():
    start_time = datetime.now()
    input_sequences_list = 'test1.csv'
    input_database = 'Summary_test10.csv'
    vsearch_output = 'bob.csv'
    curated_vsearch_output = 'Query_vs_All.csv'
    insert = '_vs_'
    output_directory = input_sequences_list.rstrip('.csv') + insert + input_database.rstrip('.csv')
    output_file = 'Summary_' + output_directory + '.csv'
    if os.path.exists(output_directory):    # temp lines for testing
        shutil.rmtree(output_directory)     # temp
        print(f'The previous {output_directory} folder has been deleted')
    os.makedirs(output_directory, exist_ok=True)
    shutil.copy(input_sequences_list, output_directory)
    shutil.copy(input_database, output_directory)
    shutil.copy(vsearch_output, output_directory)   #dummy line to make sure dummy vsearch output is in working directory
    os.chdir(output_directory)
    import_csv(input_sequences_list, input_database, output_file)
    query_df, database_df, length_of_query, length_of_database, query_fastaname, database_fastaname = import_csv(input_sequences_list, input_database, output_file)
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

def import_csv(input_sequences_list, input_database, output_file):
    query_df = pd.read_csv(input_sequences_list)
    database_df = pd.read_csv(input_database)
    length_of_query = len(query_df.index.tolist())
    length_of_database = len(database_df.index.tolist())
    query_base_name = os.path.splitext(input_sequences_list)[0]
    query_fastaname = query_base_name + '.fa'
    database_base_name = os.path.splitext(input_database)[0]
    database_fastaname = database_base_name + '.fa'
    query_length_as_list = list(range(length_of_query))
    input_list = [x + 1 for x in query_length_as_list]
    query_df.insert(0, 'Query#', input_list)
    # query_df.to_csv(output_file, sep=',', index=False)
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

def xlookup(lookup_value, lookup_array, return_array, if_not_found: str = ''):
    match_value = return_array.loc[lookup_array == lookup_value]
    if match_value.empty:
        return f'"{lookup_value}" not found!' if if_not_found == '' else if_not_found

    else:
        return match_value.tolist()[0]

def run_vsearch(query_fastaname, database_fastaname, vsearch_output):
    cmd = 'vsearch --usearch_global ' + query_fastaname + ' --db ' + database_fastaname + ' --blast6out ' + vsearch_output + ' --acceptall --id 0.0 --maxaccepts 0'
    # os.system(cmd)
    # print(cmd)
    print ('A pre-prepared vsearch output has been used here')

def parse_vsearch_usearchglobal_blast6(vsearch_output, curated_vsearch_output):
    query_vs_targetdb_df = pd.read_csv(vsearch_output, usecols=[0, 1, 2], sep='\t',
                                            names=["Query#", "DB#", "Similarity(%)"])
    df_length = len(query_vs_targetdb_df)
    vsearch_output_order = pd.Series(np.arange(1, df_length + 1, 1))
    query_vs_targetdb_df.insert(0, 'Vsearch Order', vsearch_output_order, True)
    query_vs_targetdb_df['Q'] = query_vs_targetdb_df['Query#'].str[6:].astype(int)
    query_vs_targetdb_df['T'] = query_vs_targetdb_df['DB#'].str[4:].astype(int)
    query_vs_targetdb_df = query_vs_targetdb_df.sort_values(['Q', 'T'])
    query_vs_targetdb_df.to_csv(curated_vsearch_output, sep=',', index=False)

def add_similarity_columns(curated_vsearch_output):
    curated_df = pd.read_csv(curated_vsearch_output)

    def closest_match_to_query(input):
        selected_df = curated_df[:]
        selected_df = selected_df.loc[selected_df['Query#'] == input]
        selected_df = selected_df.loc[selected_df['Q'] != selected_df['T']]
        match_row = selected_df.loc[selected_df['Similarity(%)'].idxmax()]
        closest_match = match_row['DB#']
        return closest_match

    def highest_similarity(input):
        selected_df = curated_df[:]
        selected_df = selected_df.loc[selected_df['Query#'] == input]
        selected_df = selected_df.loc[selected_df['Q'] != selected_df['T']]
        match_row = selected_df.loc[selected_df['Similarity(%)'].idxmax()]
        highest_similarity = match_row['Similarity(%)']
        return highest_similarity

    def number_of_equal_highest_similarities(input):
        selected_df = curated_df[:]
        selected_df = selected_df.loc[selected_df['Query#'] == input]
        selected_df = selected_df.loc[selected_df['Q'] != selected_df['T']]
        match_row = selected_df.loc[selected_df['Similarity(%)'].idxmax()]
        highest_similarity = match_row['Similarity(%)']
        number_of_equals = (selected_df['Similarity(%)'].values == highest_similarity).sum()
        return number_of_equals

    curated_df['Closest DB Match'] = curated_df['Query#'].apply(closest_match_to_query)
    curated_df['Highest Similarity(%)'] = curated_df['Query#'].apply(highest_similarity)
    curated_df['Alternative Matches'] = curated_df['Query#'].apply(number_of_equal_highest_similarities) - 1
    curated_df['M'] = curated_df['Closest DB Match'].str[4:].astype(int)
    curated_df = curated_df.iloc[:, [0, 1, 2, 3, 6, 7, 8, 4, 5, 9]]
    curated_df.to_csv(curated_vsearch_output, sep=',', index=False)

def create_output_file(query_df, database_df, curated_vsearch_output, output_file):
    curated_df = pd.read_csv(curated_vsearch_output)
    curated_df = curated_df.drop_duplicates(subset=['Q'])
    curated_df = curated_df.reset_index(drop=True)
    columns_to_be_inserted = ['Q', 'M', 'Highest Similarity(%)', 'Alternative Matches']
    insert = curated_df.loc[:, columns_to_be_inserted]
    query_df = pd.concat([query_df, insert], 1)
    query_df['Vs Name'] = query_df['M'].map(database_df.set_index('DB#')['DB Name'])
    query_df['Vs ID'] = query_df['M'].map(database_df.set_index('DB#')['DB Accession Number'])
    query_df = query_df.loc [:, ['Query#', 'Hun#', 'Fas#', 'DB Name', 'M', 'Vs Name', 'Highest Similarity(%)', 'Alternative Matches', 'DB Accession Number', 'Vs ID', 'DB Found', '16S rRNA sequence']]
    query_df = query_df.rename(columns={'M': 'Vs DB#', 'Highest Similarity(%)': 'Similarity(%)'})
    query_df.to_csv(output_file, sep=',', index=False)

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
