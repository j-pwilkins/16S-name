import pandas as pd
import numpy as np
from datetime import datetime

def main():
    start_time = datetime.now()
    shared_matches_df = produce_shared_matches_df()

    summary_statement(start_time)

# produce a df containing all the possible alternative mapping options for each query
def produce_shared_matches_df():
    curated_vsearch_output, database_df, query_df, summary_df, output_columns = read_inputs()
    query_plus_vsearch_output_df = merge_query_with_vsearch_output(curated_vsearch_output, query_df)
    shared_matches_df = select_only_shared_matches(query_plus_vsearch_output_df)
    shared_matches_df = create_selected_column(shared_matches_df, summary_df)
    shared_matches_df = add_database_detail(shared_matches_df, database_df, output_columns)
    shared_matches_production_statement(summary_df, shared_matches_df)
    return shared_matches_df

## functions for produce_shared_matches_df
# read the inputs produced from mapping - these will be passed to the function when integrating into the the new mapping.py
def read_inputs():
    curated_vsearch_output = read_csv('Query_vs_All.csv')
    database_df = read_csv('db_497.csv')
    query_df = read_csv('V1V3.csv')
    summary_df = read_csv('Summary_V1V3_vs_db_497.csv')
    output_columns = ['Query#', 'Hun#', 'Fas#', 'Accession Number', 'Name', 'Shared Matches', 'Selected', 'Similarity(%)',
                    'Vs DB#',
                    'DB Name', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
                    'DB Accession Number', '16S rRNA sequence']
    return curated_vsearch_output, database_df, query_df, summary_df, output_columns

# adds input sequence detail to the curated vsearch output df
def merge_query_with_vsearch_output(vsearch_df, query_df):
    query_numbers = list(range(1, len(query_df) + 1))
    query_df.insert(0, 'Query#', query_numbers)
    vsearch_plus_df = vsearch_df.merge(query_df, left_on='Q', right_on='Query#')
    vsearch_plus_df = vsearch_plus_df.rename(columns={'DB Name': 'Name', 'DB#': 'Vs DB#', 'DB Accession Number': 'Accession Number', 'Alternative Matches': 'Shared Matches'})
    return vsearch_plus_df

# large, upgraded curated vsearch output vsearch_plus_df is trimmed down to only include those queries with shared matches
def select_only_shared_matches(vsearch_plus_df):
    vsearch_plus_df = vsearch_plus_df.reset_index(drop=True)
    vsearch_plus_df = vsearch_plus_df[vsearch_plus_df['Q'] != vsearch_plus_df['T']]
    vsearch_plus_df = vsearch_plus_df.loc[vsearch_plus_df['Shared Matches'] > 0, :]
    vsearch_plus_df['max_similarity'] = vsearch_plus_df.groupby('Query')['Similarity(%)'].transform('max')
    vsearch_plus_df = vsearch_plus_df.loc[vsearch_plus_df['Similarity(%)'] == vsearch_plus_df['max_similarity'], :]
    vsearch_plus_df = vsearch_plus_df[['Query#', 'Hun#', 'Fas#', 'Accession Number', 'Name', 'Similarity(%)', 'Shared Matches', 'Vs DB#',
             '16S rRNA sequence']]
    vsearch_plus_df['Shared Matches'] = vsearch_plus_df['Shared Matches'] + 1
    return vsearch_plus_df

# A column is added, showing which match was selected by mapping program, and which are additional matches
def create_selected_column(shared_matches_df, summary_df):
    shared_matches_df['DB#'] = shared_matches_df['Vs DB#'].str.replace('DB#_', '').astype(int)
    shared_matches_df['Selected'] = 'N'
    for index, row in shared_matches_df.iterrows():
        query_number = row['Query#']
        db_number = row['DB#']
        if summary_df[
            (summary_df['Query#'].isin([query_number])) & (summary_df['Vs DB#'].isin([db_number]))].any().any():
            shared_matches_df.loc[index, 'Selected'] = 'Y'
    return shared_matches_df

# Adds detail to the shared matches df from the database of full length 16S rRNA gene sequences
def add_database_detail(shared_matches_df, database_df, output_columns):
    database_df = database_df.drop(columns=['Hun#', 'Fas#', 'Tree#', 'Vs DB#', 'Vs Name', 'Similarity (%)', 'Alternative Matches', 'Vs ID', 'DB Found', '16S rRNA sequence'])
    shared_matches_df = shared_matches_df.merge(database_df, on='DB#')
    shared_matches_df = shared_matches_df.sort_values(by=['Query#', 'DB#'], ascending=True)     # this line might be temp to help coding
    shared_matches_df = shared_matches_df.reindex(columns=output_columns)
    return shared_matches_df

# output statement for shared_matches_df
def shared_matches_production_statement(summary_df, shared_matches_df):
    summary_queries = count_unique_values(summary_df, 'Query#')
    shared_matches_queries = count_unique_values(shared_matches_df, 'Query#')
    print(f"{shared_matches_queries} of the {summary_queries} queries have shared matches.")

## finishing function
def summary_statement(start_time):
    program_length = datetime.now() - start_time
    print(f'Program finished in {program_length}')

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

def count_unique_values(df, column_header):
    unique_values = df[column_header].unique().tolist()
    number_of_unique_values = len(unique_values)
    return number_of_unique_values


main()