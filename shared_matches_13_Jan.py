import pandas as pd
import numpy as np
from datetime import datetime

def main():
    start_time = datetime.now()
    curated_vsearch_output = read_csv('Query_vs_All.csv')
    database_df = read_csv('db_497.csv')
    query_df = read_csv('V4.csv')
    summary_df = read_csv('Summary_V4_vs_db_497.csv')
    shared_matches_df = produce_shared_matches_csv(curated_vsearch_output, database_df, query_df, summary_df)
    # shared_matches_df = read_csv('shared_matches.csv')
    add_consensus_sequences(shared_matches_df)
    summary_statement(start_time)

def produce_shared_matches_csv(curated_vsearch_output, database_df, query_df, summary_df):
    query_plus_vsearch_output_df = merge_query_with_vsearch_output(curated_vsearch_output, query_df)
    shared_matches_df = select_only_shared_matches(query_plus_vsearch_output_df)
    shared_matches_df = create_selected_column(shared_matches_df, summary_df)
    merged_df = merge_query_with_database(shared_matches_df, database_df)
    write_csv(merged_df, 'shared_matches.csv')
    return merged_df

def add_consensus_sequences(shared_matches_df):
    strict_df = add_strict_consensus_sequences(shared_matches_df)
    plurality_df = add_plurality_consensus_sequences(shared_matches_df)
    write_csv(strict_df, 'strict.csv')
    write_csv(plurality_df, 'plurality.csv')
    write_csv(shared_matches_df, 'shared_matches_df.csv')
    shared_matches_with_consensus_df = join_dfs(strict_df, plurality_df, shared_matches_df)
    write_csv(shared_matches_with_consensus_df, 'shared_matches_with_consensus.csv')

def add_plurality_consensus_sequences(shared_matches_df):
    plurality_df = add_plurality_column_a(shared_matches_df)
    plurality_df = add_plurality_consensus_sequences_b_to_n(plurality_df)
    # write_csv(plurality_df, 'Quer#.csv')
    plurality_df = process_df(plurality_df)
    return plurality_df

def add_plurality_column_a(shared_matches_df):
    # Create a copy of the dataframe
    df = shared_matches_df.copy()
    # group the rows in the dataframe by the 'Query#' column
    grouped_df = df.groupby('Query#')
    # get the counts of each unique value in the 'A' column for each group
    consensus_values = grouped_df['A'].value_counts()
    # iterate through each row in the dataframe, and assign a value to the 'a' column based on the conditions
    df['a'] = df.apply(lambda x:
        # for modal 'A' frequencies - if the modal frequency in 'A' is NOT unique
        '-' if (consensus_values.loc[x['Query#']] == consensus_values.loc[x['Query#']].max()).sum() > 1
        # for non-modal 'A' frequencies
        else 'x' if x['A'] != consensus_values.loc[x['Query#']].idxmax()
        #for modal 'A' frequencies
        else x['A'], axis=1)
    return df

def add_plurality_consensus_sequences_b_to_n(shared_matches_df):
    def update_column(col_name, new_col_name):
        grouped_df = shared_matches_df.groupby('Query#')
        consensus_values = grouped_df[col_name].value_counts()
        shared_matches_df[new_col_name] = shared_matches_df.apply(lambda x: '-' if x['a'] == '-' or (consensus_values.loc[x['Query#']] == consensus_values.loc[x['Query#']].max()).sum() > 1
                                                              else 'x' if x[col_name] != consensus_values.loc[x['Query#']].idxmax()
                                                              else consensus_values.loc[x['Query#']].idxmax(), axis=1)
        shared_matches_df.drop(shared_matches_df[shared_matches_df[new_col_name] == 'x'].index, inplace=True)
    # add columns 'b' - 'n'
    update_column('B', 'b')
    update_column('C', 'c')
    update_column('D', 'd')
    update_column('E', 'e')
    update_column('F', 'f')
    update_column('G', 'g')
    update_column('H', 'h')
    update_column('I', 'i')
    update_column('J', 'j')
    update_column('K', 'k')
    update_column('L', 'l')
    update_column('M', 'm')
    update_column('N', 'n')
    return shared_matches_df

def check_barcode_consistency(df, column):
    # Get unique values in specified column
    unique_query_nums = df[column].unique()
    # Iterate through each unique value
    for query_num in unique_query_nums:
        # Get all rows with the current unique value
        query_rows = df[df[column] == query_num]
        # Check if entries in column 'a' are consistent
        if query_rows['a'].nunique() > 1:
            print('The unique values in the df are NOT consistent')
            return False
    return True

def process_df(plurality_df):
    if check_barcode_consistency(plurality_df, 'Query#'):
        print('The unique values in the shared_and_plurality_df are consistent')
        # Drop old barcode columns & rename new barcode columns
        df_grouped = plurality_df.drop(columns=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N'])
        df_grouped = df_grouped.rename(
            columns={'a': 'A', 'b': 'B', 'c': 'C', 'd': 'D', 'e': 'E', 'f': 'F', 'g': 'G', 'h': 'H', 'i': 'I', 'j': 'J',
                     'k': 'K', 'l': 'L', 'm': 'M', 'n': 'N'})
        # Group the rows by the unique values of the 'Query#' column
        df_grouped = df_grouped.groupby('Query#', as_index=False).first()
        df_grouped['Selected'] = 'Cp'
        df_grouped['Similarity(%)'] = 'N/A'
        df_grouped['Vs DB#'] = 'N/A'
        df_grouped['DB Name'] = 'N/A'
        df_grouped['DB Accession Number'] = 'N/A'
        # df_grouped['Similarity(%)', 'Vs DB#', 'DB Name', 'DB Accession Number'] = 'N/A'
        df_grouped = df_grouped.reindex(
            columns=['Query#', 'Hun#', 'Fas#', 'Accession Number', 'Name', 'Shared Matches', 'Selected',
                     'Similarity(%)', 'Vs DB#', 'DB Name', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L',
                     'M', 'N', 'DB Accession Number', '16S rRNA sequence'])
        return df_grouped
    else:
        return

def add_strict_consensus_sequences(shared_matches_df):
    consensus_df = produce_unique_query_values_column(shared_matches_df)
    add_strict_consensus_barcodes(consensus_df, shared_matches_df)
    # write_csv(consensus_df, 'strict.csv')
    consensus_df = add_detail_for_merge(consensus_df, shared_matches_df)
    consensus_df = complete_df_for_merge(consensus_df)
    column_headers = list(consensus_df. columns)
    print(column_headers)
    # write_csv(consensus_df, 'premerge.csv')
    return consensus_df


def merge_query_with_vsearch_output(df1, df2):
    query_numbers = list(range(1, len(df2) + 1))
    df2.insert(0, 'Query#', query_numbers)
    df3 = df1.merge(df2, left_on='Q', right_on='Query#')
    df3 = df3.rename(columns={'DB Name': 'Name', 'DB#': 'Vs DB#', 'DB Accession Number': 'Accession Number', 'Alternative Matches': 'Shared Matches'})
    return df3

def select_only_shared_matches(df):
    df = df.reset_index(drop=True)
    df = df[df['Q'] != df['T']]
    df = df.loc[df['Shared Matches'] > 0, :]
    df['max_similarity'] = df.groupby('Query')['Similarity(%)'].transform('max')
    df = df.loc[df['Similarity(%)'] == df['max_similarity'], :]
    df = df[['Query#', 'Hun#', 'Fas#', 'Accession Number', 'Name', 'Similarity(%)', 'Shared Matches', 'Vs DB#',
             '16S rRNA sequence']]
    df['Shared Matches'] = df['Shared Matches'] + 1
    return df

def create_selected_column(reduced_df, summary_df):
    reduced_df['DB#'] = reduced_df['Vs DB#'].str.replace('DB#_', '').astype(int)
    reduced_df['Selected'] = 'N'
    for index, row in reduced_df.iterrows():
        query_number = row['Query#']
        db_number = row['DB#']
        if summary_df[
            (summary_df['Query#'].isin([query_number])) & (summary_df['Vs DB#'].isin([db_number]))].any().any():
            reduced_df.loc[index, 'Selected'] = 'Y'
    return reduced_df

def merge_query_with_database(reduced_df, database_df):
    database_df = database_df.drop(columns=['Hun#', 'Fas#', 'Tree#', 'Vs DB#', 'Vs Name', 'Similarity (%)', 'Alternative Matches', 'Vs ID', 'DB Found', '16S rRNA sequence'])
    merged_df = reduced_df.merge(database_df, on='DB#')
    merged_df = merged_df.sort_values(by=['Query#', 'DB#'], ascending=True)     # this line might be temp to help coding
    column_order = ['Query#', 'Hun#', 'Fas#', 'Accession Number', 'Name', 'Shared Matches', 'Selected', 'Similarity(%)', 'Vs DB#',
                    'DB Name', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'DB Accession Number', '16S rRNA sequence']
    merged_df = merged_df.reindex(columns=column_order)
    return merged_df

def produce_unique_query_values_column(merged_df):
    unique_values = merged_df['Query#'].unique()
    df = pd.DataFrame(unique_values, columns=['Query#'])
    df['Query#'] = df['Query#'].sort_values()
    return df

def add_strict_consensus_barcodes(new_df, output_df):
    def consensus_barcode_a(query):
        rows = output_df[output_df['Query#'] == query]
        values = rows['A'].values
        if len(set(values)) == 1:
            return values[0]
        return '-'

    def consensus_barcode(query, curr_column, prev_column):
        row = new_df[new_df['Query#'] == query]
        curr_value = row[prev_column].iloc[0]
        if curr_value == '-':
            return '-'
        rows = output_df[output_df['Query#'] == query]
        values = rows[curr_column].values
        if len(set(values)) == 1:
            return values[0]
        return '-'

    new_df['A'] = new_df['Query#'].apply(consensus_barcode_a)

    for curr_column, prev_column in zip(['B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N'], ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']):
        new_df[curr_column] = new_df['Query#'].apply(lambda query: consensus_barcode(query, curr_column, prev_column))

def add_detail_for_merge(df1, df2):
    df2 = df2.drop_duplicates(subset='Query#')
    df2 = df2[['Query#', 'Hun#', 'Fas#', 'Accession Number', 'Name', 'Shared Matches', '16S rRNA sequence']]
    df3 = df1.merge(df2, on='Query#')
    return df3

def complete_df_for_merge(df):
    # Reorder columns in the combined_df dataframe
    complete_df = df.reindex(
        columns=['Query#', 'Hun#', 'Fas#', 'Accession Number', 'Name', 'Shared Matches', 'A', 'B',
                 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', '16S rRNA sequence'])

    # Insert new columns with constant values
    complete_df.insert(6, 'Selected', 'Cs')
    complete_df.insert(7, 'Similarity(%)', 'N/A')
    complete_df.insert(8, 'Vs DB#', 'N/A')
    complete_df.insert(9, 'DB Name', 'N/A')
    complete_df.insert(24, 'DB Accession Number', 'N/A')
    return complete_df

def join_dfs(df1, df2, df3):
    df = pd.concat([df1, df2, df3])
    df.sort_values(by=['Query#', 'Selected'], ascending=[True, False], inplace=True)
    return df


def summary_statement(start_time):
    program_length = datetime.now() - start_time
    print(f'Program finished in {program_length}')

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

main()