import pandas as pd
import numpy as np
from datetime import datetime

def main():
    start_time = datetime.now()
    shared_matches_df = produce_shared_matches_df()
    add_consensus_sequences(shared_matches_df)
    summary_statement(start_time)

# Level 1 # first shared_matches sub-program
# produce a df containing all the possible alternative mapping options for each query
def produce_shared_matches_df():
    curated_vsearch_output, database_df, query_df, summary_df, output_columns = read_inputs()
    query_plus_vsearch_output_df = merge_query_with_vsearch_output(curated_vsearch_output, query_df)
    shared_matches_df = select_only_shared_matches(query_plus_vsearch_output_df)
    shared_matches_df = create_selected_column(shared_matches_df, summary_df)
    shared_matches_df = add_database_detail(shared_matches_df, database_df, output_columns)
    shared_matches_production_statement(summary_df, shared_matches_df)
    return shared_matches_df

# Level 1 # second shared_matches sub-program
# create three (strict, majority, plurality) consensus dfs, then join them to shared_matches
def add_consensus_sequences(shared_matches_df):
    pre_consensus_df = shared_matches_df.copy()
    barcode_columns = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N']
    majority_consensus_df = create_majority_consensus_df(pre_consensus_df, shared_matches_df, barcode_columns, minimum_threshold=0.5, selected_filter='Cm')
    shared_matches_with_consensus_df = concatenate_dfs(shared_matches_df, majority_consensus_df)
    shared_matches_with_consensus_df = shared_matches_with_consensus_df.sort_values(by='Query#')
    write_csv(shared_matches_with_consensus_df, 'smwc.csv')

# Level 2 # create majority consensus df
def create_majority_consensus_df(pre_consensus_df, shared_matches_df, barcode_columns, minimum_threshold, selected_filter):
    pre_consensus_df, consensus_df = add_majority_barcodes(pre_consensus_df, barcode_columns, minimum_threshold)
    pre_consensus_df, matching_df = additional_consensus_barcode_filtering(pre_consensus_df, barcode_columns)
    consensus_df = concatenate_dfs(consensus_df, matching_df)
    consensus_df = fill_consensus_df_missing_values(pre_consensus_df, consensus_df, shared_matches_df, barcode_columns)
    check_consensus_df_correct (consensus_df, shared_matches_df)
    consensus_df = prepare_consensus_df_for_merge_with_shared_matches_df(consensus_df, barcode_columns, selected_filter)
    return consensus_df

# Level 3 # filters shared matches and adds majority consensus barcodes. Returns processed queries in consensus_df, and unprocessed queries in pre_consensus_df
def add_majority_barcodes(df, columns, minimum_threshold):
    consensus_df = pd.DataFrame(columns=df.columns)     # initialise empty df for transfer of processed queries
    grouped = df.groupby('Query#')
    df_filtered = df.copy()
    # iterate over each column
    for column in columns:
        new_column = column.lower()
        df[new_column] = '_'        # new barcode column, to be created from old
        # iterate over each group
        for name, group in grouped:
            # calculate modal_value, frequency_of_modal_value and group_size variables
            modal_value = group[column].mode()[0]
            frequency_of_modal_value = group[column].value_counts()[modal_value]
            group_size = len(group)
            # calculate threshold for group in this column
            threshold, test_result = calculate_threshold(group_size, frequency_of_modal_value, minimum_threshold)
            # if group > threshold, then keep those rows with modal value and drop rest
            if test_result:
                df.loc[(df['Query#'] == name) & (df[column] == modal_value), new_column] = modal_value
                df.loc[(df['Query#'] == name) & (df[column] != modal_value), new_column] = 'x'
                df = df[df[new_column] != 'x']
                df_filtered = df[df[new_column] != 'x']
                grouped = df_filtered.groupby('Query#')  # this line re-group the dataset for the next column iteration, delete?
            # if group < threshold, copy first_row to consensus_df. Then drop rows
            if not test_result:
                first_row = group.head(1).copy()
                first_row.loc[:, new_column] = '-'
                consensus_df = pd.concat([consensus_df, first_row])
                df = df[df['Query#'] != name]
        grouped = df_filtered.groupby('Query#')
    return df, consensus_df

# Level 4 # calculates threshold for each 'Query#' group in each barcode column. Tests if group passes threshold
def calculate_threshold(group_size, frequency_of_modal_value, minimum_threshold):
    threshold = group_size * minimum_threshold
    test_result = frequency_of_modal_value / group_size > minimum_threshold
    return threshold, test_result

# Level 3 # where query does not end with unique barcode, this function checks all rows for each group match, then combines into one row. Combined rows are returned in matching_df for concatenation with consenus_df. pre_consensus_df is returned to check is empty
def additional_consensus_barcode_filtering(pre_consensus_df, columns):
    new_columns = [x.lower() for x in columns]
    matching_df = pd.DataFrame(columns=pre_consensus_df.columns)
    groups = pre_consensus_df.groupby('Query#')
    match_counter = 0
    for group_name, group_df in groups:
        first_row = group_df.iloc[0]
        match = group_df[new_columns] == first_row[new_columns].values
        if match.all().all():
            matching_df.loc[match_counter] = first_row
            match_counter += 1
            pre_consensus_df = pre_consensus_df[pre_consensus_df['Query#'] != group_name]
    return pre_consensus_df, matching_df

# L3 # completes filtering of shared matches df. Output
def fill_consensus_df_missing_values(pre_consensus_df, consensus_df, shared_matches_df, barcode_columns):
    if check_queries_complete_one(pre_consensus_df):
        check_queries_complete_two(consensus_df, shared_matches_df)
        consensus_df = remove_duplicated_processed_queries(consensus_df)
        consensus_df, new_barcode_columns = add_missing_columns(consensus_df, barcode_columns)
        consensus_df = replace_nans(consensus_df, new_barcode_columns)
        return consensus_df
    else:
        print('stopping at process_consensus_df')

# L4 # check that all queries have been processed
def check_queries_complete_one(pre_consensus_df):
    if pre_consensus_df.empty:
        return True
    else:
        unprocessed_groups = pre_consensus_df['Query#'].nunique()
        print(f'There are {unprocessed_groups} unprocessed groups still to be processed')
        return False

# L4 # Check consensus df queries match input df queries
def check_queries_complete_two(consensus_df, shared_matches_df):
    consensus_queries = consensus_df['Query#'].unique()
    shared_queries = shared_matches_df['Query#'].unique()

    missing_queries_list = list(set(shared_queries) - set(consensus_queries))
    if missing_queries_list:
        print(f"Some queries have not been processed, they are {missing_queries_list}")

    bogus_queries_list = list(set(consensus_queries) - set(shared_queries))
    if bogus_queries_list:
        print(
            f"There are queries that have been processed that weren't in the shared_matches_df, they are {bogus_queries_list}")

#L4 # This is an extra function added because of bug in testing where a query could be processed twice, in the immediately following barcode column
# It checks for duplicated queries in the consensus_df, keeping only the first entry. Statement printed of any deletions.
def remove_duplicated_processed_queries(consensus_df):
    unique_df = consensus_df.drop_duplicates(subset='Query#')   # groups by unique Query# values
    duplicated_df = consensus_df[~consensus_df.index.isin(unique_df.index)]  # creates df of duplicated values
    if duplicated_df.empty == False:                            # checks if empty
        rows_deleted = len(duplicated_df)
        duplicated_queries = list(set(duplicated_df['Query#']))
        print(f"{rows_deleted} rows have been deleted from the consensus_df. The following queries were duplicated: {duplicated_queries}.")
    check_unique_df = list(set(consensus_df['Query#'])) == list(set(unique_df['Query#']))
    if check_unique_df:
        return unique_df

# L4 # columns are only added to consensus_df when they are used. This function adds any extra columns (a-n) that are missing
def add_missing_columns(consensus_df, barcode_columns):
    new_barcode_columns = [x.lower() for x in barcode_columns]
    consensus_columns = list(consensus_df.columns.values)
    if consensus_columns != new_barcode_columns:
        missing_columns = (set(new_barcode_columns).difference(consensus_columns))
        for str in missing_columns:
            consensus_df = consensus_df.assign(**{str: '-'})
    return consensus_df, new_barcode_columns

#L4 # When consensus barcode is added it stops when threshold failed. This function extends barcode ('-' values) to n column
def replace_nans(consensus_df, new_barcode_columns):
    consensus_df[new_barcode_columns] = consensus_df[new_barcode_columns].fillna('-')
    return consensus_df

#L3 # Checks consensus df to make sure that queries from input and output dfs match, and that output queries are not duplicated
def check_consensus_df_correct (consensus_df, shared_matches_df):
    consensus_df.name = 'consensus_df'
    shared_matches_df.name = 'shared_matches_df'
    check_unique_values(consensus_df, shared_matches_df, 'Query#')
    check_df_column_unique_entries_only(consensus_df, 'Query#')

#L3 # prepares consensus_df for merging with shared_matches_df
def prepare_consensus_df_for_merge_with_shared_matches_df(consensus_df, barcode_columns, selected_filter):
    consensus_df = consensus_df.drop(columns=barcode_columns)       # drops old barcode columns
    consensus_df = consensus_df.rename(
        columns={'a': 'A', 'b': 'B', 'c': 'C', 'd': 'D', 'e': 'E', 'f': 'F', 'g': 'G', 'h': 'H', 'i': 'I', 'j': 'J',
                 'k': 'K', 'l': 'L', 'm': 'M', 'n': 'N'})   # changes new barcode columns names to old
    consensus_df = consensus_df.sort_values(by='Query#')
    consensus_df['Selected'] = selected_filter
    consensus_df['Similarity(%)'] = 'N/A'
    consensus_df['Vs DB#'] = 'N/A'
    consensus_df['DB Name'] = 'N/A'
    consensus_df['DB Accession Number'] = 'N/A'
    curated_vsearch_output, database_df, query_df, summary_df, output_columns = read_inputs()
    if set(output_columns).issubset(set(consensus_df.columns)):
        consensus_df = consensus_df.reindex(output_columns, axis=1)
    else:
        print("Error: Not all columns in 'output_columns' list are present in 'consensus_df' DataFrame.")
    return consensus_df

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
def check_df_column_unique_entries_only(df, column_header):
    if df[column_header].duplicated().any():
        duplicated_rows = df[column_header].duplicated().sum()
        print(f"{duplicated_rows} queries in the {df.name} '{column_header}' column are duplicated")
    else:
        print(f"All entries in the {df.name} '{column_header}' column are unique")

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

def concatenate_dfs(*dfs):
    column_headers = [list(df.columns) for df in dfs]
    if all(headers == column_headers[0] for headers in column_headers):
        concatenated_df = pd.concat(dfs)
        return concatenated_df
    else:
        df_names = ', '.join([str(df) for df in dfs])
        print(f"{df_names} could not be concatenated as the column headers do not match. Please rectify")

def count_unique_values(df, column_header):
    unique_values = df[column_header].unique().tolist()
    number_of_unique_values = len(unique_values)
    return number_of_unique_values

def read_csv(filename):
    try:
        df = pd.read_csv(filename)
        return df
    except FileNotFoundError:
        print("The file '{}' could not be found. Make sure the file is in the correct location and try again.".format(
            filename))
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