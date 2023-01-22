import pandas as pd
import numpy as np
import os
import shutil
import sys
from datetime import datetime

def main():
    shared_matches_df = read_csv('shared_matches.csv')
    add_consensus_sequences(shared_matches_df)

def add_consensus_sequences(shared_matches_df):
    consensus_df = add_majority_consensus_sequences(shared_matches_df)
    write_csv(consensus_df, 'consensus.csv')

def add_majority_consensus_sequences(shared_matches_df):
    pre_consensus_df = shared_matches_df.copy()
    columns = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N']
    pre_consensus_df, consensus_df = add_majority_columns(shared_matches_df, columns)
    write_csv(pre_consensus_df, 'majority.csv')
    write_csv(consensus_df, 'rejected.csv')
    consensus_df = process_consensus_df(pre_consensus_df, consensus_df, shared_matches_df, columns)
    consensus_df.name = 'consensus_df'
    shared_matches_df.name = 'shared_matches_df'
    check_unique_values(consensus_df, shared_matches_df, 'Query#')
    check_df_column_unique_entries_only(consensus_df, 'Query#')
    return consensus_df

def calculate_threshold(group_size, frequency_of_modal_value, minimum_threshold=0.5):
    threshold = group_size * minimum_threshold
    test_result = frequency_of_modal_value / group_size > minimum_threshold
    return threshold, test_result

def add_majority_columns(df, columns):
    # initialize an empty dataframe to store rejected rows
    consensus_df = pd.DataFrame(columns=df.columns)
    # group by 'Query#'
    grouped = df.groupby('Query#')
    df_filtered = df.copy()
    # iterate over each column
    for column in columns:
        # create a new column 'a' with default value '_'
        new_column = column.lower()
        df[new_column] = '_'
        # iterate over each group
        for name, group in grouped:
            # calculate modal_value, frequency_of_modal_value and group_size
            modal_value = group[column].mode()[0]
            # print(group[column])
            frequency_of_modal_value = group[column].value_counts()[modal_value]
            group_size = len(group)
            # calculate threshold
            threshold, test_result = calculate_threshold(group_size, frequency_of_modal_value)
            group_name = column + str(name)
            if test_result:
                df.loc[(df['Query#'] == name) & (df[column] == modal_value), new_column] = modal_value
                df.loc[(df['Query#'] == name) & (df[column] != modal_value), new_column] = 'x'
                df = df[df[new_column] != 'x']
                # print(f"{group_name} - {test_result} . {group_size} rows. Frequency - {frequency_of_modal_value}. Threshold - {threshold}. Mode - {modal_value}.")
                df_filtered = df[df[new_column] != 'x']
                grouped = df_filtered.groupby('Query#')  # this line re-group the dataset for the next column iteration
            if not test_result:
                first_row = group.head(1).copy()
                first_row.loc[:, new_column] = '-'
                consensus_df = pd.concat([consensus_df, first_row])
                df = df[df['Query#'] != name]
                # print(f"{group_name} - {test_result}. {group_size} rows. Frequency - {frequency_of_modal_value}. Threshold - {threshold}. Mode - {modal_value}.")
        grouped = df_filtered.groupby('Query#')

    return df, consensus_df

def process_consensus_df(pre_consensus_df, consensus_df, shared_matches_df, columns):
    pre_consensus_df, matching_df = combine_duplicates(pre_consensus_df, columns)
    consensus_df = concatenate_dfs(consensus_df, matching_df)
    if check_queries_complete_one(pre_consensus_df):
        check_queries_complete_two(consensus_df, shared_matches_df)
        consensus_df = remove_duplicated_processed_queries(consensus_df)
        consensus_df, new_columns = add_missing_columns(consensus_df, columns)
        consensus_df = replace_nans(consensus_df, new_columns)
        # continue processing
        return consensus_df
    else:
        print('stopping at process_consensus_df')
        # stop program

def combine_duplicates(pre_consensus_df, columns):
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

def check_queries_complete_one(pre_consensus_df):
    if pre_consensus_df.empty:
        return True
    else:
        unprocessed_groups = pre_consensus_df['Query#'].nunique()
        print(f'There are {unprocessed_groups} unprocessed groups still to be processed')
        return False

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

def remove_duplicated_processed_queries(consensus_df):
    unique_df = consensus_df.drop_duplicates(subset='Query#')
    duplicated_df = consensus_df[~consensus_df.index.isin(unique_df.index)]
    if duplicated_df.empty == False:
        rows_deleted = len(duplicated_df)
        duplicated_queries = list(set(duplicated_df['Query#']))
        print(f"{rows_deleted} rows have been deleted from the consensus_df. The following queries were duplicated: {duplicated_queries}.")
    check_unique_df = list(set(consensus_df['Query#'])) == list(set(unique_df['Query#']))
    if check_unique_df:
        return unique_df

def add_missing_columns(consensus_df, columns):
    new_columns = [x.lower() for x in columns]
    consensus_columns = list(consensus_df.columns.values)
    if consensus_columns != new_columns:
        missing_columns = (set(new_columns).difference(consensus_columns))
        for str in missing_columns:
            consensus_df = consensus_df.assign(**{str: '-'})
    return consensus_df, new_columns

def replace_nans(consensus_df, new_columns):
    consensus_df[new_columns] = consensus_df[new_columns].fillna('-')
    return consensus_df

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

def concatenate_dfs(*dfs):
    column_headers = [list(df.columns) for df in dfs]
    if all(headers == column_headers[0] for headers in column_headers):
        concatenated_df = pd.concat(dfs)
        return concatenated_df
    else:
        df_names = ', '.join([str(df) for df in dfs])
        print(f"{df_names} could not be concatenated as the column headers do not match. Please rectify")

main()