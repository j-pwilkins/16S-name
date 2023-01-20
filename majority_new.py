import pandas as pd
import numpy as np
import os
import shutil
import sys
from datetime import datetime

def main():
    shared_matches_df = read_csv('test.csv')
    add_consensus_sequences(shared_matches_df)


def add_consensus_sequences(shared_matches_df):
    pre_consensus_df = add_majority_consensus_sequences(shared_matches_df)

def add_majority_consensus_sequences(shared_matches_df):
    pre_consensus_df = shared_matches_df.copy()
    columns = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N']
    pre_consensus_df, consensus_df = add_majority_column_a(shared_matches_df, columns)
    write_csv(pre_consensus_df, 'majority.csv')
    write_csv(consensus_df, 'rejected.csv')
    consensus_df = process_consensus_df(pre_consensus_df, consensus_df, shared_matches_df, columns)
    return pre_consensus_df

def calculate_threshold(group_size, frequency_of_modal_value, minimum_threshold=0.5):
    threshold = group_size * minimum_threshold
    test_result = frequency_of_modal_value / group_size > minimum_threshold
    return threshold, test_result

def add_majority_column_a(df, columns):
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
                print(f"{group_name} - {test_result}. {group_size} rows. Frequency - {frequency_of_modal_value}. Threshold - {threshold}. Mode - {modal_value}.")
        grouped = df_filtered.groupby('Query#')

    return df, consensus_df

def process_consensus_df(pre_consensus_df, consensus_df, shared_matches_df, columns):
    if check_queries_complete_one(pre_consensus_df):
        check_queries_complete_two(consensus_df, shared_matches_df)
        remove_duplicated_processed_queries(consensus_df)
        # continue processing
    else:
        print('stopping at process_consensus_df')
        # stop program

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