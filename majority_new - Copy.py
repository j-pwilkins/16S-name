import pandas as pd
import numpy as np
import os
import shutil
import sys
from datetime import datetime

def main():
    # Create new 'data' folder to house test files
    if os.path.exists("data"):
        shutil.rmtree("data")
    os.mkdir("data")
    shared_matches_df = read_csv('test.csv')
    add_consensus_sequences(shared_matches_df)
    # Move .csv files that match the specified format to the "data" folder
    for file in os.listdir():
        if file.endswith(".csv") and file[0].isupper() and file[1].isdigit():
            shutil.move(file, "data")


def add_consensus_sequences(shared_matches_df):
    majority_df = add_majority_consensus_sequences(shared_matches_df)

def add_majority_consensus_sequences(shared_matches_df):
    majority_df = shared_matches_df.copy()
    columns = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N']
    majority_df, rejected_df = add_majority_column_a(shared_matches_df, columns)
    write_csv(majority_df, 'majority.csv')
    write_csv(rejected_df, 'rejected.csv')
    return majority_df

def calculate_threshold(group_size, frequency_of_modal_value, minimum_threshold=0.5):
    threshold = group_size * minimum_threshold
    test_result = frequency_of_modal_value / group_size > minimum_threshold
    return threshold, test_result

def add_majority_column_a(df, columns):
    # initialize an empty dataframe to store rejected rows
    rejected_df = pd.DataFrame(columns=df.columns)
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
            # print(group_name)
            # if group_name == 'B4':
            #     print(group)
            #     print(df)
            group_name_csv = group_name + '.csv'
            group_name_2_csv = group_name + '_2.csv'
            group_name_3_csv = group_name + '_3.csv'
            write_csv(df, group_name_csv)
            if test_result:
                # df = df.copy()
                df.loc[(df['Query#'] == name) & (df[column] == modal_value), new_column] = modal_value
                df.loc[(df['Query#'] == name) & (df[column] != modal_value), new_column] = 'x'
                write_csv(df, group_name_2_csv)
                df = df[df[new_column] != 'x']
                write_csv(df, group_name_3_csv)
                print(f"{group_name} - {test_result} . {group_size} rows. Frequency - {frequency_of_modal_value}. Threshold - {threshold}. Mode - {modal_value}.")
                df_filtered = df[df[new_column] != 'x']
                grouped = df_filtered.groupby('Query#')  # this line re-group the dataset for the next column iteration
            if not test_result:
                first_row = group.head(1).copy()
                # print(type(first_row))
                first_row.loc[:, new_column] = '-'
                rejected_df = pd.concat([rejected_df, first_row])
                # print(f"group name is {name}")
                df = df[df['Query#'] != name]
                # print(group)
                # print(df)
                print(f"{group_name} - {test_result}. {group_size} rows. Frequency - {frequency_of_modal_value}. Threshold - {threshold}. Mode - {modal_value}.")
        grouped = df_filtered.groupby('Query#')

    return df, rejected_df




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