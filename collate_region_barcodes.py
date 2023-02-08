import pandas as pd
import numpy as np
import os
import shutil
import sys
from datetime import datetime
pd.options.display.max_columns = None

def main():
    db_df, V1V3_df, V3V4_df, V4_df, full_df, column_headers = read_inputs()
    # tree_df = produce_tree_df(Full_df)
    V1V3_df = update_regional_df(V1V3_df, column_headers)
    V3V4_df = update_regional_df(V3V4_df, column_headers)
    V4_df = update_regional_df(V4_df, column_headers)
    full_df = update_regional_df(full_df, column_headers)
    db_df = format_db_df(db_df, column_headers)
    combined_df = pd.concat([db_df, V1V3_df, V3V4_df, V4_df, full_df], axis=0, ignore_index=True)
    check_combination(db_df, V1V3_df, V3V4_df, V4_df, full_df, combined_df)
    combined_df = format_combined_df(combined_df)
    write_csvs(V1V3_df, V3V4_df, V4_df, full_df, db_df, combined_df)

def produce_tree_df(Full_df):
    tree_df = Full_df[['Accession Number', 'Tree#']]
    return tree_df

def update_regional_df(updated_df, column_headers):
    # updated_df = pd.merge(V1V3_df, tree_df, on='Accession Number')
    updated_df['Seq Length'] = updated_df['16S rRNA sequence'].str.len()
    updated_df = updated_df.rename(
        columns={'Q Accession Number': 'Accession Number', 'Q Name': 'Name', 'Similarity(%)': 'Similarity (%)', 'DB Name': 'Vs Name', 'DB Accession Number': 'Vs ID'})
    write_csv(updated_df, 'updatetest.csv')
    updated_df = updated_df[column_headers]
    return updated_df

def format_db_df(db_df, column_headers):
    db_df = db_df.rename(
        columns={'DB#': 'Query#'})
    db_df = db_df[column_headers]
    return db_df

def check_combination(db_df, V1V3_df, V3V4_df, V4_df, full_df, combined_df):
    db_length = len(db_df)
    V1V3_length = len(V1V3_df)
    V3V4_length = len(V3V4_df)
    V4_length = len(V4_df)
    full_length = len(full_df)
    Combined_length = len(combined_df)
    print(f'The Database df is {db_length} rows long, V1V3 df is {V1V3_length} rows, V3V4 is {V3V4_length}, V4 is {V4_length}, Full is {full_length} and the Combined df is {Combined_length}.')


def format_combined_df(combined_df):
    # combined_df['ID?'] = combined_df['Region'].apply(lambda x: '16S' if x == 'Full' else 'N')
    # combined_df['ID Order'] = combined_df['ID?'].apply(lambda x: 1 if x == '16S' else 2)
    combined_df['Selected Order'] = combined_df['Selected'].apply(lambda x: 1 if x == 'DB' else 2 if x == 'Y' else 3 if x == 'Y.' else 4 if x == 'N' else 5 if x == 'Cs' else 6 if x == 'Cm' else 7 if x == 'Cp' else 8 if x == 'Y-na' else 'NA')
    combined_df['Region Order'] = np.where(combined_df['Selected'] == 'DB', 1,
                                           np.where(combined_df['Region'] == 'Full', 2,
                                                    np.where(combined_df['Region'] == 'V1V3', 3,
                                                             np.where(combined_df['Region'] == 'V3V4', 4,
                                                                      np.where(combined_df['Region'] == 'V4', 5, 'NA')))))
    combined_df['Tree Order'] = combined_df['Tree#']
    # combined_df = combined_df.sort_values(by=['Query#', 'Tree Order', 'Region Order', 'Selected Order', 'ID Order'], ascending=[True, True, True, True, True])
    combined_df = combined_df.sort_values(by=['Tree Order', 'Region Order', 'Selected Order'],
                                          ascending=[True, True, True])
    combined_df['JP Order'] = np.arange(1, len(combined_df) + 1)
    return combined_df



def write_csvs(V1V3_df, V3V4_df, V4_df, full_df, db_df, combined_df):
    write_csv(V1V3_df, 'V1V3.csv')
    write_csv(V3V4_df, 'V3V4.csv')
    write_csv(V4_df, 'V4.csv')
    write_csv(full_df, 'Full.csv')
    write_csv(db_df, 'DB.csv')
    write_csv(combined_df, 'combined.csv')

## functions for produce_shared_matches_df
# read the inputs produced from mapping - these will be passed to the function when integrating into the the new mapping.py
def read_inputs():
    db_df = pd.read_csv('full_16S_database.csv', na_filter=False)
    # V1V3_df = read_csv('mapped_barcodes_for_V1V3.csv')
    V1V3_df = pd.read_csv('Complete_list_of_sequences_V1V3_vs_db_497.csv', na_filter=False)
    V3V4_df = pd.read_csv('Complete_list_of_sequences_V3V4_vs_db_497.csv', na_filter=False)
    V4_df = pd.read_csv('Complete_list_of_sequences_V4_vs_db_497.csv', na_filter=False)
    full_df = pd.read_csv('Complete_list_of_sequences_Full_vs_db_497.csv', na_filter=False)
    column_headers = ['Query#', 'Hun#', 'Fas#', 'Tree#', 'Accession Number', 'Name', 'Region', 'Shared Matches', 'Selected', 'Similarity (%)', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'Vs DB#', 'Vs Name', 'Vs ID', 'Seq Length', '16S rRNA sequence']
    output_columns = ['Query#', 'Hun#', 'Fas#', 'Accession Number', 'Name', 'Shared Matches', 'Selected', 'Similarity(%)',
                    'Vs DB#',
                    'DB Name', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
                    'DB Accession Number', '16S rRNA sequence']
    return db_df, V1V3_df, V3V4_df, V4_df, full_df, column_headers

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



main()
