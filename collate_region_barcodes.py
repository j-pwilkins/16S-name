import pandas as pd
import numpy as np
from datetime import datetime
pd.options.display.max_columns = None

def main():
    Full_df, V1V3_df, column_headers = read_inputs()
    tree_df = produce_tree_df(Full_df)
    updated_df = update_regional_df(V1V3_df, tree_df, column_headers)
    Full_df = format_full_df(Full_df, column_headers)
    combined_df = pd.concat([Full_df, updated_df], axis=0, ignore_index=True)
    check_combination(Full_df, updated_df, combined_df)
    combined_df = format_combined_df(combined_df)
    write_csv(combined_df, 'combined.csv')

def produce_tree_df(Full_df):
    tree_df = Full_df[['Accession Number', 'Tree#']]
    return tree_df

def update_regional_df(V1V3_df, tree_df, column_headers):
    merged_df = pd.merge(V1V3_df, tree_df, on='Accession Number')
    merged_df['Seq Length'] = merged_df['16S rRNA sequence'].str.len()
    merged_df = merged_df.rename(
        columns={'Similarity(%)': 'Similarity (%)', 'DB Name': 'Vs Name', 'DB Accession Number': 'Vs ID'})
    merged_df = merged_df[column_headers]
    return merged_df

def format_full_df(Full_df, column_headers):
    Full_df = Full_df.rename(
        columns={'DB#': 'Query#'})
    Full_df = Full_df[column_headers]
    return Full_df

def check_combination(Full_df, updated_df, combined_df):
    Full_length = len(Full_df)
    Updated_length = len(updated_df)
    Combined_length = len(combined_df)
    print(f'The Full df is {Full_length} rows long, Updated df is {Updated_length} rows, and the Combined df is {Combined_length}.')

def format_combined_df(combined_df):
    combined_df['ID?'] = combined_df['Region'].apply(lambda x: '16S' if x == 'Full' else 'N')
    combined_df['ID Order'] = combined_df['ID?'].apply(lambda x: 1 if x == '16S' else 2)
    combined_df['Selected Order'] = combined_df['Selected'].apply(lambda
                                                                      x: 1 if x == '16S' else 2 if x == 'Y' else 3 if x == 'Y.' else 4 if x == 'N' else 5 if x == 'Cs' else 6 if x == 'Cm' else 7 if x == 'Cp' else 8 if x == 'Y-na' else 'NA')
    combined_df['Region Order'] = combined_df['Region'].apply(
        lambda x: 1 if x == 'Full' else 2 if x == 'V1V3' else 3 if x == 'V3V4' else 4 if x == 'V4' else 'NA')
    combined_df['Tree Order'] = combined_df['Tree#']
    combined_df = combined_df.sort_values(by=['Query#', 'Tree Order', 'Region Order', 'Selected Order', 'ID Order'], ascending=[True, True, True, True, True])
    combined_df['JP Order'] = np.arange(1, len(combined_df) + 1)
    return combined_df



## functions for produce_shared_matches_df
# read the inputs produced from mapping - these will be passed to the function when integrating into the the new mapping.py
def read_inputs():
    Full_df = pd.read_csv('full_16S_database.csv', na_filter=False)
    # V1V3_df = read_csv('mapped_barcodes_for_V1V3.csv')
    V1V3_df = pd.read_csv('mapped_barcodes_for_V1V3.csv', na_filter=False)
    column_headers = ['Query#', 'Hun#', 'Fas#', 'Tree#', 'Accession Number', 'Name', 'Region', 'Shared Matches', 'Selected', 'Similarity (%)', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'Vs DB#', 'Vs Name', 'Vs ID', 'Seq Length', '16S rRNA sequence']
    output_columns = ['Query#', 'Hun#', 'Fas#', 'Accession Number', 'Name', 'Shared Matches', 'Selected', 'Similarity(%)',
                    'Vs DB#',
                    'DB Name', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
                    'DB Accession Number', '16S rRNA sequence']
    return Full_df, V1V3_df, column_headers

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