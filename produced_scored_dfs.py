import pandas as pd
import numpy as np
import os
import shutil
import sys
from datetime import datetime
pd.options.display.max_columns = None

def main():
    euc_input_df = read_inputs()
    euc_input_df = adjust_xls_order(euc_input_df)
    euc_input_df = reformat_scoring(euc_input_df)
    strict_decimal_matrix_df, strict_decimal_ascending_df, r_format_strict_df = produce_matrix(euc_input_df,
                                                                                               selection='Cs',
                                                                                               scoring_column='Decimal Scoring')
    strict_ordinal_matrix_df, strict_ordinal_ascending_df, r_format_strict_df = produce_matrix(euc_input_df,
                                                                                               selection='Cs',
                                                                                               scoring_column='Ordinal Scoring')
    write_csvs(euc_input_df, strict_decimal_matrix_df,
               strict_decimal_ascending_df, strict_ordinal_matrix_df, strict_ordinal_ascending_df)  # allows copy and paste into .xls file

# read inputs required for program
def read_inputs():
    euc_input_df = read_csv('euclidean_scores_input.csv')
    # euc_input_df = pd.read_csv('euclidean_scores_input.csv', na_filter=False)
    return euc_input_df

# This function adjusts the inputted .csv file into the format that it should be at this stage
# This function should not be needed in future when earlier programs are performing as desired
def adjust_xls_order(df):
    # reenters the 'N/A' values (if I use the na false read_csv version I turn all the columns into objects)
    df[['Similarity (%)', 'Vs DB#', 'DB Name', 'DB Accession Number']] = df[
        ['Similarity (%)', 'Vs DB#', 'DB Name', 'DB Accession Number']].fillna('N/A')
    # modifies reference rows so that they don't have values
    df.loc[df['Selected'] == 'DB', ['Decimal Scoring', 'Reverse Scoring', 'Ordinal Scoring']] = '/'
    # orders the df so that the 'Selected' column is in the preferred order
    df['Selected Order'] = df['Selected'].apply(lambda
                                                    x: 1 if x == 'DB' else 2 if x == 'Y' else 3 if x == 'Y.' else 4 if x == 'N' else 5 if x == 'Cs' else 6 if x == 'Cm' else 7 if x == 'Cp' else 8 if x == 'Y-na' else 'NA')
    # orders the df so that the regions are in the preferred order
    df['Region Order'] = np.where(df['Selected'] == 'DB', 1,
                                  np.where(df['Region'] == 'Full', 2,
                                           np.where(df['Region'] == 'V1V3', 3,
                                                    np.where(df['Region'] == 'V3V4', 4,
                                                             np.where(df['Region'] == 'V4', 5, 'NA')))))
    # Copies the tree order column (taken from a db merge in mapping.py) so that the adjustable columns are next to each other
    df['Tree Order'] = df['Tree #']
    # Order the columns - has same effect as performing this ascending order in reverse in Excel
    df = df.iloc[1:, :].sort_values(by=['Tree Order', 'Region Order', 'Selected Order'],
                                    ascending=[True, True, True])
    # add column that can be ordered at later date to return to this order (after any fiddling by user)
    df['JP Order'] = np.arange(1, len(df) + 1)
    return df

def reformat_scoring(df):
    # reverse the ordinal list
    ordinal_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    reverse_ordinal_list = [15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
    df['Ordinal Scoring'] = pd.to_numeric(df['Ordinal Scoring'], errors='coerce')
    df['Ordinal Scoring'].replace(ordinal_list, reverse_ordinal_list, inplace=True)
    df.loc[df['DB/Query #'].notna(), 'Ordinal Scoring'] = df['Ordinal Scoring'].fillna(
        value='/')
    # Change the values in 'Decimal Scoring' column that are 1 to 0.5
    df['Decimal Scoring'] = df['Decimal Scoring'].apply(lambda x: 0.5 if x == 1 else x)
    return df

def produce_matrix(df, selection, scoring_column):
    # keep only selected columns
    reduced_df = df.loc[:, ['DB/Query #', 'Tree #', 'Region', 'Selected', scoring_column]]
    # keep only rows matching certain values
    matrix_selection = ['Y', selection]
    reduced_df = reduced_df[reduced_df['Selected'].isin(matrix_selection)]
    # produce trimmed regional dfs
    regional_columns = ['DB/Query #', 'Tree #', scoring_column]
    full_df = reduced_df[(reduced_df['Region'] == 'Full')]
    full_df = full_df.loc[:, regional_columns]
    full_df.rename(columns={scoring_column: 'Full'}, inplace=True)
    v1v3_df = reduced_df[(reduced_df['Region'] == 'V1V3')]
    v1v3_df = v1v3_df.loc[:, regional_columns]
    v1v3_df.rename(columns={scoring_column: 'V1V3', 'DB/Query #': 'DB/Q'}, inplace=True)
    v3v4_df = reduced_df[(reduced_df['Region'] == 'V3V4')]
    v3v4_df = v3v4_df.loc[:, regional_columns]
    v3v4_df.rename(columns={scoring_column: 'V3V4', 'DB/Query #': 'DB/Q'}, inplace=True)
    v4_df = reduced_df[(reduced_df['Region'] == 'V4')]
    v4_df = v4_df.loc[:, regional_columns]
    v4_df.rename(columns={scoring_column: 'V4', 'DB/Query #': 'DB/Q'}, inplace=True)
    # merge the regional dfs, and drop duplicated columns
    matrix_df = full_df.merge(v1v3_df, on='Tree #')
    matrix_df = matrix_df.merge(v3v4_df, on='Tree #')
    matrix_df = matrix_df.merge(v4_df, on='Tree #')
    matrix_df = matrix_df.loc[:, ['DB/Query #', 'Tree #', 'Full', 'V1V3', 'V3V4', 'V4']]
    # create df of ascending  values for each region
    full_df = full_df.loc[:, 'Full'].sort_values(ascending=True).reset_index(drop=True)
    v1v3_df = v1v3_df.loc[:, 'V1V3'].sort_values(ascending=True).reset_index(drop=True)
    v3v4_df = v3v4_df.loc[:, 'V3V4'].sort_values(ascending=True).reset_index(drop=True)
    v4_df = v4_df.loc[:, 'V4'].sort_values(ascending=True).reset_index(drop=True)
    ascending_df = pd.concat([full_df, v1v3_df, v3v4_df, v4_df], axis=1)
    # name dfs, for writing to .csv files
    matrix_df.name = selection + '_' + scoring_column.strip(' Scoring') + '_Matrix.csv'
    ascending_df.name = selection + '_' + scoring_column.strip(' Scoring') + '_Ascending.csv'

    return matrix_df, ascending_df, reduced_df

# write out .csv files. They can then be combined into one .xls worksheet. combined is copy-pasted (values only) into output sheet
def write_csvs(euc_input_df, strict_decimal_matrix_df,
               strict_decimal_ascending_df, strict_ordinal_matrix_df, strict_ordinal_ascending_df):
    write_csv(euc_input_df, 'euctest.csv')
    write_csv(strict_decimal_matrix_df, strict_decimal_matrix_df.name)
    write_csv(strict_decimal_ascending_df, strict_decimal_ascending_df.name)
    write_csv(strict_ordinal_matrix_df, strict_ordinal_matrix_df.name)
    write_csv(strict_ordinal_ascending_df, strict_ordinal_ascending_df.name)

## helper functions
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