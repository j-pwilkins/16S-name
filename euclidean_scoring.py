import pandas as pd
import numpy as np
from datetime import datetime
pd.options.display.max_columns = None

def main():
    complete_df = read_inputs()
    eu_score_df = create_scoring_column(complete_df)
    strict_decimal_matrix_df, strict_decimal_ascending_df, r_format_strict_df = produce_matrix(eu_score_df,
                                                                                               selection='Cs',
                                                                                               scoring_column='Decimal Scoring')
    strict_ordinal_matrix_df, strict_ordinal_ascending_df, r_format_strict_df = produce_matrix(eu_score_df,
                                                                                               selection='Cs',
                                                                                               scoring_column='Ordinal Scoring')
    write_csvs(eu_score_df, strict_decimal_matrix_df,
               strict_decimal_ascending_df, strict_ordinal_matrix_df, strict_ordinal_ascending_df)

##### L1
def create_scoring_column(complete_df):
    scoring_df = add_ordinal_column(complete_df)
    write_csv(scoring_df, 'test1.csv')
    scoring_df = score_sequences_ordinal(scoring_df)
    write_csv(scoring_df, 'test2.csv')
    scoring_df = add_decimal_scoring_column(scoring_df)
    scoring_df = replace_nan_values(scoring_df)
    return scoring_df

### L2
def add_ordinal_column(df):
    # add ordinal column directly after N
    df.insert(df.columns.get_loc('N')+1, 'Ordinal Scoring', [np.nan]+['x']*(len(df)-1))
    # Add '/' value for reference row
    df["Ordinal Scoring"] = np.where(df["Selected"] == "DB", '/', df["Ordinal Scoring"])
    return df

def score_sequences_ordinal(scoring_df):
    groups = scoring_df.groupby('DB/Query #')   # Group dataframe by "DB/Query #"
    scoring_cols = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N']
    for name, group in groups:      # Iterate through each group
        # Get the reference row for each group (the row with "Selected" value as "DB")
        ref_row = group.loc[group['Selected'] == 'DB']
        ref_values = ref_row[scoring_cols].values[0]
        # Iterate through each row in the group
        for i, row in group.iterrows():
            # Skip the reference row
            if row['Selected'] != 'DB':
                # Keep track of the number of columns that agree with the reference values
                num_agreed = 0
                # Iterate through each scoring column
                for j, col in enumerate(scoring_cols):
                    # Check if the value in the current column is equal to the reference value
                    if row[col] == ref_values[j]:
                        num_agreed += 1
                    else:
                        # Assign the ordinal score
                        scoring_df.at[i, 'Ordinal Scoring'] = str(15 - num_agreed)
                        # Exit the loop if a difference is found
                        break
                # If no differences are found, assign the minimum ordinal score
                else:
                    scoring_df.at[i, 'Ordinal Scoring'] = '1'
    # Convert the 'Ordinal Scoring' column to numeric
    scoring_df['Ordinal Scoring'] = pd.to_numeric(scoring_df['Ordinal Scoring'], errors='coerce')
    # Replace nan values in 'Ordinal Scoring' column with '/' for rows with 'DB' value in 'Selected' column
    scoring_df.loc[(scoring_df['Selected'] == 'DB') & (scoring_df['Ordinal Scoring'].isnull()), 'Ordinal Scoring'] = '/'
    # Return the updated dataframe
    return scoring_df

### L2
def add_decimal_scoring_column(scoring_df):
    # create new column as copy of ordinal
    scoring_df.insert(loc=scoring_df.columns.get_loc("Ordinal Scoring"), column="Decimal Scoring",
                      value=scoring_df['Ordinal Scoring'])
    reverse_ordinal_list = [15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
    decimal_list = [0.6, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.02, 0.01, 0.005, 0.004, 0.003, 0.002, 0.001, 0]
    # switch the above lists for new decimal scoring
    scoring_df['Decimal Scoring'] = pd.to_numeric(scoring_df['Decimal Scoring'], errors='coerce')
    scoring_df['Decimal Scoring'].replace(reverse_ordinal_list, decimal_list, inplace=True)
    scoring_df.loc[scoring_df['DB/Query #'].notna(), 'Decimal Scoring'] = scoring_df['Decimal Scoring'].fillna(
        value='/')
    return scoring_df

### L2
# This function removes any nan values that have appeared due to the failure of the 'N/A' values being read in properly
# If I do read them in, the score_sequences_ordinal is currently breaking, so using this as a workaround
def replace_nan_values(scoring_df):
    scoring_df.loc[scoring_df['DB/Query #'].notna(), 'Similarity (%)'] = scoring_df['Similarity (%)'].fillna(
        value='N/A')
    scoring_df.loc[scoring_df['DB/Query #'].notna(), 'Vs DB#'] = scoring_df['Vs DB#'].fillna(
        value='N/A')
    scoring_df.loc[scoring_df['DB/Query #'].notna(), 'DB Name'] = scoring_df['DB Name'].fillna(
        value='N/A')
    scoring_df.loc[scoring_df['DB/Query #'].notna(), 'DB Accession Number'] = scoring_df['DB Accession Number'].fillna(
        value='N/A')
    return scoring_df

##### L1
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

##### L1
def read_inputs():
    complete_df = read_csv('collected_barcodes.csv')
    return complete_df

def write_csvs(eu_score_df, strict_decimal_matrix_df,
               strict_decimal_ascending_df, strict_ordinal_matrix_df, strict_ordinal_ascending_df):
    write_csv(eu_score_df,'euclidean_scoring.csv')
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