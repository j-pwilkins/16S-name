import pandas as pd
import numpy as np
from datetime import datetime
pd.options.display.max_columns = None

def main():
    complete_df = read_inputs()
    eu_score_df = create_scoring_column(complete_df)
    write_csvs(eu_score_df)

def create_scoring_column(complete_df):
    scoring_df = add_column(complete_df)
    scoring_df = add_16S_values(scoring_df)
    scoring_df = score_sequences_ordinal(scoring_df)
    scoring_df = add_scoring_columns(scoring_df)
    scoring_df = rescore_decimal_scoring_column(scoring_df)
    scoring_df = rescore_reverse_scoring_column(scoring_df)
    scoring_df = remove_nan_values(scoring_df)
    return scoring_df

def add_column(complete_df):
    complete_df.insert(complete_df.columns.get_loc('N')+1, 'Ordinal Scoring', [np.nan]+['x']*(len(complete_df)-1))
    return complete_df

def add_16S_values(scoring_df):
    scoring_df["Ordinal Scoring"] = np.where(scoring_df["Selected"] == "16S", '/', scoring_df["Ordinal Scoring"])
    return scoring_df

def score_sequences_ordinal(scoring_df):
    groups = scoring_df.groupby('DB/Query #')
    scoring_cols = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N']
    for name, group in groups:
        ref_row = group.loc[group['Selected'] == 'DB']
        ref_values = ref_row[scoring_cols].values[0]
        for i, row in group.iterrows():
            if row['Selected'] != 'DB':
                for j, col in enumerate(scoring_cols):
                    if row[col] != ref_values[j]:
                        scoring_df.at[i, 'Ordinal Scoring'] = str(j+1)
                        break
                else:
                    scoring_df.at[i, 'Ordinal Scoring'] = '15'
    return scoring_df

def add_scoring_columns(scoring_df):
    scoring_df.insert(loc=scoring_df.columns.get_loc("Ordinal Scoring"), column="Decimal Scoring", value=scoring_df['Ordinal Scoring'])
    scoring_df.insert(loc=scoring_df.columns.get_loc("Ordinal Scoring"), column="Reverse Scoring", value=scoring_df['Ordinal Scoring'])
    return scoring_df

def rescore_decimal_scoring_column(scoring_df):
    ordinal_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    decimal_list = [1, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.02, 0.01, 0.005, 0.004, 0.003, 0.002, 0.001, 0]
    scoring_df['Decimal Scoring'] = pd.to_numeric(scoring_df['Decimal Scoring'], errors='coerce')
    scoring_df['Decimal Scoring'].replace(ordinal_list, decimal_list, inplace=True)
    scoring_df.loc[scoring_df['DB/Query #'].notna(), 'Decimal Scoring'] = scoring_df['Decimal Scoring'].fillna(
        value='/')
    return scoring_df

def rescore_reverse_scoring_column(scoring_df):
    ordinal_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    reverse_list = [0, 600, 700, 800, 850, 900, 950, 980, 990, 995, 996, 997, 998, 999, 1000]
    scoring_df['Reverse Scoring'] = pd.to_numeric(scoring_df['Reverse Scoring'], errors='coerce')
    scoring_df['Reverse Scoring'].replace(ordinal_list, reverse_list, inplace=True)
    scoring_df.loc[scoring_df['DB/Query #'].notna(), 'Reverse Scoring'] = scoring_df['Reverse Scoring'].fillna(
        value='/')
    return scoring_df


# This function removes any nan values that have appeared due to the failure of the 'N/A' values being read in properly
# If I do read them in, the score_sequences_ordinal is currently breaking, so using this as a workaround
def remove_nan_values(scoring_df):
    scoring_df.loc[scoring_df['DB/Query #'].notna(), 'Similarity (%)'] = scoring_df['Similarity (%)'].fillna(
        value='N/A')
    scoring_df.loc[scoring_df['DB/Query #'].notna(), 'Vs DB#'] = scoring_df['Vs DB#'].fillna(
        value='N/A')
    scoring_df.loc[scoring_df['DB/Query #'].notna(), 'DB Name'] = scoring_df['DB Name'].fillna(
        value='N/A')
    scoring_df.loc[scoring_df['DB/Query #'].notna(), 'DB Accession Number'] = scoring_df['DB Accession Number'].fillna(
        value='N/A')
    return scoring_df

def read_inputs():
    complete_df = read_csv('collected_barcodes.csv')
    # complete_df = pd.read_csv('collected_barcodes.csv', na_filter=False)
    return complete_df

def write_csvs(eu_score_df):
    write_csv(eu_score_df,'euclidean_scoring.csv')

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
