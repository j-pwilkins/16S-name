import pandas as pd
import numpy as np
from datetime import datetime


def main():
    start_time = datetime.now()
    region = 'test'
    database = 'db_497.csv'
    database_df = read_csv(database)
    summary_df = pd.read_csv('summary_test.csv', na_filter=False)
    shared_df = pd.read_csv('shared_test.csv', na_filter=False)
    barcoded_df, missing_df = add_barcode(summary_df, database_df)
    # write_csv(barcoded_df, 'barcoded.csv')
    merged_df = merge_shared_matches_to_summary(barcoded_df, shared_df, region, missing_df)
    write_csv(merged_df, 'merged.csv')
    summary_statement(start_time)

def add_barcode(summary_df, database_df):
    summary_df, missing_df = remove_na_rows(summary_df)
    write_csv(missing_df, 'mis.csv')
    summary_df = reorganise_summary_df(summary_df)
    database_df = reorganise_database_df(database_df)
    barcoded_df = merge_summary_with_db(summary_df, database_df)
    # write_csv(barcoded_df, 'b1.csv')
    barcoded_df = check_merge_correct(barcoded_df)
    # write_csv(barcoded_df, 'b2.csv')
    barcoded_df = reorganise_barcoded_df(barcoded_df)
    write_csv(barcoded_df, 'b3.csv')
    return barcoded_df, missing_df

def merge_shared_matches_to_summary(barcoded_df, shared_df, region, missing_df):
    edit_shared_df(shared_df)
    missing_df = format_missing_df(missing_df)
    write_csv(missing_df, 'mis2.csv')
    if check_queries_match(barcoded_df, shared_df):
        edit_shared_df(shared_df)
        barcoded_df = edit_barcoded_df(barcoded_df)
        write_csv(barcoded_df, 'barcoded2.csv')
        merged_df = merge_dfs(barcoded_df, shared_df, missing_df)
        add_identifying_column(merged_df, region)
    else:
        print("Error: Query# values in barcoded_df and shared_df do not match.")
        exit()

    return merged_df

def remove_na_rows(summary_df):
    # create an empty dataframe to store the rows from summary_df that contain 'N/A' in the 'Vs DB#' column
    missing_df = pd.DataFrame(columns=summary_df.columns)

    # counter to track the number of rows removed
    count = 0

    # iterate through the rows of summary_df
    for index, row in summary_df.iterrows():
        # if the row contains the string 'N/A' in the 'Vs DB#' column
        if row['Vs DB#'] == 'N/A':
            # create a new dataframe from row
            df = pd.DataFrame(row).T
            # append df to missing_df
            missing_df = pd.concat([missing_df, df], ignore_index=True)
            # drop the row from summary_df
            summary_df.drop(index, inplace=True)
            # increment the counter
            count += 1

    # print the statement
    print(f"{count} rows that could not be compared to the DB have been removed, they will be re-integrated later")

    # return the modified summary_df and missing_df
    return summary_df, missing_df


def reorganise_summary_df(summary_df):
    # Rename columns
    summary_df = summary_df.rename(columns={'DB Name': 'Name', 'DB Accession Number': 'Accession Number', 'Vs Name': 'DB Name', 'Alternative Matches': 'Shared Matches' })

    # Convert columns to the types they would have been without the 'N/A' rows
    summary_df['Vs DB#'] = summary_df['Vs DB#'].astype('int64')
    summary_df['Shared Matches'] = summary_df['Shared Matches'].astype('int64')
    summary_df['Similarity(%)'] = summary_df['Similarity(%)'].astype('float64')

    # Add 1 to every value in the 'Shared Matches' column
    summary_df['Shared Matches'] = summary_df['Shared Matches']+1

    # Add a 'Selected' column and give every row the string 'Y'
    summary_df = summary_df.assign(Selected=lambda x: 'Y')

    # Reorder columns
    summary_df = summary_df[['Query#', 'Hun#', 'Fas#', 'Accession Number', 'Name', 'Shared Matches', 'Selected', 'Similarity(%)', 'Vs DB#', 'DB Name', '16S rRNA sequence']]

    return summary_df


def reorganise_database_df(database_df):
    # Keep only specified columns in specified order
    database_df = database_df[['DB#', 'DB Name', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'DB Accession Number']]

    # Rename 'DB#' column as 'Vs DB#'
    database_df = database_df.rename(columns={'DB#': 'Vs DB#', 'DB Name': 'dbDB Name'})

    return database_df

def merge_summary_with_db(summary_df, database_df):
    # Rename the 'DB#' column in the database_df to 'Vs DB#' to match the column name in the summary_df
    database_df = database_df.rename(columns={'DB#': 'Vs DB#'})

    # Perform the merge using an inner join to keep only common rows
    barcoded_df = summary_df.merge(database_df, on='Vs DB#', how='inner')

    return barcoded_df

def check_merge_correct(barcoded_df):
    # Iterate through each row and check if 'DB Name' and 'dbDB Name' match
    match = True
    for index, row in barcoded_df.iterrows():
        if row['DB Name'] != row['dbDB Name']:
            match = False
            break

    # Print the result and delete 'dbDB Name' if all rows match
    if match:
        print("The merged dfs match")
        del barcoded_df['dbDB Name']
    else:
        print("The merged dfs DO NOT match")

    return barcoded_df

def reorganise_barcoded_df(barcoded_df):
    # Get the index of the '16S rRNA sequence' column
    index = barcoded_df.columns.get_loc('16S rRNA sequence')

    # Create a new list of column names with '16S rRNA sequence' at the end
    new_columns = list(barcoded_df.columns[:index]) + list(barcoded_df.columns[index+1:]) + ['16S rRNA sequence']

    # Reorder the columns
    barcoded_df = barcoded_df[new_columns]

    return barcoded_df

def remove_db_string(s):
    if s == 'N/A':
        return s
    return int(s.replace('DB#_', ''))

def format_missing_df(missing_df):
    missing_df = missing_df.rename(
        columns={'DB Name': 'Name', 'DB Accession Number': 'Accession Number', 'Vs Name': 'DB Name',
                 'Alternative Matches': 'Shared Matches'}
    )
    missing_df['DB Accession Number'] = 'N/A'
    missing_df['Selected'] = 'Y-na'
    missing_df = missing_df.reindex(
        columns=['Query#', 'Hun#', 'Fas#', 'Accession Number', 'Name', 'Shared Matches', 'Selected', 'Similarity(%)',
                 'Vs DB#', 'DB Name', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
                 'DB Accession Number', '16S rRNA sequence'],
        fill_value='-'
    )
    return missing_df


def edit_shared_df(shared_df):
    shared_df.loc[shared_df['Selected'] == 'Y', 'Selected'] = 'Y.'
    shared_df['Vs DB#'] = shared_df['Vs DB#'].astype(str)
    shared_df['Vs DB#'] = shared_df['Vs DB#'].apply(remove_db_string)
    return shared_df

def edit_barcoded_df(barcoded_df):
    barcoded_df = barcoded_df[barcoded_df['Shared Matches'] <= 1]
    return barcoded_df

def check_queries_match(barcoded_df, shared_df):
    # Get a list of Query# values from barcoded_df where the value in the Shared Matches column is >1
    queries = barcoded_df[barcoded_df['Shared Matches'] > 1]['Query#'].tolist()

    # Check if all values in 'queries' are also in the 'Query#' column of 'shared_df'
    if not all(query in shared_df['Query#'].tolist() for query in queries):
        print('shared_df is missing a/some queries')
        return False

    # Check if there are any values in the 'Query#' column of 'shared_df' that are not in 'queries'
    not_in_queries = [query for query in shared_df['Query#'].tolist() if query not in queries]
    if not_in_queries:
        print('shared_df has more shared queries than barcoded_df')
        return False

    print('The barcoded_df and shared_dfs agree')
    return True

def merge_dfs(barcoded_df, shared_df, missing_df):
    merged_df = pd.concat([barcoded_df, shared_df, missing_df])
    # Re-order the merged dataframe by the 'Query#' column in ascending order
    merged_df = merged_df.sort_values(by='Query#', ascending=True)
    return merged_df

def add_identifying_column(merged_df, region):
    merged_df.insert (5, "Region", region)
    return merged_df

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