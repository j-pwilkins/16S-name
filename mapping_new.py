import pandas as pd
import numpy as np
import os
import shutil
import sys
from datetime import datetime
pd.options.display.max_columns = None

def main():
    start_time = datetime.now()
    input_sequences_list, input_database, vsearch_output, curated_vsearch_output, output_directory, output_filename, summary_filename, output_columns, region, barcode_columns = read_inputs()
    organise_directories(output_directory, input_sequences_list, input_database, vsearch_output)
    database_df, query_df, length_of_query, query_fastaname, database_fastaname = prepare_then_run_vsearch(input_sequences_list, input_database, vsearch_output, region)
    curated_df, summary_df = map_vsearch_output_to_database(curated_vsearch_output, database_df, query_df,
                                                            vsearch_output)
    shared_matches_df, unique_df, sm_unique_queries = produce_shared_matches_df(curated_df, output_columns, region)
    shared_matches_with_consensus_df = add_consensus_sequences(shared_matches_df, sm_unique_queries, barcode_columns,
                                                               output_columns)
    final_df = combine_shared_matches_and_unique_dfs(shared_matches_with_consensus_df, unique_df, database_df)
    length_of_output = file_outputs(final_df, summary_df, output_filename, summary_filename, query_df, curated_vsearch_output, query_fastaname, database_fastaname, vsearch_output, input_sequences_list, input_database)
    summary_statement(length_of_output, start_time)

##### L1
def summary_statement(length_of_output, start_time):
    program_length = datetime.now() - start_time
    print(f'Mapping has completed in {program_length}. {length_of_output} queries have been processed.')

##### L1
def read_inputs():
    input_sequences_list = sys.argv[1]
    input_database = sys.argv[2]
    vsearch_output = 'vsearch_blast6_output.csv'
    curated_vsearch_output = 'Query_vs_All.csv'
    text_insert = '_vs_'
    region = input_sequences_list.rstrip('.csv')
    output_directory = input_sequences_list.rstrip('.csv') + text_insert + input_database.rstrip('.csv')
    output_filename = 'Complete_list_of_sequences_' + output_directory + '.csv'
    summary_filename = 'Short_summary_' + output_directory + '.csv'
    output_columns = [
        'Query#', 'Hun#', 'Fas#', 'Q Accession Number', 'Q Name', 'Region', 'Shared Matches', 'Selected',
        'Similarity(%)', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'Vs DB#', 'DB Name',
        'DB Accession Number', 'Seq Length', '16S rRNA sequence']
    barcode_columns = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N']
    return input_sequences_list, input_database, vsearch_output, curated_vsearch_output, output_directory, output_filename, summary_filename, output_columns, region, barcode_columns

##### L1
def organise_directories(output_directory, input_sequences_list, input_database, vsearch_output):
    if os.path.exists(output_directory):
        shutil.rmtree(output_directory)
        print(f'The previous {output_directory} folder has been deleted')
    os.makedirs(output_directory, exist_ok=True)
    shutil.copy(input_sequences_list, output_directory)
    shutil.copy(input_database, output_directory)
    # shutil.copy(vsearch_output, output_directory)  # temp file until in Kelvin
    os.chdir(output_directory)

##### L1 # Read in .csv files, create fasta, run vsearch
def prepare_then_run_vsearch(input_sequences_list, input_database, vsearch_output, region):
    query_df, database_df, length_of_query, length_of_database, query_fastaname, database_fastaname = import_and_prepare_files(
        input_sequences_list, input_database, region)
    for i in range(length_of_query):
        convert_query_to_fasta(query_df, query_fastaname, i)
    for i in range(length_of_database):
        convert_db_to_fasta(database_df, database_fastaname, i)
    run_vsearch(query_fastaname, database_fastaname, vsearch_output)
    return database_df, query_df, length_of_query, query_fastaname, database_fastaname

### L2
def import_and_prepare_files(input_sequences_list, input_database, region):
    query_df = read_csv(input_sequences_list)
    query_df = query_df.replace(np.nan, '', regex=True)
    query_df.name = input_sequences_list[:-4]  # remove .csv from the string
    database_df = read_csv(input_database)
    length_of_query = len(query_df.index.tolist())
    length_of_database = len(database_df.index.tolist())
    query_fastaname, database_fastaname = create_fasta_filenames(input_sequences_list, input_database)
    query_df.insert(0, 'Query#', np.arange(1, len(query_df) + 1))
    print(
        f'The {input_sequences_list} of {length_of_query} {region} region 16S rRNA sequences has been read in. It will be compared against the database {input_database} of {length_of_database} full length 16S rRNA gene sequences. ')
    return query_df, database_df, length_of_query, length_of_database, query_fastaname, database_fastaname

## L3
def create_fasta_filenames(input_sequences_list, input_database):
    query_base_name = os.path.splitext(input_sequences_list)[0]
    query_fastaname = query_base_name + '.fa'
    database_base_name = os.path.splitext(input_database)[0]
    database_fastaname = database_base_name + '.fa'
    return query_fastaname, database_fastaname

### L2
def convert_query_to_fasta(query_df, query_fastaname, i):
    firstline = '>Query_' + str(query_df['Query#'].values[i]) + ' ' + str(
        query_df['DB Accession Number'].values[i]) + ' ' + str(query_df['DB Name'].values[i])
    secondline = str(query_df['16S rRNA sequence'].values[i])
    with open(query_fastaname, 'a') as f:
        f.write(f"{firstline}\n{secondline}\n")

### L2
def convert_db_to_fasta(database_df, database_fastaname, i):
    firstline = '>DB#_' + str(database_df['DB#'].values[i]) + ' ' + str(
        database_df['DB Accession Number'].values[i]) + " " + str(database_df['DB Name'].values[i])
    secondline = str(database_df['16S rRNA sequence'].values[i])
    with open(database_fastaname, 'a') as f:
        f.write(f"{firstline}\n{secondline}\n")

### L2 # Run vsearch - might need further function for --minwordmatches decision
def run_vsearch(query_fastaname, database_fastaname, vsearch_output):
    cmd = 'vsearch --usearch_global ' + query_fastaname + ' --db ' + database_fastaname + ' --blast6out ' + vsearch_output + ' --acceptall --id 0.0 --maxaccepts 0 --minwordmatches 3'
    os.system(cmd)
    print('vsearch has finished')

##### L1 # uses the vsearch output to map each regional query to full 16S sequence(s)
def map_vsearch_output_to_database(curated_vsearch_output, database_df, query_df, vsearch_output):
    query_vs_target_db_df = parse_vsearch_usearchglobal_blast6(vsearch_output, curated_vsearch_output)
    curated_df, summary_df = map_regions_to_full_sequences(database_df, query_df, query_vs_target_db_df)
    return curated_df, summary_df

### L2 # Take vsearch output file and make it useable
def parse_vsearch_usearchglobal_blast6(vsearch_output, curated_vsearch_output):
    query_vs_target_db_df = pd.read_csv(vsearch_output, usecols=[0, 1, 2], sep='\t',
                                        names=["Query", "DB#", "Similarity(%)"])
    df_length = len(query_vs_target_db_df)
    vsearch_output_order = pd.Series(np.arange(1, df_length + 1, 1))
    query_vs_target_db_df.insert(0, 'Vsearch Order', vsearch_output_order, True)
    query_vs_target_db_df['Q'] = query_vs_target_db_df['Query'].str[6:].astype(int)
    query_vs_target_db_df['T'] = query_vs_target_db_df['DB#'].str[4:].astype(int)
    query_vs_target_db_df = query_vs_target_db_df.sort_values(['Q', 'T'])
    query_vs_target_db_df.to_csv(curated_vsearch_output, sep=',', index=False)
    return query_vs_target_db_df

### L2
def map_regions_to_full_sequences(database_df, query_df, query_vs_target_db_df):
    curated_df = purge_self_references(query_vs_target_db_df)
    curated_df = retain_only_highest_similarities(curated_df)
    curated_df = add_alternative_matches_column(curated_df)
    curated_df = add_database_detail_for_matches(curated_df, database_df)
    curated_df = add_query_detail(curated_df, query_df)
    summary_df = create_short_summary_output(curated_df)
    check_short_summary_complete(query_df, summary_df)
    return curated_df, summary_df

## L3
def purge_self_references(query_vs_target_db_df):
    mask = query_vs_target_db_df['Q'] != query_vs_target_db_df['T']
    return query_vs_target_db_df[mask]

## L3
def retain_only_highest_similarities(curated_df):
    grouped = curated_df.groupby('Query')
    highest_similarity_df = grouped.apply(
        lambda x: x.nlargest(len(x[x['Similarity(%)'] == x['Similarity(%)'].max()]), 'Similarity(%)'))
    return highest_similarity_df.reset_index(drop=True)

## L3
def add_alternative_matches_column(curated_df):
    curated_df["Shared Matches"] = curated_df["Query"].map(curated_df["Query"].value_counts())
    return curated_df

## L3
def add_database_detail_for_matches(curated_df, database_df):
    db_columns_to_be_added_to_df = ['DB#', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
                                    'DB Name', 'DB Accession Number', 'DB Found']
    database_df = database_df[db_columns_to_be_added_to_df]
    database_df = database_df.rename(columns={'DB#': 'T'})
    curated_df = curated_df[['Query', 'DB#', 'Similarity(%)', 'Shared Matches', 'Q', 'T']]
    curated_df = pd.merge(curated_df, database_df, on='T')
    curated_df = curated_df.rename(columns={'T': 'Vs DB#', 'DB Name': 'Vs Name'})
    return curated_df

## L3
def add_query_detail(curated_df, query_df):
    query_df = query_df.rename(
        columns={'Query#': 'Q', 'DB Accession Number': 'Q Accession Number', 'DB Name': 'Q Name'})
    query_columns_to_be_added_to_df = ['Q', 'Hun#', 'Fas#', 'Q Name', 'Q Accession Number', '16S rRNA sequence']
    query_df = query_df.loc[:, query_columns_to_be_added_to_df]
    curated_df = pd.merge(curated_df, query_df, on='Q')
    curated_df['Seq Length'] = curated_df['16S rRNA sequence'].str.len()
    return curated_df

## L3
def create_short_summary_output(curated_df):
    summary_columns = ['Query#', 'Hun#', 'Fas#', 'Q Name', 'Q Accession Number', 'Similarity(%)', 'Vs Name',
                       'DB Accession Number', 'Alternative Matches']
    curated_df = curated_df.rename(columns={'Q': 'Query#'})
    curated_df['Alternative Matches'] = curated_df['Shared Matches'] - 1
    summary_df = curated_df.loc[:, summary_columns]
    summary_df = summary_df.sort_values(by='Query#', ascending=True)
    check_groups_match_before_drop(summary_df, 'Query#', 'Similarity(%)')
    summary_df = summary_df.drop_duplicates(subset='Query#', keep='first')
    return summary_df

## L3
def check_short_summary_complete(query_df, summary_df):
    print(f'A summary of {len(summary_df)} mapped sequences has been created.')
    summary_df.name = 'the summary'
    check_unique_values(query_df, summary_df, 'Query#')

###### L1 # Take the curated_df and split for further analysis
def produce_shared_matches_df(curated_df, output_columns, region):
    output_df = produce_output_df(curated_df, output_columns, region)
    shared_matches_df, unique_df, sm_unique_queries, length_of_shared_matches_df = split_output_df(output_df)
    print(f'{length_of_shared_matches_df} of the queries have multiple possible matches. Consensus sequences will be produced for them.')
    return shared_matches_df, unique_df, sm_unique_queries

### L2 # The curated output from vsearch is formatted ready for the output csv
def produce_output_df(curated_df, output_columns, region):
    curated_df = curated_df.rename(columns={'Q': 'Query#', 'Vs Name': 'DB Name'})
    curated_df['Region'] = region
    curated_df = add_selected_column(curated_df)
    curated_df = curated_df.sort_values(by=['Query#', 'Vs DB#']).reset_index(drop=True)
    curated_df.name = 'curated_df'
    check_df_contains_required_columns(curated_df, output_columns)
    output_df = curated_df[output_columns]
    return output_df

## L3 # Adds the 'Selected' column that can be later used to filter the output
def add_selected_column(df):
    grouped = df.groupby('Query#')
    df['Selected'] = 'x'
    for name, group in grouped:
        if len(group) == 1:
            df.loc[group.index, 'Selected'] = 'Y'  # Identifier when only one match created
        else:
            df.loc[group.index[0], 'Selected'] = 'Y.'  # Identifier for 1st row
            df.loc[group.index[1:], 'Selected'] = 'N'  # Identifier for subsequent (non-selected) rows
    return df

## L3 # A check to show whch columns are missing before formatting
def check_df_contains_required_columns(df, column_headers):
    missing_columns = set(column_headers) - set(df.columns)
    if missing_columns:
        print(f"The {df.name} can not be reordered as columns {missing_columns} are not present.")
    else:
        return

### L2 # Splits the formatted output_df. Shared_matches will be processed to produce consensus sequence, then concatenated with the unique_df later
def split_output_df(output_df):
    shared_matches_df = output_df[(output_df['Selected'] == 'Y.') | (output_df['Selected'] == 'N')]
    unique_df = output_df[output_df['Selected'] == 'Y']
    rump_df = output_df[~output_df.index.isin(shared_matches_df.index.union(unique_df.index))]
    length_of_shared_matches_df = shared_matches_df['Query#'].nunique()
    if len(rump_df) > 0:
        print('Not all of the rows from the output_df have been processed.')
    elif len(shared_matches_df) + len(unique_df) != len(output_df):
        print('Something has gone wrong in splitting the output_df')
    else:
        sm_unique_queries = list(set(shared_matches_df['Query#'].tolist()))
        un_unique_queries = list(set(unique_df['Query#'].tolist()))
        if any(val in sm_unique_queries for val in un_unique_queries):
            print("The output_df has not been split correctly")
        unique_df.name = 'Non-Shared Matches'
        check_df_column_unique_entries_only(unique_df, 'Query#')
    return shared_matches_df, unique_df, sm_unique_queries, length_of_shared_matches_df

###### L1 # Adds 3 sets of consensus sequences to shared_matches_df
def add_consensus_sequences(shared_matches_df, sm_unique_queries, barcode_columns, output_columns):
    pre_consensus_df = shared_matches_df.copy()
    strict_consensus_df = create_majority_or_strict_consensus_df(pre_consensus_df, shared_matches_df, barcode_columns,
                                                                 output_columns,
                                                                 minimum_threshold=0.99, selected_filter='Cs')
    majority_consensus_df = create_majority_or_strict_consensus_df(pre_consensus_df, shared_matches_df, barcode_columns,
                                                                   output_columns, minimum_threshold=0.5,
                                                                   selected_filter='Cm')
    plurality_consensus_df = create_plurality_consensus_df(pre_consensus_df, shared_matches_df, barcode_columns,
                                                           output_columns, selected_filter='Cp')
    shared_matches_with_consensus_df = concatenate_dfs(shared_matches_df, majority_consensus_df, strict_consensus_df,
                                                       plurality_consensus_df)
    shared_matches_with_consensus_df = shared_matches_with_consensus_df.sort_values(by='Query#')
    return shared_matches_with_consensus_df

### L2 # creates new dfs to be concatenated with the shared matches df
def create_majority_or_strict_consensus_df(pre_consensus_df, shared_matches_df, barcode_columns, output_columns,
                                           minimum_threshold, selected_filter):
    # pre_consensus_df = pre_consensus_df.copy()
    pre_consensus_df, consensus_df = add_barcodes(pre_consensus_df, barcode_columns, minimum_threshold)
    pre_consensus_df, matching_df = additional_consensus_barcode_filtering(pre_consensus_df, barcode_columns)
    consensus_df = concatenate_dfs(consensus_df, matching_df)
    consensus_df = fill_consensus_df_missing_values(pre_consensus_df, consensus_df, shared_matches_df, barcode_columns)
    check_consensus_df_correct(consensus_df, shared_matches_df)
    consensus_df = prepare_consensus_df_for_merge_with_shared_matches_df(consensus_df, barcode_columns, selected_filter,
                                                                         output_columns)
    return consensus_df

# L2 # alternative version of adding consensus sequences, used for the plurality
def create_plurality_consensus_df(pre_consensus_df, shared_matches_df, barcode_columns, output_columns,
                                  selected_filter):
    pre_consensus_df, consensus_df = add_plurality_barcodes(pre_consensus_df, barcode_columns)
    pre_consensus_df, matching_df = additional_consensus_barcode_filtering(pre_consensus_df, barcode_columns)
    consensus_df = concatenate_dfs(consensus_df, matching_df)
    consensus_df = fill_consensus_df_missing_values(pre_consensus_df, consensus_df, shared_matches_df, barcode_columns)
    check_consensus_df_correct(consensus_df, shared_matches_df)
    consensus_df = prepare_consensus_df_for_merge_with_shared_matches_df(consensus_df, barcode_columns, selected_filter,
                                                                         output_columns)
    return consensus_df

## L3 # Adds barcodes for majority and strict consensus dfs
def add_barcodes(df, columns, minimum_threshold):
    consensus_df = pd.DataFrame(columns=df.columns)  # initialise empty df for transfer of processed queries
    grouped = df.groupby('Query#')
    df_filtered = df.copy()
    # iterate over each column
    for column in columns:
        new_column = column.lower()
        df[new_column] = '_'  # new barcode column, to be created from old
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
                grouped = df_filtered.groupby(
                    'Query#')  # this line re-group the dataset for the next column iteration, delete?
            # if group < threshold, copy first_row to consensus_df. Then drop rows
            if not test_result:
                first_row = group.head(1).copy()
                first_row.loc[:, new_column] = '-'
                consensus_df = pd.concat([consensus_df, first_row])
                df = df[df['Query#'] != name]
        grouped = df_filtered.groupby('Query#')
    return df, consensus_df

## L3 # plurality version of add_majority_barcodes
def add_plurality_barcodes(df, columns):
    consensus_df = pd.DataFrame(columns=df.columns)  # initialise empty df for transfer of processed queries
    grouped = df.groupby('Query#')
    df_filtered = df.copy()
    for column in columns:
        new_column = column.lower()
        df[new_column] = '_'  # new barcode column, to be created from old
        for name, group in grouped:
            modal_value = group[column].mode()[0]
            frequency_of_modal_value = group[column].value_counts()[modal_value]
            group_size = len(group)
            number_of_values_mode_is_true = len([value for value in group[column].value_counts().index if
                                                 value in group[column].value_counts() and group[column].value_counts()[
                                                     value] == frequency_of_modal_value])
            if frequency_of_modal_value > 1 and number_of_values_mode_is_true == 1:
                df_filtered.loc[
                    (df_filtered['Query#'] == name) & (df_filtered[column] == modal_value), new_column] = modal_value
                df_filtered.loc[
                    (df_filtered['Query#'] == name) & (df_filtered[column] != modal_value), new_column] = 'x'
                df_filtered = df_filtered[df_filtered[new_column] != 'x']
            elif frequency_of_modal_value == 1:
                first_row = group.head(1).copy()
                first_row[new_column] = '-'
                consensus_df = pd.concat([consensus_df, first_row])
                df_filtered = df_filtered.drop(group.index)
            elif frequency_of_modal_value > 1 and number_of_values_mode_is_true > 1:
                first_row = group.head(1).copy()
                first_row[new_column] = '-'
                consensus_df = pd.concat([consensus_df, first_row])
                df_filtered = df_filtered.drop(group.index)
        grouped = df_filtered.groupby('Query#')
    return df_filtered, consensus_df

# L4 # calculates threshold for each 'Query#' group in each barcode column. Tests if group passes threshold
def calculate_threshold(group_size, frequency_of_modal_value, minimum_threshold):
    threshold = group_size * minimum_threshold
    test_result = frequency_of_modal_value / group_size > minimum_threshold
    return threshold, test_result

## L3 # where query does not end with unique barcode, this function checks all rows for each group match, then combines into one row. Combined rows are returned in matching_df for concatenation with consenus_df. pre_consensus_df is returned to check is empty
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
        consensus_df = consensus_df.copy()  # added in to prevent Setting with Copy Warning
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

# L4 # This is an extra function added because of bug in testing where a query could be processed twice, in the immediately following barcode column
# It checks for duplicated queries in the consensus_df, keeping only the first entry. Statement printed of any deletions.
def remove_duplicated_processed_queries(consensus_df):
    unique_df = consensus_df.drop_duplicates(subset='Query#')  # groups by unique Query# values
    duplicated_df = consensus_df[~consensus_df.index.isin(unique_df.index)]  # creates df of duplicated values
    if duplicated_df.empty == False:  # checks if empty
        rows_deleted = len(duplicated_df)
        duplicated_queries = list(set(duplicated_df['Query#']))
        print(
            f"{rows_deleted} rows have been deleted from the consensus_df. The following queries were duplicated: {duplicated_queries}.")
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

# L4 # When consensus barcode is added it stops when threshold failed. This function extends barcode ('-' values) to n column
def replace_nans(consensus_df, new_barcode_columns):
    consensus_df[new_barcode_columns] = consensus_df[new_barcode_columns].fillna('-')
    return consensus_df

## L3 # Checks consensus df to make sure that queries from input and output dfs match, and that output queries are not duplicated
def check_consensus_df_correct(consensus_df, shared_matches_df):
    consensus_df.name = 'consensus_df'
    shared_matches_df.name = 'shared_matches_df'
    check_unique_values(consensus_df, shared_matches_df, 'Query#')
    check_df_column_unique_entries_only(consensus_df, 'Query#')

## L3 # prepares consensus_df for merging with shared_matches_df
def prepare_consensus_df_for_merge_with_shared_matches_df(consensus_df, barcode_columns, selected_filter,
                                                          output_columns):
    consensus_df = consensus_df.drop(columns=barcode_columns)  # drops old barcode columns
    consensus_df = consensus_df.rename(
        columns={'a': 'A', 'b': 'B', 'c': 'C', 'd': 'D', 'e': 'E', 'f': 'F', 'g': 'G', 'h': 'H', 'i': 'I', 'j': 'J',
                 'k': 'K', 'l': 'L', 'm': 'M', 'n': 'N'})  # changes new barcode columns names to old
    consensus_df = consensus_df.sort_values(by='Query#')
    consensus_df['Selected'] = selected_filter
    consensus_df['Similarity(%)'] = 'N/A'
    consensus_df['Vs DB#'] = 'N/A'
    consensus_df['DB Name'] = 'N/A'
    consensus_df['DB Accession Number'] = 'N/A'
    check_df_contains_required_columns(consensus_df, output_columns)
    # curated_vsearch_output, database_df, query_df, summary_df, output_columns = read_inputs()
    if set(output_columns).issubset(set(consensus_df.columns)):
        consensus_df = consensus_df.reindex(output_columns, axis=1)
    else:
        print("Error: Not all columns in 'output_columns' list are present in 'consensus_df' DataFrame.")
    return consensus_df

##### L1
def combine_shared_matches_and_unique_dfs(shared_matches_with_consensus_df, unique_df, database_df):
    final_df = concatenate_dfs(shared_matches_with_consensus_df, unique_df)
    final_df.sort_values(by='Query#', ascending=True, inplace=True)
    final_df = add_detail_for_xls_file(final_df, database_df)
    return final_df

### L2
def add_detail_for_xls_file(final_df, database_df):
    final_df = add_tree_order_column(final_df, database_df)
    return final_df

def add_tree_order_column(final_df, database_df):
    tree_df = database_df[['DB#', 'Tree#']]
    final_df = final_df.merge(tree_df, left_on='Query#', right_on='DB#', how='left')
    final_df.insert(final_df.columns.get_loc('Query#')+3, 'Tree#', final_df.pop('Tree#'))
    final_df.drop(columns=['DB#'], inplace=True)
    return final_df


##### L1 # Move files and make sure things are where they need to be
def file_outputs(final_df, summary_df, output_filename, summary_filename, query_df, curated_vsearch_output, query_fastaname, database_fastaname, vsearch_output, input_sequences_list, input_database):
    length_of_output = output_checks(final_df, query_df, output_filename)
    tidy_files(curated_vsearch_output, query_fastaname, database_fastaname, vsearch_output, input_sequences_list, input_database)
    write_csv(final_df, output_filename)
    write_csv(summary_df, summary_filename)
    return length_of_output

### L2
def output_checks(final_df, query_df, output_filename):
    length_of_output = final_df['Query#'].nunique()
    final_df.name = os.path.splitext(output_filename)[0]
    check_unique_values(query_df, final_df, 'Query#')
    return length_of_output

### L2
def tidy_files(curated_vsearch_output, query_fastaname, database_fastaname, vsearch_output, input_sequences_list, input_database):
    os.mkdir('data')
    shutil.move(curated_vsearch_output, 'data')
    shutil.move(vsearch_output, 'data')
    shutil.move(query_fastaname, 'data')
    shutil.move(database_fastaname, 'data')
    shutil.move(input_sequences_list, 'data')
    shutil.move(input_database, 'data')

####### Helper functions

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
        # print(f"{df1.name} and {df2.name} are compatible. They contain {number_of_shared_unique_values} unique values in the '{column_header}' column.")
        return
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
        return

def check_groups_match_before_drop(df, col1, col2):
    groups = df.groupby(col1)   # group by 1st column
    all_groups_match = True # Initialize a variable to store col 2 info
    for group in groups:
        unique_values = group[1][col2].unique()
        if len(unique_values) != 1:
            all_groups_match = False
            break
    if all_groups_match:
        return
    else:
        print('There is a discrepancy in the grouped similarities for the summary_df')

def concatenate_dfs(*dfs):
    column_headers = [list(df.columns) for df in dfs]
    if all(headers == column_headers[0] for headers in column_headers):
        concatenated_df = pd.concat(dfs)
        return concatenated_df
    else:
        df_names = ', '.join([str(df) for df in dfs])
        print(f"{df_names} could not be concatenated as the column headers do not match. Please rectify")

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