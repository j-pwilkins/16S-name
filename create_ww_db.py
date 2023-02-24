import pandas as pd
import numpy as np
import os
import shutil
import sys
from datetime import datetime

def main():
    start_time = datetime.now()
    input_sequences_list, vsearch_output, curated_vsearch_output, output_directory, output_file, alt_matches_file, detail_database = set_variables()
    organise_directories(output_directory, input_sequences_list, detail_database, vsearch_output)
    input_df, length_of_input = prepare_and_run_vsearch(input_sequences_list)
    processed_vsearch_df, alternative_matches_df, AvA_df = process_vsearch_output(vsearch_output, curated_vsearch_output, input_df, length_of_input)
    barcoded_df = create_barcodes(processed_vsearch_df)
    barcoded_df = add_input_detail(barcoded_df, input_df)
    summary_df = add_database_detail(barcoded_df, detail_database)
    tidy_files(barcoded_df, summary_df, vsearch_output, input_sequences_list, AvA_df, alternative_matches_df, alt_matches_file, curated_vsearch_output)

def tidy_files(barcoded_df, summary_df, vsearch_output, input_sequences_list, AvA_df, alternative_matches_df, alt_matches_file, curated_vsearch_output):
    os.mkdir('data')
    shutil.move(vsearch_output, 'data')
    shutil.move(input_sequences_list, 'data')
    write_csv(AvA_df, curated_vsearch_output)
    shutil.move(curated_vsearch_output, 'data')
    # write_csv(alternative_matches_df, alt_matches_file)
    # shutil.move(alt_matches_file, 'data')
    write_csv(barcoded_df, 'barcoded_strains.csv')
    write_csv(summary_df, 'barcoded_database.csv')

def set_variables():
    # input_sequences_list = sys.argv[1]
    input_sequences_list = '10seqSample.fas'         # temp line as sys.argv not available
    vsearch_output = 'vsearch_output.csv'
    curated_vsearch_output = 'All_vs_All.csv'
    detail_database = 'NIRE_metadata_gisaid_2023_02_18.tsv'
    default_output_directory = input_sequences_list.rstrip('.csv')
    if len(sys.argv) > 2:
        output_directory = sys.argv[2]
    else:
        output_directory = default_output_directory
    output_file = 'Summary_' + output_directory + '.csv'
    alt_matches_file = 'Alternative_Matches_for_' + output_directory + '.csv'
    return input_sequences_list, vsearch_output, curated_vsearch_output, output_directory, output_file, alt_matches_file, detail_database

def organise_directories(output_directory, input_sequences_list, detail_database, vsearch_output):
    if os.path.exists(output_directory):    # temp lines for testing
        shutil.rmtree(output_directory)     # removes previous directory if it exists
        print(f'The previous {output_directory} folder has been deleted')
    os.makedirs(output_directory, exist_ok=True)
    print(f'{output_directory} has been created')
    shutil.copy(input_sequences_list, output_directory)
    shutil.copy(vsearch_output, output_directory)           # temp line to allow running
    shutil.copy(detail_database, output_directory)         # temp line to allow running
    os.chdir(output_directory)
    # return input_df, input_filename, number_of_sequences_entered

##### L1
def prepare_and_run_vsearch(input_sequences_list):
    # Read the fasta file into a string
    with open(input_sequences_list) as f:
        fasta_string = f.read()
    # Split the fasta string into a list of sequences
    input_sequences_list = fasta_string.split(">")[1:]
    # Call the function to create the DataFrame
    input_df, length_of_input = produce_input_df_from_fasta(input_sequences_list)
    # Add 'DB#' column
    input_df.insert(0, 'DB #', range(1, len(input_df) + 1))
    # Add column with the length of the sequence
    input_df.insert(input_df.columns.get_loc("Sequence"), "Length", input_df["Sequence"].apply(len))
    # run_vsearch(input_sequences_list)
    return input_df, length_of_input

#####L2 # This function takes a list of fasta sequences and returns a pandas DataFrame with two columns: 'Name' and 'Sequence'.
def produce_input_df_from_fasta(input_sequences_list):
    # Initialize empty lists to store headers and sequences
    headers = []
    sequences = []
    # Iterate over input sequences
    for sequence in input_sequences_list:
        lines = sequence.split("\n")    # Split the sequence into lines
        # Get the header (first line) and sequence (second line)
        header = lines[0][0:]
        sequence = "".join(lines[1:])
        # Add the header and sequence to the respective lists
        headers.append(header)
        sequences.append(sequence)
    # Create a dictionary to store the data
    data_dict = {"Name": headers, "Sequence": sequences}
    # Convert the dictionary to a DataFrame
    df = pd.DataFrame.from_dict(data_dict)
    length_of_input = len(df)
    print(f"vsearch will process {length_of_input} sequences")
    return df, length_of_input

### L2
def run_vsearch(input_sequences_list):
    cmd = 'vsearch --allpairs_global ' + input_sequences_list + ' --blast6out vsearch_output.csv --acceptall'
    os.system(cmd)
    if "cmd" not in locals():
        print('A pre-prepared vsearch output has been used here')

###### L1
def process_vsearch_output(vsearch_output, curated_vsearch_output, input_df, length_of_input):
    AvA_df = parse_vsearch_allpairsglobal_blast6(vsearch_output)
    check_length_of_parsed_output(AvA_df, length_of_input)
    similarity_df, alternative_matches_df = retain_only_highest_similarities(AvA_df)
    return similarity_df, alternative_matches_df, AvA_df

### L2 # Convert the vsearchoutput .csv file to a useable df
def parse_vsearch_allpairsglobal_blast6(vsearch_output):
    # Take required detail from vsearch output
    all_vs_all_forward_df = pd.read_csv(vsearch_output, usecols=[0, 1, 2], sep='\t',
                                        names=["Query Header", "Target Header", "Sim to Target"])
    # Double & Reverse the columns to include all pairwise comparisons
    all_vs_all_df = pd.concat([all_vs_all_forward_df, all_vs_all_forward_df.rename(
        columns={"Target Header": "Query Header", "Query Header": "Target Header"})], ignore_index=True)
    # Add column specifying vsearch order - helps troubleshoot
    all_vs_all_df.insert(0, "Vsearch Order", np.arange(1, len(all_vs_all_df) + 1), True)
    all_vs_all_df.insert(4, 'F/R', np.concatenate(
        [np.array(["F"] * len(all_vs_all_forward_df)), np.array(["R"] * len(all_vs_all_forward_df))]))
    # Split 'Query Header' and 'Target Header' column into 4 new columns
    all_vs_all_df[['Query #', 'Query Name']] = all_vs_all_df['Query Header'].str.extract('(Q\d+)\s+(.*)', expand=True)
    all_vs_all_df[['Vs #', 'Vs Name']] = all_vs_all_df['Target Header'].str.extract('(Q\d+)\s+(.*)', expand=True)
    all_vs_all_df['Query #'] = all_vs_all_df['Query #'].str.replace('Q', '') # Strip Q from # columns
    all_vs_all_df['Vs #'] = all_vs_all_df['Vs #'].str.replace('Q', '')
    all_vs_all_df['Query #'] = all_vs_all_df['Query #'].astype('int64')     # convert to int to be able to manipulate
    all_vs_all_df['Vs #'] = all_vs_all_df['Vs #'].astype('int64')
    # Reorder columns and drop 'Header' columns
    all_vs_all_df = all_vs_all_df[['Vsearch Order', 'Query #', 'Query Name', 'Vs #', 'Vs Name', 'F/R', 'Sim to Target']]
    return all_vs_all_df

### L2 # check that all of the sequences have been parsed with the correct number of targets
def check_length_of_parsed_output(AvA_df, length_of_input):
    expected_length = length_of_input * (length_of_input - 1)
    actual_length = len(AvA_df)
    if expected_length == actual_length:
        print(f'{length_of_input} sequences were processed by vsearch correctly')
    else:
        print(f"Parsed output length check failed. Expected length: {expected_length}. Actual length: {actual_length}.")

### L2 # Find highest similaruty for each query & any alternate matches
def retain_only_highest_similarities(curated_df):
    dummy_first_row = create_dummy_first_row(curated_df)
    # Only sequences lower in the input order are eligible for comparison - code drops later sequences
    filtered_df = curated_df[curated_df['Query #'] > curated_df['Vs #']]
    grouped = filtered_df.groupby('Query #')
    filtered_df = filtered_df.copy()
    # create new column that finds the highest similarity for each group
    filtered_df['Highest Similarity'] = grouped['Sim to Target'].transform(max)
    # create new column that indicates which sequences had alternative matching possibilities
    filtered_df['Alternative Matches'] = filtered_df.apply(lambda row: sum(
        (row['Highest Similarity'] == filtered_df['Sim to Target']) & (
                filtered_df['Query #'] == row['Query #'])) - 1, axis=1)
    # Only those rows with the highest similarity are kept
    filtered_df = filtered_df[filtered_df['Sim to Target'] == filtered_df['Highest Similarity']]
    # a new df containing all sequences with alternative matches is kept
    alternative_matches_df = filtered_df[filtered_df['Alternative Matches'] > 0]
    # Concatenate the dummy row with the filtered dataframe
    filtered_df = pd.concat([dummy_first_row, filtered_df])
    # Calculate and print the number of sequences processed and the number of sequences with alternate matches
    number_of_sequences_processed = len(filtered_df)
    number_of_sequences_with_alternate_matches = len(alternative_matches_df['Query #'].unique())
    print(f"{number_of_sequences_with_alternate_matches} of {number_of_sequences_processed} sequences have alternate matches")
    return filtered_df, alternative_matches_df

## L3 # Adds dummy 1st row for Query 1 which otherwise has no match
def create_dummy_first_row(curated_df):
    # Create a new dataframe with just the first row where Query # = 1
    dummy_first_row = curated_df.loc[curated_df['Query #'] == 1].iloc[:1].copy()
    # Format dummy_first_row ready for concatenation
    dummy_first_row.loc[:, 'Vs #'] = 'N/A'
    dummy_first_row.loc[:, 'Vs Name'] = 'N/A'
    dummy_first_row.loc[:, 'Sim to Target'] = 'N/A'
    dummy_first_row['Highest Similarity'] = 'N/A'
    dummy_first_row['Alternative Matches'] = 0
    return dummy_first_row

##### L1
def create_barcodes(processed_vsearch_df):
    processed_vsearch_df = processed_vsearch_df.rename(columns={'Query #': 'DB #', 'Vs #': 'Vs DB #', 'Sim to Target': 'Similarity(%)', 'Query Name': 'DB Name'})
    barcode_df = pd.DataFrame({'DB #': [], 'A': [], 'B': [], 'C': [], 'D': [], 'E': [], 'F': [], 'G': [], 'H': [], 'I': [], 'J': [], 'K': [], 'L': [], 'M': [], 'N': [], 'Vs': [], 'Similarity(%)': []})
    for i in range(len(processed_vsearch_df)):
        barcode_df = name_each_sequence(processed_vsearch_df, barcode_df, i)
    result, message = check_barcode_df(processed_vsearch_df, barcode_df)
    if result:
        barcoded_df = merge_vsearch_with_barcode(processed_vsearch_df, barcode_df)
        barcoded_df = add_tree_order_column(barcoded_df)
    else:
        print(message)
    return barcoded_df

### L2
def name_each_sequence(input_df, barcode_df, i, A=60, B=70, C=80, D=85, E=90, F=95, G=98, H=99, I=99.5, J=99.6, K=99.7, L=99.8, M=99.9, N=100):
    row = i + 1
    if row == 1:
        barcode_df = pd.DataFrame({'DB #': [row], 'A': [0], 'B': [0], 'C': [0], 'D': [0],
                                   'E': [0], 'F': [0], 'G': [0], 'H': [0], 'I': [0], 'J': [0],
                                   'K': [0], 'L': [0], 'M': [0], 'N': [0], 'Vs #': ['N/A'], 'Similarity(%)': ['N/A']})

    else:
        matched_df = input_df.loc[input_df['DB #'] == row]
        closest_current_sequence = int(matched_df['Vs DB #'].values[0])
        similarity_to_closest_sequence = float(matched_df['Similarity(%)'].values[0])

        Ai = int(barcode_df.loc[barcode_df['DB #'] == closest_current_sequence, 'A'].values[0])
        Bi = int(barcode_df.loc[barcode_df['DB #'] == closest_current_sequence, 'B'].values[0])
        Ci = int(barcode_df.loc[barcode_df['DB #'] == closest_current_sequence, 'C'].values[0])
        Di = int(barcode_df.loc[barcode_df['DB #'] == closest_current_sequence, 'D'].values[0])
        Ei = int(barcode_df.loc[barcode_df['DB #'] == closest_current_sequence, 'E'].values[0])
        Fi = int(barcode_df.loc[barcode_df['DB #'] == closest_current_sequence, 'F'].values[0])
        Gi = int(barcode_df.loc[barcode_df['DB #'] == closest_current_sequence, 'G'].values[0])
        Hi = int(barcode_df.loc[barcode_df['DB #'] == closest_current_sequence, 'H'].values[0])
        Ii = int(barcode_df.loc[barcode_df['DB #'] == closest_current_sequence, 'I'].values[0])
        Ji = int(barcode_df.loc[barcode_df['DB #'] == closest_current_sequence, 'J'].values[0])
        Ki = int(barcode_df.loc[barcode_df['DB #'] == closest_current_sequence, 'K'].values[0])
        Li = int(barcode_df.loc[barcode_df['DB #'] == closest_current_sequence, 'L'].values[0])
        Mi = int(barcode_df.loc[barcode_df['DB #'] == closest_current_sequence, 'M'].values[0])
        Ni = int(barcode_df.loc[barcode_df['DB #'] == closest_current_sequence, 'N'].values[0])

        if similarity_to_closest_sequence < A:
            Ao = int(max(barcode_df['A'].values)) + 1
            Bo = Co = Do = Eo = Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = 0
        else:
            Ao = Ai
            if similarity_to_closest_sequence < B:
                df_b = barcode_df[barcode_df['A'] == Ai]
                Bo = int(max(df_b['B'].values)) + 1
                Co = Do = Eo = Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = 0
            else:
                Bo = Bi
                if similarity_to_closest_sequence < C:
                    df_c = barcode_df[(barcode_df['A'] == Ai) & (barcode_df['B'] == Bi)]
                    Co = int(max(df_c['C'].values)) + 1
                    Do = Eo = Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = 0
                else:
                    Co = Ci
                    if similarity_to_closest_sequence < D:
                        df_d = barcode_df[(barcode_df['A'] == Ai) & (barcode_df['B'] == Bi) & (barcode_df['C'] == Ci)]
                        Do = int(max(df_d['D'].values)) + 1
                        Eo = Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = 0
                    else:
                        Do = Di
                        if similarity_to_closest_sequence < E:
                            df_e = barcode_df[
                                (barcode_df['A'] == Ai) & (barcode_df['B'] == Bi) & (barcode_df['C'] == Ci) & (
                                            barcode_df['D'] == Di)]
                            Eo = int(max(df_e['E'].values)) + 1
                            Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = 0
                        else:
                            Eo = Ei
                            if similarity_to_closest_sequence < F:
                                df_f = barcode_df[
                                    (barcode_df['A'] == Ai) & (barcode_df['B'] == Bi) & (barcode_df['C'] == Ci) & (
                                                barcode_df['D'] == Di) & (barcode_df['E'] == Ei)]
                                Fo = int(max(df_f['F'].values)) + 1
                                Go = Ho = Io = Jo = Ko = Lo = Mo = No = 0
                            else:
                                Fo = Fi
                                if similarity_to_closest_sequence < G:
                                    df_g = barcode_df[
                                        (barcode_df['A'] == Ai) & (barcode_df['B'] == Bi) & (barcode_df['C'] == Ci) & (
                                                barcode_df['D'] == Di) & (barcode_df['E'] == Ei) & (
                                                barcode_df['F'] == Fi)]
                                    Go = int(max(df_g['G'].values)) + 1
                                    Ho = Io = Jo = Ko = Lo = Mo = No = 0
                                else:
                                    Go = Gi
                                    if similarity_to_closest_sequence < H:
                                        df_h = barcode_df[(barcode_df['A'] == Ai) & (barcode_df['B'] == Bi) & (
                                                    barcode_df['C'] == Ci) & (barcode_df['D'] == Di) & (
                                                                      barcode_df['E'] == Ei) & (
                                                                      barcode_df['F'] == Fi) & (
                                                                      barcode_df['G'] == Gi)]
                                        Ho = int(max(df_h['H'].values)) + 1
                                        Io = Jo = Ko = Lo = Mo = No = 0
                                    else:
                                        Ho = Hi
                                        if similarity_to_closest_sequence < I:
                                            df_i = barcode_df[
                                                (barcode_df['A'] == Ai) & (barcode_df['B'] == Bi) & (
                                                            barcode_df['C'] == Ci) & (barcode_df['D'] == Di) & (
                                                            barcode_df['E'] == Ei) & (barcode_df['F'] == Fi) & (
                                                            barcode_df['G'] == Gi) & (barcode_df['H'] == Hi)]
                                            Io = int(max(df_i['I'].values)) + 1
                                            Jo = Ko = Lo = Mo = No = 0
                                        else:
                                            Io = Ii
                                            if similarity_to_closest_sequence < J:
                                                df_j = barcode_df[(barcode_df['A'] == Ai) & (barcode_df['B'] == Bi) & (
                                                            barcode_df['C'] == Ci) & (barcode_df['D'] == Di) & (
                                                                              barcode_df['E'] == Ei) & (
                                                                              barcode_df['F'] == Fi) & (
                                                                              barcode_df['G'] == Gi) & (
                                                                              barcode_df['H'] == Hi)]
                                                Jo = int(max(df_j['J'].values)) + 1
                                                Ko = Lo = Mo = No = 0
                                            else:
                                                Jo = Ji
                                                if similarity_to_closest_sequence < K:
                                                    df_k = barcode_df[
                                                        (barcode_df['A'] == Ai) & (barcode_df['B'] == Bi) & (
                                                                    barcode_df['C'] == Ci) & (barcode_df['D'] == Di) & (
                                                                    barcode_df['E'] == Ei) & (barcode_df['F'] == Fi) & (
                                                                    barcode_df['G'] == Gi) & (barcode_df['H'] == Hi) & (
                                                                    barcode_df['I'] == Ii)]
                                                    Ko = int(max(df_k['K'].values)) + 1
                                                    Lo = Mo = No = 0
                                                else:
                                                    Ko = Ki
                                                    if similarity_to_closest_sequence < L:
                                                        df_l = barcode_df[
                                                            (barcode_df['A'] == Ai) & (barcode_df['B'] == Bi) & (
                                                                        barcode_df['C'] == Ci) & (
                                                                        barcode_df['D'] == Di) & (
                                                                        barcode_df['E'] == Ei) & (
                                                                        barcode_df['F'] == Fi) & (
                                                                        barcode_df['G'] == Gi) & (
                                                                        barcode_df['H'] == Hi) & (
                                                                        barcode_df['I'] == Ii) & (
                                                                        barcode_df['J'] == Ji)]
                                                        Lo = int(max(df_l['L'].values)) + 1
                                                        Mo = No = 0
                                                    else:
                                                        Lo = Li
                                                        if similarity_to_closest_sequence < M:
                                                            df_m = barcode_df[
                                                                (barcode_df['A'] == Ai) & (barcode_df['B'] == Bi) & (
                                                                            barcode_df['C'] == Ci) & (
                                                                            barcode_df['D'] == Di) & (
                                                                            barcode_df['E'] == Ei) & (
                                                                            barcode_df['F'] == Fi) & (
                                                                            barcode_df['G'] == Gi) & (
                                                                            barcode_df['H'] == Hi) & (
                                                                            barcode_df['I'] == Ii) & (
                                                                            barcode_df['J'] == Ji) & (
                                                                            barcode_df['K'] == Ki)]
                                                            Mo = int(max(df_m['M'].values)) + 1
                                                            No = 0
                                                        else:
                                                            Mo = Mi
                                                            if similarity_to_closest_sequence < N:
                                                                df_m = barcode_df[(barcode_df['A'] == Ai) & (
                                                                            barcode_df['B'] == Bi) & (
                                                                                              barcode_df['C'] == Ci) & (
                                                                                              barcode_df['D'] == Di) & (
                                                                                              barcode_df['E'] == Ei) & (
                                                                                              barcode_df['F'] == Fi) & (
                                                                                              barcode_df['G'] == Gi) & (
                                                                                              barcode_df['H'] == Hi) & (
                                                                                              barcode_df['I'] == Ii) & (
                                                                                              barcode_df['J'] == Ji) & (
                                                                                              barcode_df['K'] == Ki) & (
                                                                                              barcode_df['L'] == Li)]
                                                                No = int(max(df_m['N'].values)) + 1
                                                            else:
                                                                No = Ni


        new_line = {'DB #': [row], 'A': [Ao], 'B': [Bo], 'C': [Co], 'D': [Do], 'E': [Eo], 'F': [Fo], 'G': [Go], 'H': [Ho],
                    'I': [Io], 'J': [Jo], 'K': [Ko], 'L': [Lo], 'M': [Mo], 'N': [No],
                    'Vs #': [closest_current_sequence], 'Similarity(%)': [similarity_to_closest_sequence]}
        new_line_df = pd.DataFrame(new_line)
        barcode_df = pd.concat([barcode_df, new_line_df], ignore_index=True, axis=0)
    return barcode_df

### L2
def check_barcode_df(processed_vsearch_df, barcode_df):
    # Check for unique values in 'DB #' column
    if processed_vsearch_df['DB #'].nunique() != processed_vsearch_df.shape[0]:
        return False, "processed_vsearch_df has duplicate values in 'DB #' column"
    if barcode_df['DB #'].nunique() != barcode_df.shape[0]:
        return False, "barcode_df has duplicate values in 'DB #' column"

    # Check for same unique values in 'DB#' column
    if set(processed_vsearch_df['DB #'].unique()) != set(barcode_df['DB #'].unique()):
        return False, "processed_vsearch_df and barcode_df have different values in 'DB #' column"

    # Check for same values in 'Similarity(%)' column
    for db_num in processed_vsearch_df['DB #'].unique():
        processed_similarity = processed_vsearch_df.loc[processed_vsearch_df['DB #'] == db_num, 'Similarity(%)'].iloc[0]
        barcode_similarity = barcode_df.loc[barcode_df['DB #'] == db_num, 'Similarity(%)'].iloc[0]
        if processed_similarity != barcode_similarity:
            return False, f"processed_vsearch_df and barcode_df have different values in 'Similarity(%)' for DB # {db_num}"

    return True, print("All checks passed meaning barcodes can be attached.")

### L2
def merge_vsearch_with_barcode(processed_vsearch_df, barcode_df):
    barcode_df = barcode_df.drop(columns=['Similarity(%)', 'Vs #'])
    merged_df = pd.merge(processed_vsearch_df, barcode_df, on='DB #', how='left')
    cols = ['DB #', 'DB Name', 'Alternative Matches',  'Similarity(%)', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'Vs DB #', 'Vs Name']
    return merged_df[cols]

def add_tree_order_column(barcoded_df):
    columns = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N']
    # sort the df rows using the barcode columns, from right to left
    barcoded_df.sort_values(by=columns, axis=0, ascending=True, inplace=True)
    # insert a new column with the Tree Order so that it can be selected by user
    barcoded_df.insert(1, 'Tree #', range(1, len(barcoded_df) + 1))
    # order df so that it will be returned in same order as it came in
    barcoded_df.sort_values(by='DB #', axis=0, ascending=True, inplace=True)
    return barcoded_df

##### L1 # merge input_df to barcoded_df to add input detail
def add_input_detail(barcoded_df, input_df):
    # Merge the two dataframes using the 'DB #' column
    merged_df = barcoded_df.merge(input_df, on='DB #')
    # Re-order columns
    merged_df = merged_df [['DB #', 'Tree #', 'DB Name', 'Alternative Matches', 'Similarity(%)', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'Vs DB #', 'Vs Name', 'Name', 'Length', 'Sequence']]
    # Drop the 'Name' column if 'DB Name' and 'Name' are the same
    merged_df.loc[merged_df['DB Name'] == merged_df['Name'], 'Name'] = None
    merged_df = merged_df.drop(columns=['Name'])
    return merged_df

def add_database_detail(barcoded_df, detail_database):
    # Read in detail database as a pandas DataFrame
    database_df = pd.read_csv(detail_database, delimiter='\t')
    # Merge the dataframes on 'DB Name' and 'strain' columns
    merged_df = pd.merge(barcoded_df, database_df, how='left', left_on='DB Name', right_on='strain')
    # Reorder columns when 'DB Name' and 'strain' match
    if all(merged_df['DB Name'] == merged_df['strain']):
        cols_to_reorder = ['DB #', 'Tree #', 'DB Name', 'Alternative Matches', 'Similarity(%)',
                           'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
                           'Vs DB #', 'Vs Name', 'Last vaccinated', 'Passage details/history',
                           'Type', 'gisaid_epi_isl', 'date', 'Additional location information',
                           'Sequence length', 'Host', 'Patient age', 'Gender', 'Clade', 'Pango lineage',
                           'Pango version', 'Variant', 'AA Substitutions', 'Submission date',
                           'Is reference?', 'Is complete?', 'Is high coverage?', 'Is low coverage?',
                           'N-Content', 'GC-Content', 'region', 'country', 'division', 'location',
                           'Length', 'Sequence']
        merged_df = merged_df[cols_to_reorder]
    return merged_df


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
print('finished')