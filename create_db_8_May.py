import pandas as pd
import numpy as np
import os
import shutil
import sys
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    start_time = datetime.now()
    input_sequences_list, vsearch_output, curated_vsearch_output, output_directory, output_file, alt_matches_file, category_thresholds = set_variables()
    input_df, input_filename, number_of_sequences_entered = organise_directories(output_directory, input_sequences_list,
                                                                                 output_file, vsearch_output)
    fastaname = prepare_and_run_vsearch(input_df, input_filename, number_of_sequences_entered)
    processed_vsearch_df, alternative_matches_df, AvA_df = process_vsearch_output(vsearch_output,
                                                                                  curated_vsearch_output, input_df)
    barcoded_df = create_barcodes(processed_vsearch_df, category_thresholds, input_sequences_list)
    tidy_files(fastaname, vsearch_output, input_sequences_list, barcoded_df, alternative_matches_df, AvA_df,
               output_file, alt_matches_file, curated_vsearch_output)
    summary_statement(output_directory, start_time)

##### L1
def set_variables():
    # input_sequences_list = sys.argv[1]
    input_sequences_list = 'Test10.csv'         # temp line
    vsearch_output = 'vsearch_output.csv'
    curated_vsearch_output = 'All_vs_All.csv'
    default_output_directory = input_sequences_list.rstrip('.csv')
    if len(sys.argv) > 2:
        output_directory = sys.argv[2]
    else:
        output_directory = default_output_directory
    output_file = 'Summary_' + output_directory + '.csv'
    alt_matches_file = 'Alternative_Matches_for_' + output_directory + '.csv'
    category_thresholds = {
        'A': 74.200,
        'B': 80.700,
        'C': 87.100,
        'D': 91.700,
        'E': 94.900,
        'F': 96.800,
        'G': 98.400,
        'H': 99.100,
        'I': 99.500,
        'J': 99.700,
        'K': 99.810,
        'L': 99.880,
        'M': 99.940,
        'N': 100
    }
    return input_sequences_list, vsearch_output, curated_vsearch_output, output_directory, output_file, alt_matches_file, category_thresholds


##### L1
def organise_directories(output_directory, input_sequences_list, output_file, vsearch_output):
    if os.path.exists(output_directory):  # temp lines for testing
        shutil.rmtree(output_directory)  # temp
        print(f'The previous {output_directory} folder has been deleted')
    os.makedirs(output_directory, exist_ok=True)
    shutil.copy(input_sequences_list, output_directory)
    shutil.copy(vsearch_output, output_directory)           # temp line
    os.chdir(output_directory)
    input_df, number_of_sequences_entered = import_info(input_sequences_list)
    check_input(input_df)
    input_filename = input_sequences_list.rsplit(".", 1)[0]
    print(f'{input_filename} contains {number_of_sequences_entered} rows and has no missing values.')
    return input_df, input_filename, number_of_sequences_entered


### L2  # Imports the input as a df, adding a column to show order of entry
def import_info(input_sequences_list):
    input_df = read_csv(input_sequences_list)
    # Add entry order column at the start of df
    input_df.insert(0, 'Entry Order', range(1, len(input_df) + 1))
    # Add column giving the length of each 16S sequence - directly before sequence column
    input_df.insert(input_df.columns.get_loc("16S rRNA sequence"), "Length", input_df["16S rRNA sequence"].apply(len))
    number_of_sequences_entered = len(input_df)
    return input_df, number_of_sequences_entered


### L2  # Function that checks the input file for missing sequences
def check_input(input_df):
    missing_values = input_df.isnull().sum()
    if missing_values.sum() > 0:
        print("The following columns contain missing values:")
        print(missing_values[missing_values > 0].index.tolist())
        exit
    else:
        return


##### L1
def prepare_and_run_vsearch(input_df, input_filename, number_of_sequences_entered):
    fastaname = input_filename + '.fa'
    for i in range(number_of_sequences_entered):
        convert_input_to_fasta(i, fastaname, input_df)
    # run_vsearch(fastaname)
    return fastaname


### L2
def convert_input_to_fasta(i, fastaname, input_df):
    firstline = '>Q' + str(input_df['Entry Order'].values[i]) + " " + str(
        input_df['DB Accession Number'].values[i]) + " " + str(input_df['DB Name'].values[i])
    secondline = str(input_df['16S rRNA sequence'].values[i])
    with open(fastaname, 'a') as f:
        f.write(f"{firstline}\n{secondline}\n")


### L2
def run_vsearch(fastaname):
    cmd = 'vsearch --allpairs_global ' + fastaname + ' --blast6out vsearch_output.csv --acceptall'
    os.system(cmd)
    if "cmd" not in locals():
        print('A pre-prepared vsearch output has been used here')


##### L1
def process_vsearch_output(vsearch_output, curated_vsearch_output, input_df):
    AvA_df = parse_vsearch_allpairsglobal_blast6(vsearch_output)
    # write_csv(AvA_df, 'AvA.csv')      # temp line to test 'Alternate Matches'
    # AvA_df = read_csv('AvA.csv')        # temp line to test 'Alternate Matches'
    # function that trims the vvsearch df to leave only eligible highest similarity matches for each inputted sequence
    # Where multiple matches exist all alternates moved to separate df. Lowest alternate input only chosen by default
    similarity_df, alternative_matches_df = retain_only_highest_similarities(AvA_df)
    # Detail from input .csv file integrated with vsearch output, ready for program output
    similarity_df = merge_detail_from_input(similarity_df, input_df)
    alternative_matches_df = merge_detail_from_input(alternative_matches_df, input_df)
    # fill in nan values for first DB # which has no comparison sequence as it is first
    similarity_df.loc[0, ['Similarity (%)', 'Alternative Matches', 'Vs DB #', 'Vs Name', 'Vs ID']] = similarity_df.loc[
        0, ['Similarity (%)', 'Alternative Matches', 'Vs DB #', 'Vs Name', 'Vs ID']].fillna('N/A')
    # This line is necessary because merge_detail produces entries for every sequence in the input -
    # could maybe be removed at some stage if merge_detail is fixed, though I think there is a benefit for similarity_df
    # Line keeps only those sequences with alternate matches
    alternative_matches_df = alternative_matches_df[alternative_matches_df['Alternative Matches'] > 0]
    # produce summary statement - allows visual check of data
    alternative_matches = alternative_matches_df['DB #'].unique()
    print(
        f'{len(similarity_df)} sequences were processed. {len(alternative_matches)} sequences had alternative matches.')
    return similarity_df, alternative_matches_df, AvA_df


### L2
def parse_vsearch_allpairsglobal_blast6(vsearch_output):
    all_vs_all_forward_df = pd.read_csv(vsearch_output, usecols=[0, 1, 2], sep='\t',
                                        names=["Query Seq", "Target Seq", "Sim to Target"])
    all_vs_all_df = pd.concat([all_vs_all_forward_df, all_vs_all_forward_df.rename(
        columns={"Target Seq": "Query Seq", "Query Seq": "Target Seq"})], ignore_index=True)
    all_vs_all_df.insert(0, "Vsearch Order", np.arange(1, len(all_vs_all_df) + 1), True)
    all_vs_all_df.insert(4, 'F/R', np.concatenate(
        [np.array(["F"] * len(all_vs_all_forward_df)), np.array(["R"] * len(all_vs_all_forward_df))]))
    all_vs_all_df[['Q', 'T']] = all_vs_all_df[['Query Seq', 'Target Seq']].apply(lambda x: x.str[1:].astype(int))
    all_vs_all_df = all_vs_all_df.sort_values(['Q', 'T'])
    return all_vs_all_df


### L2
def retain_only_highest_similarities(curated_df):
    # curated_df = read_csv('AvA.csv')    # temp line to test 'ALternative Matches'
    # Only sequences lower in the input order are eligible for comparison - code drops later sequences
    filtered_df = curated_df[curated_df['Q'] > curated_df['T']]
    grouped = filtered_df.groupby('Query Seq')
    filtered_df = filtered_df.copy()
    # create new column that finds the highest similarity for each group
    filtered_df['Highest Similarity'] = grouped['Sim to Target'].transform(max)
    # create new column that indicates which sequences had alternative matching possibilities
    filtered_df['Alternative Matches'] = filtered_df.apply(lambda row: sum(
        (row['Highest Similarity'] == filtered_df['Sim to Target']) & (
                filtered_df['Query Seq'] == row['Query Seq'])) - 1, axis=1)
    # Only those rows with the highest similarity are kept
    filtered_df = filtered_df[filtered_df['Sim to Target'] == filtered_df['Highest Similarity']]
    # a new df containing all sequences with alternative matches is kept
    alternative_matches_df = filtered_df[filtered_df['Alternative Matches'] > 0]
    # only the first alternate match is selected (by default the lowest target sequence)
    filtered_df = filtered_df.groupby('Query Seq').first().reset_index()
    return filtered_df, alternative_matches_df


### L2
def merge_detail_from_input(vsearch_df, input_df):
    vsearch_detail_df = vsearch_df[['Q', 'Sim to Target', 'Alternative Matches', 'T']]
    input_detail_df = pd.merge(input_df, vsearch_detail_df, left_on='Entry Order', right_on='Q', how='left')
    versus_df = input_df[['Entry Order', 'DB Name', 'DB Accession Number']]
    versus_df = versus_df.rename(
        columns={'Entry Order': 'Vs DB #', 'DB Name': 'Vs Name', 'DB Accession Number': 'Vs ID', 'T': 'Vs DB #'})
    input_detail_df = pd.merge(input_detail_df, versus_df, left_on='T', right_on='Vs DB #', how='left')
    input_detail_df = input_detail_df.rename(columns={'Entry Order': 'DB #', 'Sim to Target': 'Similarity (%)'})
    selected_columns = ['DB #', 'Hun#', 'Fas#', 'DB Accession Number', 'DB Name', 'Similarity (%)',
                        'Alternative Matches', 'Vs DB #', 'Vs Name', 'Vs ID', 'DB Found', 'Date Accessed',
                        'Length', '16S rRNA sequence']
    return input_detail_df[selected_columns]


##### L1
def create_barcodes(processed_vsearch_df, category_thresholds, input_sequences_list):
    barcode_df = pd.DataFrame(
        {'DB #': [], 'A': [], 'B': [], 'C': [], 'D': [], 'E': [], 'F': [], 'G': [], 'H': [], 'I': [], 'J': [], 'K': [],
         'L': [], 'M': [], 'N': [], 'Vs': [], 'Similarity (%)': []})
    barcode_columns = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N']
    for i in range(len(processed_vsearch_df)):
        barcode_df = name_each_sequence(processed_vsearch_df, barcode_df, i, category_thresholds)
    result, message = check_barcode_df(processed_vsearch_df, barcode_df)
    if result:
        barcoded_df = merge_vsearch_with_barcode(processed_vsearch_df, barcode_df)
        barcoded_df = add_tree_order_column(barcoded_df, barcode_columns)
        barcoded_df = add_distance_score_column(barcoded_df, input_sequences_list, category_thresholds, barcode_columns)
    else:
        print(message)
    return barcoded_df


### L2
def name_each_sequence(input_df, barcode_df, i, category_thresholds):
    # Unpack category thresholds
    A, B, C, D, E, F, G, H, I, J, K, L, M, N = category_thresholds.values()
    row = i + 1
    if row == 1:
        barcode_df = pd.DataFrame({'DB #': [row], 'NT Diff Raw':[0], 'Cat Changed': ['A'], 'A': [0], 'B': [0], 'C': [0], 'D': [0],
                                   'E': [0], 'F': [0], 'G': [0], 'H': [0], 'I': [0], 'J': [0],
                                   'K': [0], 'L': [0], 'M': [0], 'N': [0], 'Vs': ['N/A'], 'Similarity (%)': ['N/A']})

    else:
        matched_df = input_df.loc[input_df['DB #'] == row]
        closest_current_sequence = int(matched_df['Vs DB #'].values[0])
        similarity_to_closest_sequence = float(matched_df['Similarity (%)'].values[0])

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
            category_changed = 'A'
            nt_diff = (100 - A)
        else:
            Ao = Ai
            if similarity_to_closest_sequence < B:
                df_b = barcode_df[barcode_df['A'] == Ai]
                Bo = int(max(df_b['B'].values)) + 1
                Co = Do = Eo = Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = 0
                category_changed = 'B'
                nt_diff = (100 - B)
            else:
                Bo = Bi
                if similarity_to_closest_sequence < C:
                    df_c = barcode_df[(barcode_df['A'] == Ai) & (barcode_df['B'] == Bi)]
                    Co = int(max(df_c['C'].values)) + 1
                    Do = Eo = Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = 0
                    category_changed = 'C'
                    nt_diff = (100 - C)
                else:
                    Co = Ci
                    if similarity_to_closest_sequence < D:
                        df_d = barcode_df[(barcode_df['A'] == Ai) & (barcode_df['B'] == Bi) & (barcode_df['C'] == Ci)]
                        Do = int(max(df_d['D'].values)) + 1
                        Eo = Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = 0
                        category_changed = 'D'
                        nt_diff = (100 - D)
                    else:
                        Do = Di
                        if similarity_to_closest_sequence < E:
                            df_e = barcode_df[
                                (barcode_df['A'] == Ai) & (barcode_df['B'] == Bi) & (barcode_df['C'] == Ci) & (
                                        barcode_df['D'] == Di)]
                            Eo = int(max(df_e['E'].values)) + 1
                            Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = 0
                            category_changed = 'E'
                            nt_diff = (100 - E)
                        else:
                            Eo = Ei
                            if similarity_to_closest_sequence < F:
                                df_f = barcode_df[
                                    (barcode_df['A'] == Ai) & (barcode_df['B'] == Bi) & (barcode_df['C'] == Ci) & (
                                            barcode_df['D'] == Di) & (barcode_df['E'] == Ei)]
                                Fo = int(max(df_f['F'].values)) + 1
                                Go = Ho = Io = Jo = Ko = Lo = Mo = No = 0
                                category_changed = 'F'
                                nt_diff = (100 - F)
                            else:
                                Fo = Fi
                                if similarity_to_closest_sequence < G:
                                    df_g = barcode_df[
                                        (barcode_df['A'] == Ai) & (barcode_df['B'] == Bi) & (barcode_df['C'] == Ci) & (
                                                barcode_df['D'] == Di) & (barcode_df['E'] == Ei) & (
                                                barcode_df['F'] == Fi)]
                                    Go = int(max(df_g['G'].values)) + 1
                                    Ho = Io = Jo = Ko = Lo = Mo = No = 0
                                    category_changed = 'G'
                                    nt_diff = (100 - G)
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
                                        category_changed = 'H'
                                        nt_diff = (100 - H)
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
                                            category_changed = 'I'
                                            nt_diff = (100 - I)
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
                                                category_changed = 'J'
                                                nt_diff = (100 - J)
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
                                                    category_changed = 'K'
                                                    nt_diff = (100 - A)
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
                                                        category_changed = 'L'
                                                        nt_diff = (100 - A)
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
                                                            category_changed = 'M'
                                                            nt_diff = (100 - A)
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
                                                                category_changed = 'N'
                                                                nt_diff = (100 - A)
                                                            else:
                                                                No = Ni
                                                                category_changed = '='
                                                                nt_diff = (100 - 100)

        new_line = {'DB #': [row], 'NT Diff Raw': [nt_diff], 'Cat Changed': [category_changed], 'A': [Ao], 'B': [Bo],
                    'C': [Co], 'D': [Do], 'E': [Eo], 'F': [Fo], 'G': [Go],'H': [Ho],
                    'I': [Io], 'J': [Jo], 'K': [Ko], 'L': [Lo], 'M': [Mo], 'N': [No],
                    'Vs': [closest_current_sequence], 'Similarity (%)': [similarity_to_closest_sequence]}
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

    # Check for same unique values in 'DB #' column
    if set(processed_vsearch_df['DB #'].unique()) != set(barcode_df['DB #'].unique()):
        return False, "processed_vsearch_df and barcode_df have different values in 'DB #' column"

    # Check for same values in 'Similarity (%)' column
    for db_num in processed_vsearch_df['DB #'].unique():
        processed_similarity = processed_vsearch_df.loc[processed_vsearch_df['DB #'] == db_num, 'Similarity (%)'].iloc[0]
        barcode_similarity = barcode_df.loc[barcode_df['DB #'] == db_num, 'Similarity (%)'].iloc[0]
        if processed_similarity != barcode_similarity:
            return False, f"processed_vsearch_df and barcode_df have different values in 'Similarity (%)' for DB # {db_num}"

    return True, print("All checks passed meaning barcodes can be attached.")


### L2
def merge_vsearch_with_barcode(processed_vsearch_df, barcode_df):
    barcode_df = barcode_df.drop(columns=['Similarity (%)', 'Vs'])
    merged_df = pd.merge(processed_vsearch_df, barcode_df, on='DB #', how='left')
    cols = ['DB #', 'Hun#', 'Fas#', 'DB Accession Number', 'DB Name', 'Alternative Matches', 'Similarity (%)', 'NT Diff Raw', 'Cat Changed', 'A', 'B',
            'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'Vs DB #', 'Vs Name', 'Vs ID', 'DB Found',
            'Date Accessed', 'Length', '16S rRNA sequence']
    write_csv(merged_df, 'mergetest.csv')
    return merged_df[cols]


def add_tree_order_column(barcoded_df, barcode_columns):
    # sort the df rows using the barcode columns, from right to left
    barcoded_df.sort_values(by=barcode_columns, axis=0, ascending=True, inplace=True)
    # insert a new column with the Tree Order so that it can be selected by user
    barcoded_df.insert(3, 'Tree#', range(1, len(barcoded_df) + 1))
    # order df so that it will be returned in same order as it came in
    barcoded_df.sort_values(by='DB #', axis=0, ascending=True, inplace=True)
    return barcoded_df

### L2
def add_distance_score_column(barcoded_df, input_sequences_list, category_thresholds, barcode_columns):
    barcoded_df = above_column(barcoded_df, barcode_columns)
    barcoded_df = below_column(barcoded_df, barcode_columns)
    barcoded_df = nearest_score(barcoded_df)
    barcoded_df = combine_scoring_columns(barcoded_df, input_sequences_list, category_thresholds, barcode_columns)
    write_csv(barcoded_df,'distance.csv')
    barcoded_df = add_nt_difference_columns(barcoded_df)
    return barcoded_df

## L3 # Adds column showing relatedness to the column above
def above_column(barcoded_df, barcode_columns):
    scoring_cols = barcode_columns
    scores = [15] + [15 - i for i in range(1, len(scoring_cols) - 1)] + [2, 1]  # scores for each column

    above_col = []  # initialize new column with empty values
    prev_row = None  # initialize previous row
    for i, row in barcoded_df.iterrows():
        if prev_row is None:
            above_col.append(15)  # first row gets a score of 15
        else:
            score = 0
            for j, col in enumerate(scoring_cols):
                if row[col] != prev_row[col]:
                    score = scores[j]
                    break
            if score == 0:
                score = 1  # if all columns match, assign score 1
            above_col.append(score)
        prev_row = row
    barcoded_df.insert(barcoded_df.columns.get_loc('N') + 1, 'Above Score', above_col)  # insert new column after 'N'
    return barcoded_df

## L3 # Adds column showing relatedness to the column below
def below_column(barcoded_df, barcode_columns):
    scoring_cols = barcode_columns
    scores = [15] + [15 - i for i in range(1, len(scoring_cols) - 1)] + [2, 1]  # scores for each column

    below_col = [15] * len(barcoded_df)  # initialize new column with 15
    next_row = None  # initialize next row
    for i in range(len(barcoded_df) - 1, -1, -1):
        row = barcoded_df.iloc[i]
        if next_row is None:
            below_col[i] = 15  # last row gets a score of 15
        else:
            score = 0
            for j, col in enumerate(scoring_cols):
                if row[col] != next_row[col]:
                    score = scores[j]
                    break
            if score == 0:
                score = 1  # if all columns match, assign score 1
            below_col[i] = score
        next_row = row

    barcoded_df.insert(barcoded_df.columns.get_loc('Above Score') + 1, 'Below Score', below_col)  # insert new column after 'Above Score'
    return barcoded_df

# L3 # Combines the above & below columns, then drops them
def nearest_score(barcoded_df):
    below_score_col = barcoded_df['Below Score'].tolist()
    above_score_col = barcoded_df['Above Score'].tolist()
    nearest_col = [min(x, y) for x, y in zip(below_score_col, above_score_col)]
    barcoded_df.insert(barcoded_df.columns.get_loc('Below Score') + 1, 'Nearest Score', nearest_col)  # insert new column after 'Below Score'
    barcoded_df = barcoded_df.drop(
        columns=['Below Score', 'Above Score'])  # drop 'Below Score' and 'Above Score' columns
    return barcoded_df

# L3 # Creates a Disctance Score column from Nearest Score column - based on % instead of ordinal
def combine_scoring_columns(barcoded_df, input_sequences_list, category_thresholds, barcode_columns):
    # Create a new column 'Distance Score' based on the 'Nearest Score' column
    reverse_ordinal_list = [15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
    scoring_list = [round(((100 - value)/100), 5) for value in category_thresholds.values()]
    scoring_list.insert(0, 0.005)
    barcoded_df.insert(barcoded_df.columns.get_loc('N') + 1, 'Distance Score', barcoded_df['Nearest Score'].replace(reverse_ordinal_list, scoring_list))
    write_csv(barcoded_df, 'errortest.csv')
    # Calculate how many entries are duplicates
    duplicated_entries = (barcoded_df['Similarity (%)'] == 100).sum()
    unique_entries = len(barcoded_df) - duplicated_entries
    print(f'{duplicated_entries} sequences take the code of another sequence. {unique_entries} unique codes were created.')
    # Sort the dataframe and create bar charts
    barcoded_df.sort_values(by='DB #', axis=0, ascending=True, inplace=True)
    create_similarity_bar_chart(barcoded_df, input_sequences_list, scoring_list)
    create_category_bar_chart(barcoded_df, input_sequences_list, barcode_columns)
    create_category_aggregate_bar_chart(barcoded_df, input_sequences_list, barcode_columns)
    return barcoded_df

### L4
def create_similarity_bar_chart(barcoded_df, input_sequences_list, scoring_list):
    # Filter out rows with Similarity (%) = 100
    filtered_df = barcoded_df[barcoded_df['Similarity (%)'] != 100]
    # Define categories and counts
    categories = scoring_list
    counts = [filtered_df['Distance Score'].value_counts()[cat] if cat in filtered_df['Distance Score'].values else 0 for cat in categories]
    total = len(filtered_df)  # get total number of sequences
    # Calculate percentage frequency for each category
    percentages = [count/total*100 for count in counts]
    # Create a bar chart
    plt.figure(figsize=(12, 6))  # set figure size to 12x6 inches
    plt.bar(range(len(categories)), percentages, width=0.6, color='green')
    # Add labels and title
    plt.xticks(range(len(categories)), categories)
    plt.xlabel('Similarity Distance Scores')
    plt.ylabel('% Frequency')
    plt.title(f'Frequency of {input_sequences_list} Covid Full Genome Similarity Distance Scores ({len(barcoded_df)-len(filtered_df)} excluded)')
    plt.ylim([0, 30])
    # Save the chart as a PNG file
    plt.savefig(f'Similarity Bar Chart_{input_sequences_list}.png')

### L4
def create_category_bar_chart(barcoded_df, input_sequences_list, barcode_columns):
    filtered_df = barcoded_df
    # Define categories and counts
    categories = barcode_columns + ['=']
    counts = [filtered_df['Cat Changed'].value_counts()[cat] if cat in filtered_df['Cat Changed'].values else 0 for cat in categories]
    total = len(filtered_df)  # get total number of sequences
    # Calculate percentage frequency for each category
    percentages = [count/total*100 for count in counts]
    # Create a bar chart
    plt.figure(figsize=(12, 6))  # set figure size to 12x6 inches
    plt.bar(range(len(categories)), percentages, width=0.6, color='blue')
    # Add labels and title
    plt.xticks(range(len(categories)), categories)
    plt.xlabel('Barcode Category Changes')
    plt.ylabel('% Frequency')
    plt.title(f'Frequency of {input_sequences_list} Covid Full Genome Similarity Distance Scores ({len(barcoded_df)-len(filtered_df)} excluded)')
    plt.ylim([0, 30])
    # Save the chart as a PNG file
    plt.savefig(f'Category Bar Chart_{input_sequences_list}.png')

### L4
def create_category_aggregate_bar_chart(barcoded_df, input_sequences_list, barcode_columns):
    filtered_df = barcoded_df
    # Define categories and counts
    categories = barcode_columns + ['=']
    counts = [filtered_df['Cat Changed'].value_counts()[cat] if cat in filtered_df['Cat Changed'].values else 0 for cat in categories]
    cumulative_counts = np.cumsum(counts) # Get cumulative sum of counts
    total_sequences = len(filtered_df)
    # Create a bar chart
    plt.figure(figsize=(12, 6))  # set figure size to 12x6 inches
    plt.bar(range(len(categories)), cumulative_counts, width=0.6, color='lightblue')
    # Add labels and title
    plt.xticks(range(len(categories)), categories)
    plt.xlabel('Accumulated Barcode Category Changes')
    plt.ylabel('Accumulated Frequency')
    plt.title(f'Accumulated Frequency of {input_sequences_list} Covid Full Genome Similarity Distance Scores ({len(barcoded_df)-len(filtered_df)} excluded)')
    plt.ylim([0, total_sequences])
    # Save the chart as a PNG file
    plt.savefig(f'Category Aggregate Bar Chart_{input_sequences_list}.png')

## L3
def add_nt_difference_columns(barcoded_df):
    write_csv(barcoded_df, 'btest1.csv')
    # convert column to show nucleotide distance from comparison sequence (i.e. nearest below in db)
    barcoded_df['Diff to Comparison (NT)'] = ((barcoded_df['NT Diff Raw'] / 100) * barcoded_df['Length']).round(0)
    # add new column that shows nucleotide distance from nearest db sequence (i.e. above/below in tree #)
    barcoded_df['Diff to Nearest (NT)'] = (barcoded_df['Distance Score'] * barcoded_df['Length']).round(0)
    # # reorder df to account for new cols
    cols = ['DB #', 'Hun#', 'Fas#', 'Tree#', 'DB Accession Number', 'DB Name', 'Alternative Matches', 'Similarity (%)',
            'Cat Changed', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'Distance Score',
            'Nearest Score', 'Diff to Comparison (NT)', 'Diff to Nearest (NT)', 'Vs DB #', 'Vs Name', 'Vs ID',
            'DB Found', 'Date Accessed', 'Length', '16S rRNA sequence']
    barcoded_df = barcoded_df[cols]
    write_csv(barcoded_df, 'btest2.csv')
    return barcoded_df


##### L1
def tidy_files(fastaname, vsearch_output, input_sequences_list, barcoded_df, alternative_matches_df, AvA_df,
               output_file, alt_matches_file, curated_vsearch_output):
    os.mkdir('data')
    shutil.move(vsearch_output, 'data')
    shutil.move(fastaname, 'data')
    shutil.move(input_sequences_list, 'data')
    write_csv(barcoded_df, output_file)
    write_csv(alternative_matches_df, alt_matches_file)
    shutil.move(alt_matches_file, 'data')
    write_csv(AvA_df, curated_vsearch_output)
    shutil.move(curated_vsearch_output, 'data')


##### L1
def summary_statement(output_directory, start_time):
    program_length = datetime.now() - start_time
    if len(sys.argv) > 1:
        print(
            f'{sys.argv[0]} has completed in {program_length}. {sys.argv[1]} was processed and the Summary.csv file can be '
            f'found in the new {output_directory} directory.')
    else:
        print(f'Test program has completed. Test processed, and Summary.csv file can be '
            f'found in the new {output_directory} directory.')

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
print('finished')
