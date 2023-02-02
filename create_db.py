import pandas as pd
import numpy as np
import os
import shutil
import sys
from datetime import datetime

def main():
    start_time = datetime.now()
    input_sequences_list = sys.argv[1]
    # input_sequences_list = 'test50.csv'         # temp line
    vsearch_output = 'vsearch_output.csv'
    curated_vsearch_output = 'All_vs_All.csv'
    default_output_directory = input_sequences_list.rstrip('.csv')
    if len(sys.argv) > 2:
        output_directory = sys.argv[2]
    else:
        output_directory = default_output_directory
    output_directory = default_output_directory[:]  # temp line
    output_file = 'Summary_' + output_directory + '.csv'
    if os.path.exists(output_directory):    # temp lines for testing
        shutil.rmtree(output_directory)     # temp
        print(f'The previous {output_directory} folder has been deleted')
    os.makedirs(output_directory, exist_ok=True)
    shutil.copy(input_sequences_list, output_directory)
    # shutil.copy(vsearch_output, output_directory)           # temp line
    os.chdir(output_directory)
    import_csv(input_sequences_list, output_file)
    fastaname, input_csv_as_df, length_of_csv = import_csv(input_sequences_list, output_file)
    for i in range(length_of_csv):
        convert_input_to_fasta(i, fastaname, input_csv_as_df)
    run_vsearch(fastaname)
    parse_vsearch_allpairsglobal_blast6(vsearch_output, curated_vsearch_output)
    add_similarity_columns(curated_vsearch_output)
    create_output_file(output_file, curated_vsearch_output)
    name_sequences(output_file)
    join_dataframes(output_file)
    tidy_files(curated_vsearch_output, fastaname, vsearch_output)
    # summary_statement(output_directory, start_time)

def import_csv(input_sequences_list, output_file):
    input_csv_as_df = pd.read_csv(input_sequences_list)
    length_of_csv = len(input_csv_as_df.index.tolist())
    length_as_list = list(range(length_of_csv))
    input_list = [x + 1 for x in length_as_list]
    input_csv_as_df.insert(0, "DB#", input_list)
    input_csv_as_df.to_csv(output_file, sep=',', index=False)
    base_name = os.path.splitext(input_sequences_list)[0]
    fastaname = base_name + '.fa'
    return fastaname, input_csv_as_df, length_of_csv

def convert_input_to_fasta(i, fastaname, input_csv_as_df):
    firstline = '>Q' + str(input_csv_as_df['DB#'].values[i]) + " " + str(input_csv_as_df['DB Accession Number'].values[i]) + " " + str(input_csv_as_df['DB Name'].values[i])
    secondline = str(input_csv_as_df['16S rRNA sequence'].values[i])
    with open(fastaname, 'a') as f:
        f.write(f"{firstline}\n{secondline}\n")

def run_vsearch(fastaname):
    cmd = 'vsearch --allpairs_global ' + fastaname + ' --blast6out vsearch_output.csv --acceptall'
    os.system(cmd)
    # print('A pre-prepared vsearch output has been used here')       # temp file

def xlookup(lookup_value, lookup_array, return_array, if_not_found: str = ''):
    match_value = return_array.loc[lookup_array == lookup_value]
    if match_value.empty:
        return f'"{lookup_value}" not found!' if if_not_found == '' else if_not_found

    else:
        return match_value.tolist()[0]

def parse_vsearch_allpairsglobal_blast6(vsearch_output, curated_vsearch_output):
    query_vs_query_forward_df = pd.read_csv(vsearch_output, usecols=[0, 1, 2], sep='\t',
                                            names=["Query Seq", "Target Seq", "Sim to Target"])
    query_vs_query_reverse_df = query_vs_query_forward_df[:]
    query_vs_query_reverse_df = query_vs_query_reverse_df[query_vs_query_reverse_df.columns[[1, 0, 2]]]
    query_vs_query_reverse_df = query_vs_query_reverse_df.rename(
        columns={'Target Seq': 'Query Seq', 'Query Seq': 'Target Seq'})
    df_length = len(query_vs_query_forward_df)
    vsearch_output_order_forward = pd.Series(np.arange(1, df_length + 1, 1))
    vsearch_output_order_reverse = pd.Series(np.arange(df_length + 1, df_length * 2 + 1, 1))
    query_vs_query_forward_df.insert(0, "Vsearch Order", vsearch_output_order_forward, True)
    query_vs_query_reverse_df.insert(0, "Vsearch Order", vsearch_output_order_reverse, True)
    query_vs_query_forward_df.insert(4, 'F/R', 'F', True)
    query_vs_query_reverse_df.insert(4, 'F/R', 'R', True)
    query_vs_query_df = pd.concat([query_vs_query_forward_df, query_vs_query_reverse_df], ignore_index=True, axis=0)
    query_vs_query_df['Q'] = query_vs_query_df['Query Seq'].str[1:].astype(int)
    query_vs_query_df['T'] = query_vs_query_df['Target Seq'].str[1:].astype(int)
    query_vs_query_df.loc[-1] = [0, 0, 0, 40, 0, 1, 0]
    query_vs_query_df = query_vs_query_df.sort_values(['Q', 'T'])
    query_vs_query_df.to_csv(curated_vsearch_output, sep=',', index=False)

def add_similarity_columns(curated_vsearch_output):
    curated_vsearch_output = pd.read_csv('All_vs_All.csv')

    def find_highest_similarity_for_query(input):
        selected_df = curated_vsearch_output[:]
        selected_df = selected_df.loc[selected_df['Q'] == input]
        selected_df = selected_df[selected_df['Q'] >= selected_df['T']]
        match_row = selected_df.loc[selected_df['Sim to Target'].idxmax()]
        highest_similarity = float(match_row['Sim to Target'])
        return highest_similarity

    def id_highest_similarity(input):
        selected_df = curated_vsearch_output[:]
        selected_df = selected_df.loc[selected_df['Q'] == input]
        selected_df = selected_df[selected_df['Q'] >= selected_df['T']]
        match_row = selected_df.loc[selected_df['Sim to Target'].idxmax()]
        closest_match = match_row['T']
        return closest_match

    def number_of_equal_highest_similarities(input):
        selected_df = curated_vsearch_output[:]
        selected_df = selected_df.loc[selected_df['Query Seq'] == input]
        match_row = selected_df.loc[selected_df['Sim to Target'].idxmax()]
        highest_similarity = float(match_row['Sim to Target'])
        number_of_equals = (selected_df['Sim to Target'].values == highest_similarity).sum()
        return number_of_equals

    curated_vsearch_output['Highest Sim to Target'] = curated_vsearch_output['Q'].apply(find_highest_similarity_for_query)
    curated_vsearch_output['Closest Match'] = curated_vsearch_output['Q'].apply(id_highest_similarity)
    curated_vsearch_output['Alternative Matches'] = curated_vsearch_output['Query Seq'].apply(number_of_equal_highest_similarities) - 1
    curated_vsearch_output = curated_vsearch_output.loc[:, ['Vsearch Order', 'Query Seq', 'Target Seq', 'Sim to Target', 'F/R', 'Highest Sim to Target', 'Closest Match', 'Alternative Matches', 'Q', 'T']]
    curated_vsearch_output.to_csv('All_vs_All.csv', sep=',', index=False)

def create_output_file(output_file, curated_vsearch_output):
    update_df = pd.read_csv(output_file)
    lookup_df = pd.read_csv(curated_vsearch_output)
    lookup_df = lookup_df.drop_duplicates(subset=['Q'])
    lookup_df = lookup_df.reset_index(drop=True)
    update_df['Closest Match'] = update_df['DB#'].map(lookup_df.set_index('Q')['Closest Match'])
    update_df['Similarity %'] = update_df['DB#'].map(lookup_df.set_index('Q')['Highest Sim to Target'])
    update_df['Alternative Matches'] = update_df['DB#'].map(lookup_df.set_index('Q')['Alternative Matches'])
    # update_df['Match DB#'] = update_df['Closest Match']
    update_df['Match DB#'] = update_df['Closest Match']
    update_df['Vs Name'] = update_df['Match DB#'].map(update_df.set_index('DB#')['DB Name'])
    update_df['Vs ID'] = update_df['Match DB#'].map(update_df.set_index('DB#')['DB Accession Number'])


    # update_df['Vs Name'] = update_df['Match DB#'].apply(xlookup, args=(
    #     update_df['DB#'], update_df['DB Name']))
    # update_df['Vs ID'] = update_df['Match DB#'].apply(xlookup, args=(
    #     update_df['DB#'], update_df['DB Accession Number']))
    update_df['Vs DB#'] = update_df['Match DB#'].apply(xlookup, args=(
        update_df['DB#'], update_df['DB#']))

    # update_df['Vs Name'] = update_df['Match DB#'].map(update_df.set_index('DB#')['DB Name'])
    # update_df['Vs ID'] = update_df['Match DB#'].map(update_df.set_index('DB#')['DB Accession Number'])
    # update_df['Vs DB#'] = update_df['Match DB#'].map(update_df.set_index('DB#')['DB#'])


    # update_df['Vs Name'] = update_df['Match DB#'].apply(xlookup, args=(
    #     update_df['DB#'], update_df['DB Name']))
    # update_df['Vs ID'] = update_df['Match DB#'].apply(xlookup, args=(
    #     update_df['DB#'], update_df['DB Accession Number']))
    # update_df['Vs DB#'] = update_df['Match DB#'].apply(xlookup, args=(
    #     update_df['DB#'], update_df['DB#']))
    # update_df = update_df.iloc[:, [0, 1, 2, 5, 4, 13, 12, 9, 10, 3, 7]]
    update_df = update_df.iloc[:, [0, 1, 2, 5, 4, 14, 12, 9, 10, 13, 3, 7]]
    update_df.to_csv(output_file, sep=',', index=False)

def name_sequences(output_file):
    input_df = pd.read_csv(output_file)
    d1, d2, d3 = {'Entry': [0], 'Db Species Name': ['Template'], 'A': [-1]}, dict.fromkeys('BCDEFGHIJKLMN', [0]), {
        'Vs': [0], 'Similarity': [0], 'DB': 'Silva', 'Accession No.': [0], '16S rRNA sequence': ['AAAAAGGGG']}
    naming_df = pd.DataFrame(dict(d1, **d2, **d3))
    naming_df.to_csv("Naming_dataframe.csv", sep=',', index=False)
    length_of_run = len(input_df.index.tolist())
    for i in range(length_of_run):
        name_each_sequence(input_df, i)

def name_each_sequence(input_df, i):
    naming_df = pd.read_csv("Naming_dataframe.csv")
    comparison_row = i
    max_entry = max(naming_df['Entry'].values)
    new_entry = max_entry + 1
    species_name = input_df['DB Name'].values[comparison_row]
    db_found = input_df['DB Found'].values[comparison_row]
    input_db_accession_number = input_df['DB Accession Number'].values[comparison_row]
    # date_accessed = input_df['Date Accessed'].values[comparison_row]
    sequence = input_df['16S rRNA sequence'].values[comparison_row]
    input_id = int(input_df['DB#'].values[comparison_row])
    matched_similarity_df = pd.read_csv("All_vs_All.csv")
    matched_similarity_df = matched_similarity_df.loc[matched_similarity_df['Q'] == input_id]
    matched_similarity_single_row_only = matched_similarity_df.values[0]
    closest_current_sequence = matched_similarity_single_row_only[6]
    similarity_to_closest_sequence = matched_similarity_single_row_only[5]


    Ai, Bi, Ci, Di, Ei, Fi, Gi, Hi, Ii, Ji, Ki, Li, Mi, Ni = int(naming_df['A'].values[closest_current_sequence]), int(
        naming_df['B'].values[closest_current_sequence]), int(naming_df['C'].values[closest_current_sequence]), int(
        naming_df['D'].values[closest_current_sequence]), int(naming_df['E'].values[closest_current_sequence]), int(
        naming_df['F'].values[closest_current_sequence]), int(naming_df['G'].values[closest_current_sequence]), int(
        naming_df['H'].values[closest_current_sequence]), int(naming_df['I'].values[closest_current_sequence]), int(
        naming_df['J'].values[closest_current_sequence]), int(naming_df['K'].values[closest_current_sequence]), int(
        naming_df['L'].values[closest_current_sequence]), int(naming_df['M'].values[closest_current_sequence]), int(
        naming_df['N'].values[closest_current_sequence])
    A, B, C, D, E, F, G, H, I, J, K, L, M, N = 60, 70, 80, 85, 90, 95, 98, 99, 99.5, 99.6, 99.7, 99.8, 99.9, 100

    if similarity_to_closest_sequence < A:
        Ao = int(max(naming_df['A'].values)) + 1
        Bo = Co = Do = Eo = Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = 0
    else:
        Ao = Ai
        if similarity_to_closest_sequence < B:
            df_b = naming_df[naming_df['A'] == Ai]
            Bo = int(max(df_b['B'].values)) + 1
            Co = Do = Eo = Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = 0
        else:
            Bo = Bi
            if similarity_to_closest_sequence < C:
                df_c = naming_df[(naming_df['A'] == Ai) & (naming_df['B'] == Bi)]
                Co = int(max(df_c['C'].values)) + 1
                Do = Eo = Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = 0
            else:
                Co = Ci
                if similarity_to_closest_sequence < D:
                    df_d = naming_df[(naming_df['A'] == Ai) & (naming_df['B'] == Bi) & (naming_df['C'] == Ci)]
                    Do = int(max(df_d['D'].values)) + 1
                    Eo = Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = 0
                else:
                    Do = Di
                    if similarity_to_closest_sequence < E:
                        df_e = naming_df[(naming_df['A'] == Ai) & (naming_df['B'] == Bi) & (naming_df['C'] == Ci) & (
                                    naming_df['D'] == Di)]
                        Eo = int(max(df_e['E'].values)) + 1
                        Fo = Go = Ho = Io = Jo = Ko = Lo = Mo = No = 0
                    else:
                        Eo = Ei
                        if similarity_to_closest_sequence < F:
                            df_f = naming_df[(naming_df['A'] == Ai) & (naming_df['B'] == Bi) & (naming_df['C'] == Ci) &
                                        (naming_df['D'] == Di) & (naming_df['E'] == Ei)]
                            Fo = int(max(df_f['F'].values)) + 1
                            Go = Ho = Io = Jo = Ko = Lo = Mo = No = 0
                        else:
                            Fo = Fi
                            if similarity_to_closest_sequence < G:
                                df_g = naming_df[
                                    (naming_df['A'] == Ai) & (naming_df['B'] == Bi) & (naming_df['C'] == Ci) & (
                                                naming_df['D'] == Di) & (naming_df['E'] == Ei) & (naming_df['F'] == Fi)]
                                Go = int(max(df_g['G'].values)) + 1
                                Ho = Io = Jo = Ko = Lo = Mo = No = 0
                            else:
                                Go = Gi
                                if similarity_to_closest_sequence < H:
                                    df_h = naming_df[
                                        (naming_df['A'] == Ai) & (naming_df['B'] == Bi) & (naming_df['C'] == Ci) & (
                                                    naming_df['D'] == Di) & (naming_df['E'] == Ei) & (
                                                    naming_df['F'] == Fi) & (naming_df['G'] == Gi)]
                                    Ho = int(max(df_h['H'].values)) + 1
                                    Io = Jo = Ko = Lo = Mo = No = 0
                                else:
                                    Ho = Hi
                                    if similarity_to_closest_sequence < I:
                                        df_i = naming_df[
                                            (naming_df['A'] == Ai) & (naming_df['B'] == Bi) & (naming_df['C'] == Ci) & (
                                                        naming_df['D'] == Di) & (naming_df['E'] == Ei) & (
                                                        naming_df['F'] == Fi) & (naming_df['G'] == Gi) & (
                                                        naming_df['H'] == Hi)]
                                        Io = int(max(df_i['I'].values)) + 1
                                        Jo = Ko = Lo = Mo = No = 0
                                    else:
                                        Io = Ii
                                        if similarity_to_closest_sequence < J:
                                            df_j = naming_df[
                                                (naming_df['A'] == Ai) & (naming_df['B'] == Bi) & (naming_df['C'] == Ci) & (
                                                            naming_df['D'] == Di) & (naming_df['E'] == Ei) & (
                                                            naming_df['F'] == Fi) & (naming_df['G'] == Gi) & (
                                                            naming_df['H'] == Hi) & (naming_df['I'] == Ii)]
                                            Jo = int(max(df_j['J'].values)) + 1
                                            Ko = Lo = Mo = No = 0
                                        else:
                                            Jo = Ji
                                            if similarity_to_closest_sequence < K:
                                                df_k = naming_df[(naming_df['A'] == Ai) & (naming_df['B'] == Bi) & (
                                                            naming_df['C'] == Ci) & (naming_df['D'] == Di) & (
                                                                             naming_df['E'] == Ei) & (
                                                                             naming_df['F'] == Fi) & (
                                                                             naming_df['G'] == Gi) & (
                                                                             naming_df['H'] == Hi) & (
                                                                             naming_df['I'] == Ii) & (naming_df['J'] == Ji)]
                                                Ko = int(max(df_k['K'].values)) + 1
                                                Lo = Mo = No = 0
                                            else:
                                                Ko = Ki
                                                if similarity_to_closest_sequence < L:
                                                    df_l = naming_df[(naming_df['A'] == Ai) & (naming_df['B'] == Bi) & (
                                                                naming_df['C'] == Ci) & (naming_df['D'] == Di) & (
                                                                                 naming_df['E'] == Ei) & (
                                                                                 naming_df['F'] == Fi) & (
                                                                                 naming_df['G'] == Gi) & (
                                                                                 naming_df['H'] == Hi) & (
                                                                                 naming_df['I'] == Ii) & (
                                                                                 naming_df['J'] == Ji) & (
                                                                                 naming_df['K'] == Ki)]
                                                    Lo = int(max(df_l['L'].values)) + 1
                                                    Mo = No = 0
                                                else:
                                                    Lo = Li
                                                    if similarity_to_closest_sequence < M:
                                                        df_m = naming_df[(naming_df['A'] == Ai) & (naming_df['B'] == Bi) & (
                                                                    naming_df['C'] == Ci) & (naming_df['D'] == Di) & (
                                                                                     naming_df['E'] == Ei) & (
                                                                                     naming_df['F'] == Fi) & (
                                                                                     naming_df['G'] == Gi) & (
                                                                                     naming_df['H'] == Hi) & (
                                                                                     naming_df['I'] == Ii) & (
                                                                                     naming_df['J'] == Ji) & (
                                                                                     naming_df['K'] == Ki) & (
                                                                                     naming_df['L'] == Li)]
                                                        Mo = int(max(df_m['M'].values)) + 1
                                                        No = 0
                                                    else:
                                                        Mo = Mi
                                                        if similarity_to_closest_sequence < N:
                                                            df_n = naming_df[
                                                                (naming_df['A'] == Ai) & (naming_df['B'] == Bi) & (
                                                                            naming_df['C'] == Ci) & (
                                                                            naming_df['D'] == Di) & (
                                                                            naming_df['E'] == Ei) & (
                                                                            naming_df['F'] == Fi) & (
                                                                            naming_df['G'] == Gi) & (
                                                                            naming_df['H'] == Hi) & (
                                                                            naming_df['I'] == Ii) & (
                                                                            naming_df['J'] == Ji) & (
                                                                            naming_df['K'] == Ki) & (
                                                                            naming_df['L'] == Li) & (naming_df['M'] == Mi)]
                                                            No = int(max(df_n['N'].values)) + 1
                                                        else:
                                                            No = Ni

    # new_line = {'Entry': [new_entry],'Db Species Name': [species_name], 'A': [Ao], 'B': [Bo], 'C': [Co], 'D': [Do], 'E': [Eo], 'F': [Fo], 'G': [Go], 'H': [Ho], 'I': [Io], 'J': [Jo], 'K': [Ko], 'L': [Lo], 'M': [Mo], 'N': [No], 'Vs': [closest_current_sequence], 'Similarity': [similarity_to_closest_sequence], 'DB': [db_found], 'Accession No.': [input_db_accession_number], 'Hun#': [additional_ordering_one], 'Fas#': [additional_ordering_two], '16S rRNA sequence': [sequence] }
    new_line = {'Entry': [new_entry], 'Db Species Name': [species_name], 'A': [Ao], 'B': [Bo], 'C': [Co], 'D': [Do],
                'E': [Eo], 'F': [Fo], 'G': [Go], 'H': [Ho], 'I': [Io], 'J': [Jo], 'K': [Ko], 'L': [Lo], 'M': [Mo],
                'N': [No], 'Vs': [closest_current_sequence], 'Similarity': [similarity_to_closest_sequence],
                'DB': [db_found], 'Accession No.': [input_db_accession_number], '16S rRNA sequence': [sequence]}
    new_line_df = pd.DataFrame(new_line)
    naming_df = pd.concat([naming_df, new_line_df], ignore_index=True, axis=0)
    naming_df.to_csv("Naming_dataframe.csv", sep=',', index=False)

def join_dataframes(output_file):
    input_df = pd.read_csv(output_file)
    name_df = pd.read_csv('Naming_dataframe.csv')
    name_df = name_df.drop(0).reset_index(drop=True)

    columns_to_be_inserted = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N']
    insert = name_df.loc[:, columns_to_be_inserted]
    output_df = pd.concat([input_df, insert], 1)
    # output_df = output_df.iloc[:, [0, 1, 2, 4, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 5, 6, 7, 8, 3, 9, 10]]
    output_df = output_df.iloc[:, [0, 1, 2, 4, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 5, 6, 7, 8, 9, 3, 10, 11]]

    output_df.loc[0, 'Vs ID'] = 'N/A'
    output_df.loc[0, 'Vs Name'] = 'N/A'
    output_df.loc[0, 'Similarity %'] = 'N/A'
    output_df.loc[0, 'Vs DB#'] = 'N/A'
    output_df.loc[0, 'Alternative Matches'] = 0

    output_df = output_df.sort_values(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N'])
    length_of_output_df = len(output_df.index.tolist())
    length_as_list = list(range(length_of_output_df))
    input_list = [x + 1 for x in length_as_list]
    output_df.insert(3, "Tree#", input_list)
    output_df.to_csv(output_file, sep=',', index=False)

def tidy_files(curated_vsearch_output, fastaname, vsearch_output):
    os.mkdir('data')
    shutil.move(curated_vsearch_output, 'data')
    shutil.move(vsearch_output, 'data')
    shutil.move('Naming_dataframe.csv', 'data')
    shutil.move(fastaname, 'data')

def summary_statement(output_directory, start_time):
    program_length = datetime.now() - start_time
    print(f'{sys.argv[0]} has completed in {program_length}. {sys.argv[1]} was processed and the Summary.csv file can be '
          f'found in the new {output_directory} directory.')

main()
print('finished')