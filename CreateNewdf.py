# load packages
import pandas as pd

# create Output pandas df
data = {'Entry': [1],
        'Db Species Name': ['Template'],
        'A': [0],
        'B': [0],
        'C': [0],
        'D': [0],
        'E': [0],
        'F': [0],
        'G': [0],
        'H': [0],
        'I': [0],
        'J': [0],
        'K': [0],
        'L': [0],
        'M': [0],
        'N': [0],
        'O': [0],
        'P': [0],
        'Q': [0],
        'Vs': [0],
        'Similarity': [0],
        '16S Sequence': ['AAAAAGGGG']
        }

Output_df = pd.DataFrame(data)
Output_df.to_csv("RunningDataframe.csv", sep=',', index=False)
