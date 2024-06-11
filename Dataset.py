import pandas as pd

def read_fasta(file):
    '''simple function to read FASTA file and identify the sequence using pandas'''
    sequence = ''
    for line in file:
        if line[0] == '>':
            id = line[1:-1]
        else:
            sequence += line.strip()
    sequence_df = pd.Series(sequence)
    return sequence_df

pd.set_option("display.max_colwidth", None)
dataset = read_fasta('sequence.fasta')