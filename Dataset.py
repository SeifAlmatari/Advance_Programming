import pandas as pd

def read_fasta(file):
    '''reads a FASTA file and converts it into a Pandas series
    '''
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
