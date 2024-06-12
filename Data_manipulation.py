from Dataset import *
import matplotlib.pyplot as plt
from abc import ABC, abstractmethod

class Sequence(ABC):
    def __init__(self, sequence):
        '''initialization for each sequence from FASTA file
        '''
        self.sequence = sequence

    @abstractmethod
    def get_frequencies(self):
        pass

    @abstractmethod
    def transcription(self):
        '''converts DNA to RNA sequence by replacing Ts with Us, if its already an
        RNA file, will just return itself
        '''
        pass
        
class DNA(Sequence):
    '''counts the frequencies of the bases in the DNA sequence
    '''
    def __init__(self, sequence):
        dna_sequence_list = []
        for i in sequence:
            dna_sequence_list.append(Nucleotide(i))
        self.sequence = pd.Series(dna_sequence_list)
        
    def transcription(self):
        rna_sequence_list = []
        for i in self.sequence:
            if i.base = "T"
                rna_sequence_list.append(Nucleotide("U"))
        else:
            rna_sequence.append(i)
        return rna_sequence_list

    def get_frequencies(self):
        base_list = [nucl.base for nucl in self.sequence]
        frequencies = {}
        for i in ["A", "T", "C", "G"]:
            frequencies[i] = base_list.count(i)
        return frequencies

class mRNA(Sequence):
    def __init__(self, sequence):
        self.sequence = pd.Series(sequence)
        
    def transcription(self):
        return self.sequence
        
    def get_frequencies(self):
        base_list = [nucl.base for nucl in self.sequence]
        frequencies = {}
        for i in ["A", "U", "C", "G"]:
            frequencies[i] = base_list.count(i)
        return frequencies

    def translation(self):
        '''
        converts an RNA sequence into an amino acid sequence, then calls separate
        method to identify each subsequence as either protein or oligopeptide
        ultimately returning a dictionary in this format:
        dictionary = {1:{"Protein":[],"Oligo":[]}, 2:...}
        The numbers as keys indicate the ORFs
        '''
        codons = {
            "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
            "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
            "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
            "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
            "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
            "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
            "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
            "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
            "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
            "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
            "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
            "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
            "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
            "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
            "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
            "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
        }
        orf_sequences = []
        for frame in range(3): # Generates three ORFs from the plus (+) strand
            aa_sequence = ''
            for i in range(frame, len(self.sequence) -2, 3):
                codon = self.sequence[i:i + 3]
                str_codon = str([nucl.base for nucl in codon])
                str_aa = codons.get(str_codon, '')
                aa_sequence += str_aa
            orf_sequences.append(aa_sequence)
        
        minus = self.sequence.to_string()[::-1]
        new_minus = ''
        changes = {"C":"G", "G":"C", "A":"U", "U":"A", '':'', ' ':' ','0':'0'}
        for i in range(len(minus)): # Generates complementary strand
            new_minus += changes[minus[i]] 
        
        for frame in range(3): # Generates three ORFs from the minus (-) strand # COPIARE PLUS STRAND
            aa_sequence = ''
            for i in range(frame, len(new_minus) -2, 3):
                str_codon = new_minus[i:i+3]
                str_aa = codons.get(str_codon, '')
                aa_sequence += str_aa
            orf_sequences.append(aa_sequence)
        
        '''this part calls the other method to generate all the proteins / oligo chains'''
        protein_data = {}
        for p in range(6):
            data = AminoAcidChain(orf_sequences[p])
            protein_data[p+1] = data.sequence_type()
            
        return orf_sequences, protein_data
    
class AminoAcidChain(Sequence):
    def __init__(self, sequence):
        l  = []
        for i in sequence:
            l.append(AminoAcid(i))
        self.sequence = pd.Series(l)
    
    def get_frequencies(self):
        aa_list = [aa.aa for aa in self.sequence[0]]
        frequencies = {}
        for c in ["A", "C", "D", "E","F", "G", "H", "I","J","L", "M", "N", "P","Q", "R", "S", "T","V", "W", "Y"]:
            frequencies[c] = aa_list.count(c)
        return frequencies
        
    def sequence_type(self): # SEQUENCE TYPE DA MODIFICARE
        '''
        returns a dictionary, method called during translation
        to return a dictionary with all the protein and oligopeptide
        sequences,only considers subsequences that start with M and end with *
        '''
        sequence = self.sequence
        chain = {"Protein":[], "Oligo":[]}
        idx = 0
        while idx < len(sequence):
            seq = ""
            if sequence[idx] == "M":
                while sequence[idx] != "*":
                    seq += sequence[idx]
                    idx += 1
                    if idx == len(sequence):
                        break
                if len(seq) > 20:
                    seq = str(Protein(seq))
                    chain["Protein"] += [(len(seq), seq)]
                elif len(seq) > 0:
                    seq = str(Oligopeptide(seq))
                    chain["Oligo"] += [(len(seq), seq)]
            else:
                idx += 1
        for k,v in chain.items():
            chain[k] = sorted(chain[k], reverse=True)
        return chain

class Protein(Sequence):
    '''to identify the proteins sequences as Protein'''
    def __str__(self):
        return self.sequence

class Oligopeptide(Sequence):
    '''to idenitfy the oligo sequences as Oligopeptides'''
    def __str__(self):
        return self.sequence

class AminoAcid():
    '''when converting the mRNA sequence using codons during transcription'''
    def __init__(self, aa):
        self.aa = aa

class Nucleotide():
    '''for each of the bases'''
    def __init__(self, base):
        self.base = base
