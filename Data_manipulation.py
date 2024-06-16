import pandas as pd
import matplotlib.pyplot as plt
from Dataset import *
from abc import ABC, abstractmethod

class Sequence(ABC):
    def __init__(self, sequence):
        '''initializes every sequence that can be obtained from the FASTA file
        '''
        self.__sequence = sequence

    @abstractmethod
    def get_frequencies(self):
        '''returns the occurrence count of each base in the sequence
        '''
        pass
    
    @abstractmethod
    def transcription(self):
        '''return the RNA sequence obtained from DNA by replacing Ts with Us,
        if input is already an RNA file, just returns itself
        '''
        pass
    
    @abstractmethod
    def get_sequence(self):
        '''returns the visual representation of the sequence
        '''
        pass
        
class DNA(Sequence):
    def __init__(self, sequence):
        '''initializes a DNA object
        '''
        dna_sequence_list = []
        for i in sequence[0]:
            dna_sequence_list.append(Nucleotide(i))
        self.__sequence = pd.Series(dna_sequence_list)
        
    def transcription(self):
        rna_sequence_list = []
        for i in self.__sequence:
            if i.get_base() == "T":
                rna_sequence_list.append(Nucleotide("U"))
            else:
                rna_sequence_list.append(i)
        return mRNA(rna_sequence_list)

    def get_frequencies(self):
        base_list = [nucl.get_base() for nucl in self.__sequence]
        frequencies = {}
        for i in ["A", "T", "C", "G"]:
            frequencies[i] = base_list.count(i)
        return frequencies
    
    def get_sequence(self):
        base_sequence = [nucl.get_base() for nucl in self.__sequence]
        string_seq = ''.join(base_sequence)
        return string_seq

class mRNA(Sequence):
    def __init__(self, sequence):
        '''initializes an mRNA object
        '''
        self.__sequence = pd.Series(sequence)
        
    def transcription(self):
        return self.__sequence
        
    def get_frequencies(self):
        base_list = [nucl.get_base() for nucl in self.__sequence]
        frequencies = {}
        for i in ["A", "U", "C", "G"]:
            frequencies[i] = base_list.count(i)
        return frequencies

    def get_sequence(self):
        base_sequence = [nucl.get_base() for nucl in self.__sequence]
        string_seq = ''.join(base_sequence)
        return string_seq

    def reverse_complementary(self):
        '''returns the reverse complementary of the RNA sequence
        '''
        string_seq = self.get_sequence()
        reverse = string_seq[::-1]
        rev_comp = ''
        changes = {"C": "G", "G": "C", "A": "U", "U": "A"}
        for i in reverse:
            rev_comp += changes[i]
        return rev_comp

    def translation(self):
        '''returns the amino acid sequence obtained from the RNA considering all the possible ORFs
        and, depending on the lenght of each product, stores each as protein or oligopeptide by calling
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
            for i in range(frame, len(self.__sequence) -2, 3):
                codon = self.__sequence[i:i + 3]
                str_codon = ''.join([nucl.get_base() for nucl in codon])
                str_aa = codons.get(str_codon, '')
                aa_sequence += str_aa
            orf_sequences.append(aa_sequence)
        
        minus = self.reverse_complementary()
        
        for frame in range(3): # Generates three ORFs from the minus (-) strand # COPIARE PLUS STRAND
            aa_sequence = ''
            for i in range(frame, len(minus) -2, 3):
                str_codon = minus[i:i+3]
                str_aa = codons.get(str_codon, '')
                aa_sequence += str_aa
            orf_sequences.append(aa_sequence)
        
        protein_data = {}
        for p in range(6):
            data = AminoAcidChain(orf_sequences[p])
            protein_data[p+1] = data.sequence_type()
            
        return orf_sequences, protein_data
    
class AminoAcidChain(Sequence):
    def __init__(self, sequence):
        '''initializes an amino acid chain object
        '''
        l  = []
        for i in sequence:
            l.append(AminoAcid(i))
        self.sequence = pd.Series(l)

    def transcription(self):
        return self
    
    def get_frequencies(self):
        aa_list = [aa.get_aa() for aa in self.sequence]
        frequencies = {}
        for c in ["A", "C", "D", "E","F", "G", "H", "I","J","L", "M", "N", "P","Q", "R", "S", "T","V", "W", "Y"]:
            frequencies[c] = aa_list.count(c)
        return frequencies

    def get_sequence(self):
        aa_sequence = [aa.get_aa() for aa in self.sequence]
        string_seq = ''.join(aa_sequence)
        return string_seq
        
    def sequence_type(self): 
        '''returns a dictionary storing amino acid sequences as proteins or oligopeptides
        depending on their lenght (each product starts with M and end with a stop codon (*)
        '''
        sequence = self.sequence
        chain = {"Protein":[], "Oligo":[]}
        idx = 0
        while idx < len(sequence):
            seq = ""
            if sequence[idx].get_aa() == "M":
                while sequence[idx].get_aa() != "*":
                    seq += sequence[idx].get_aa()
                    idx += 1
                    if idx == len(sequence):
                        break
                if len(seq) > 20:
                    seq = Protein(seq)
                    chain["Protein"] += [(len(seq.get_sequence()), seq.get_sequence())]
                elif len(seq) > 0:
                    seq = Oligopeptide(seq)
                    chain["Oligo"] += [(len(seq.get_sequence()), seq.get_sequence())]
            else:
                idx += 1
        for k,v in chain.items():
            chain[k] = sorted(chain[k], reverse=True)
        return chain

class Protein(AminoAcidChain):
    '''identifies a protein'''
    def __str__(self):
        return self.sequence
    
    def get_sequence(self):
        aa_sequence = [aa.get_aa() for aa in self.sequence]
        string_seq = ''.join(aa_sequence)
        return string_seq

class Oligopeptide(AminoAcidChain):
    '''identifies an Oligopeptide'''
    def __str__(self):
        return self.sequence
    
    def get_sequence(self):
        aa_sequence = [aa.get_aa() for aa in self.sequence]
        string_seq = ''.join(aa_sequence)
        return string_seq

class AminoAcid():
    def __init__(self, aa):
        '''initializes an amino acid object
        '''
        self.__aa = aa
    def get_aa(self):
        return self.__aa

class Nucleotide():
    def __init__(self, base):
        '''initializes a nucleotide object
        '''
        self.__base = base
    def get_base(self):
        return self.__base
