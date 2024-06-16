DNA Sequence Analyzer

This project is a Python tool with an HTML GUI that can help when dealing with genetic analysis. Starting from a FASTA file of a DNA sequence, it can replicate biological processes such as transcription into mRNA and translation into an amino acid sequence, as well as provide statistical information on the sequence content. The user interface was realized using the Flask module. The .png images of the frequency histograms are generated with the Matplotlib module. The sequences are stored as arrays through the Pandas module.


How to use

When you first open the program you will find yourself on the landing page. From here, by clicking on the "Choose file" button, you will be able to upload a FASTA file (.fasta) to perform some genetic analysis on its content. Clicking the "Go" button will start the analysis, taking you to the home page. Within the project files, two FASTA files ("sequence.fasta" and "sequence_short.fasta") have been provided to help you grasp the functions of the program.

On the home page, you will be able to visualize the entire DNA sequence in a text box. Below the text box, you will see the histogram describing the absolute frequency of the 4 nucleotides in your sequence. You will find a picture of the last histogram created stored as a .png file in the static project folder. Additionally, you will have the total count for each amino acid and will know the most and least frequent ones. Below the statistics, you will read the line "Click here to convert sequence to mRNA and to transcribe it," which will take you to the transcription_translation page.

On the transcription_translation page, you will be able to visualize the entire mRNA sequence in a text box. Below the text box, you will see the same graph and statistical information you got for the DNA sequence. In the Translation section, you will be able to visualize the entire amino acid chains obtained from the translation of all the possible open reading frames divided per strand. Then, for each reading frame, you will get two text boxes containing all the proteins and all the oligopeptides that have been synthesized, and you will know the length of each one. The chains can be sorted in ascending or descending order based on their length.


File content

Dataset.py
	The file is responsible for the handling of the FASTA file. Using the Pandas module, read_fasta() reads, stores, and returns the nucleotide content of a file in a one-dimensional array sequence_df.

Data_manipulation.py
	The file contains all the classes for genetic entities and their methods used to replicate biological processes or employed during statistical analysis.
	Sequence is an abstract class used to represent any genetic object defined by a unique sequence.
	DNA is a Sequence subclass characterized by a sequence of Nucleotide objects stored in a private Pandas one-dimensional array.
		The transcription() method returns an mRNA object made of a series of Nucleotide objects by changing every T (Thymine) element with a U (Uracil) element.
		The get_frequencies() method returns a dictionary having as keys each string representation of the nucleotides present in the sequence and as values their absolute frequencies.
		The get_sequence() method returns a string representation of the sequence of Nucleotide objects of its DNA object.
	mRNA is a Sequence subclass characterized by a sequence of Nucleotide objects stored in a private Pandas one-dimensional array.
		The transcription() method returns its own sequence attribute since in an mRNA the transcription has already taken place and it wouldn't cause any change in the nucleotide content.
		The get_frequencies() method is identical to the DNA one except for the keys of the frequencies dictionary, which reflect that Thymines have been substituted by Uracils.
		The get_sequence() method is identical to the DNA one and returns a string representation of the sequence of Nucleotide objects of its mRNA object.
		The reverse_complementary() method creates the hypothetical complementary sequence to the one represented by the mRNA object. It first stores as a list the reversed sequence so as to respect the 5' to 3' reading direction, then uses a dictionary to map each base character to its complementary base and returns the complementary sequence as a string.
		The translation() method returns a list orf_sequences of length 6 with the string representations of the amino acid chains translated in each reading frame (3 on the "+ strand" and 3 on the "- strand") and a dictionary protein_data having as keys indexes representing one of the 6 reading frames and as values the dictionaries obtained from the AminoAcidChain.sequence_type() method containing all the proteins and oligopeptides synthesized during translation. The method contains a codons dictionary used for translation based on the human codon table. It can be modified from the source code if working with species that have a different codon table. The method stores in a list the string representations of the amino acid chains. In order to also get the string representations of the amino acid chains of the minus strand, the reverse_complementary() method is called in order to obtain the complementary sequence. The string representations are then used to create AminoAcidChain objects from which the sequence_type() method is called in order to obtain the list of proteins and oligopeptides of each reading frame translation.
	AminoAcidChain is a Sequence subclass characterized by a sequence of AminoAcid objects stored in a Pandas one-dimensional array.
		The transcription() method returns its own sequence attribute since transcription can't be applied to an amino acid sequence.
		The get_frequencies() method is identical to the DNA one except for the keys of the frequencies dictionary, which don't represent the 4 nucleotides but the 20 amino acids instead.
		The get_sequence() method is identical to the DNA one and returns a string representation of the sequence of AminoAcid objects of its AminoAcidChain object.
		The sequence_type() method returns a dictionary chain made of 2 keys. The "Protein" key has as value all the Protein objects obtained from translation and their length, the "Oligo" key has as value all the Oligopeptide objects obtained from translation and their length. When reading the translation sequence, once an "M" (Methionine) character is found, the sequence is stored as a string character along with its length, until a "*" (stop codon) character is found. The sequence is then used to initialize a Protein or Oligopeptide object depending on its length, which is then stored as a string in the dictionary. When the entire sequence has been analyzed, the values are ordered in ascending order based on the length of each element.
	Protein is an AminoAcidChain subclass characterized by a sequence of AminoAcid objects stored in a Pandas one-dimensional array. The sequence has length > 20.
		The get_sequence() method returns a string representation of the sequence of AminoAcid objects of its Protein object.
	Oligopeptide is an AminoAcidChain subclass characterized by a sequence of AminoAcid objects stored in a Pandas one-dimensional array. The sequence has length <= 20.
		The get_sequence() method returns a string representation of the sequence of AminoAcid objects of its AminoAcid object.
	AminoAcid is a class characterized by a private string of length 1 used to represent the corresponding amino acid.
		The get_aa() method returns a string representation of the amino acid sequence.
	Nucleotide is a class characterized by a private string of length 1 used to represent the corresponding nucleotide.
		The get_base() method returns a string representation of the nucleotide.

User_interface.py
	The generate_base_frequencies_chart() method takes the frequency dictionary created by some Sequence.get_frequencies() method and, using Matplotlib, plots a histogram with the bases as its x-axis (the keys of frequency) and their absolute frequency as the y-axis (the values of frequency). The histogram is then saved as "base_frequencies.png" in the static folder.
	The upload_file() method ensures the provided file is in .fasta format and, if so, opens and stores its content in dataset by calling the read_fasta() method. dataset is then used to initialize a DNA object dna_sequence, from which an mRNA object rna_seq is also initialized by calling the DNA.transcription() method in its declaration.
	The home() method ensures the DNA sequence is not empty and, if so, retrieves the dictionary with the frequencies through the DNA.get_frequencies() method, which is then used to call generate_base_frequencies_chart() and visualize the graph. Then it finds the key (the nucleotide) of the frequencies with the highest and lowest value (their count).
	The transcribe_and_translate() method ensures the DNA sequence is not empty and, if so, like the home() method plots the histogram and provides the nucleotide with the highest and lowest frequency. Then it provides the amino acid sequence as a string protein_sequences and the Protein and Oligopeptide content in a dictionary orf_dict by calling the mRNA.translation() method. The orf_dict is then sorted in ascending or descending order based on the user preference.

The templates folder contains the 3 .html files "home.html", "landing.html" and "transcription_translation.html" which make up the graphical user interface of the program.


We hope our tool will be helpful and speed up your future research in the Genomics field.


Contributors

	Almatari Seif (https://github.com/SeifAlmatari)
	Nicolis Tommaso (https://github.com/TommasoNicolis)
	Piha Giacomo (https://github.com/GiacomoPiha)
	Todaro Adriano (https://github.com/AdrianoTodaro)
