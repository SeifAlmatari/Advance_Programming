This folder contains the following:

3 python programs:
1) Data_manipulation 	- contains all the classes and methods
2) Dataset 		- for file handling, identifying the sequence in a FASTA file
3) User_interface	- the main body which calls all the other files and templates
-----------------------------------------------------------------------------------------

templates folder contains 3 html files
1) landing - this is the first page the user sees when opening the website, allows the 
	     user to input FASTA file
2) home	   - once the file is processed, this page presents the DNA sequence along with 
	     its basic statistics, a graph, and a mode to be redirected to the next page 
             to transcribe the DNA into RNA and translate it into its amino acid sequence
3) transcription_translation - presents the RNA sequence, all 6 ORFs, and the protein and 
			       oligopeptide sequences for each ORF
-----------------------------------------------------------------------------------------

This folder also contains:
> static folder, this automatically saves a (temp) copy of the DNA statistics graph once 
  called. required for ease of file handling and to avoid recurring errors
> a folder called __pycache__ that seems to save automatically once you run the code
> 2 FASTA sequences, the one provided in the project specification, and the other being
  a shorter version for quick testing during debugging
-----------------------------------------------------------------------------------------

Instructions:
1) The FASTA file containing the sequence must be in the same directory as the python files
2) make sure that you have python and all the modules already installed (especially if you
   are using the terminal)

