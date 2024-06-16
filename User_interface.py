from Data_manipulation import *
from flask import Flask, render_template, request, redirect, url_for
import matplotlib.pyplot as plt

app = Flask(__name__)

global dna_sequence, rna_seq
dna_sequence = None

#----------------------------------------------------------------------------------

def generate_base_frequencies_chart(frequencies):
    bases = list(frequencies.keys())
    counts = list(frequencies.values())
    colors = ['red', 'orange', 'blue', 'green'] 
    plt.bar(bases, counts, color=colors)
    plt.xlabel('Base')
    plt.ylabel('Frequency')
    plt.title('Base Frequencies')
    chart_path = 'static/base_frequencies.png'
    plt.savefig(chart_path)
    plt.close()

#----------------------------------------------------------------------------------

@app.route('/')
def landing():
    return render_template('landing.html')

#----------------------------------------------------------------------------------

@app.route('/upload', methods=['POST'])
def upload_file():
    global dna_sequence, rna_seq
    fasta_file = request.files['fasta-file']
    if fasta_file.filename == '':
        return redirect(request.url)
    if fasta_file and fasta_file.filename.endswith('.fasta'):
        file_path = fasta_file.filename
        fasta_file.save(file_path)
        with open(file_path, 'r') as file:
            dataset = read_fasta(file)
        dna_sequence = DNA(dataset)
        rna_seq = dna_sequence.transcription()
        return redirect(url_for('home'))
    else:
        return "Invalid file format. Please upload a FASTA file."

#----------------------------------------------------------------------------------

@app.route('/home')
def home():
    if dna_sequence is None:
        return redirect(url_for('landing'))
    stats = dna_sequence.get_frequencies()
    generate_base_frequencies_chart(stats)
    Bases = {"A":"Adenine","T":"Thymine","G":"Guanine","C":"Cytosine"}
    maxi = max((v, k) for k, v in stats.items())
    mini = min((v, k) for k, v in stats.items())
    names = (Bases[maxi[1]],Bases[mini[1]])  
    dna_sequence_string = dna_sequence.get_sequence()
    return render_template('home.html',
                           dna_seq=dna_sequence_string,
                           stats=[stats, maxi, mini,names])

#----------------------------------------------------------------------------------

@app.route('/transcription_translation', methods=['GET', 'POST'])
def transcribe_and_translate():
    if dna_sequence is None:
        return redirect(url_for('landing'))
    stats = rna_seq.get_frequencies()
    generate_base_frequencies_chart(stats)
    Bases = {"A": "Adenine", "U": "Uracil", "G": "Guanine", "C": "Cytosine"}
    maxi = max((v, k) for k, v in stats.items())
    mini = min((v, k) for k, v in stats.items())
    names = (Bases[maxi[1]], Bases[mini[1]])
    sort_order = request.args.get('sort', 'asc')
    protein_sequences, orf_dict = rna_seq.translation()
    
    for key in orf_dict:
        proteins = orf_dict[key]['Protein']
        oligos = orf_dict[key]['Oligo']
        if sort_order == 'asc':
            orf_dict[key]['Protein'] = sorted(proteins, key=lambda x: x[0])
            orf_dict[key]['Oligo'] = sorted(oligos, key=lambda x: x[0])
        else:
            orf_dict[key]['Protein'] = sorted(proteins, key=lambda x: x[0], reverse=True)
            orf_dict[key]['Oligo'] = sorted(oligos, key=lambda x: x[0], reverse=True)
    
    return render_template('transcription_translation.html',
                           stats=[stats, maxi, mini, names],
                           rna_sequence=rna_seq.get_sequence(),
                           pseq=protein_sequences,
                           PROTEINS=orf_dict,
                           sort_order=sort_order)

#----------------------------------------------------------------------------------

if __name__ == '__main__':
    app.run(debug=True)
