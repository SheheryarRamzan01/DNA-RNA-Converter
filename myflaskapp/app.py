from flask import Flask, render_template, request

app = Flask(__name__)

codon_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def translate_to_protein(sequence, is_rna):
    if is_rna:
        sequence = sequence.replace('U', 'T')  # Convert U to T in RNA sequences
    protein_sequence = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if codon in codon_table:
            amino_acid = codon_table[codon]
            if amino_acid == 'Stop':
                break 
            protein_sequence.append(amino_acid)
        else:
            protein_sequence.append('X') 
    return protein_sequence



@app.route('/', methods=['GET', 'POST'])
def translate_sequence():
    if request.method == 'POST':
        sequence = request.form['sequence'].upper()
        if all(base in "ACGTU" for base in sequence):
            is_rna = 'U' in sequence
            protein_sequence = translate_to_protein(sequence, is_rna)
            sequence_type = "RNA" if is_rna else "DNA"
            return render_template('result.html', sequence=sequence, protein_sequence=protein_sequence, sequence_type=sequence_type)
        else:
            error = "Invalid input. Please enter a valid DNA or RNA sequence."
            return render_template('index.html', error=error)
    return render_template('index.html')

if __name__ == '__main__':
    app.run(debug=True)