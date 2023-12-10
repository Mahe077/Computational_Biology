from Bio.Seq import Seq


dna_sequence = input("Enter the DNA sequence: ")

dna_seq = Seq(dna_sequence)

protein_seq = dna_seq.translate()

print("Translated Protein Sequence:", protein_seq)
