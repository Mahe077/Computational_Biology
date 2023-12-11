import re

from Bio import Entrez, SeqIO


def fetch_dengue_sequences():
    Entrez.email = "maheshlakshan766@gmail.com"
    handle = Entrez.esearch(db="nucleotide", term="Dengue 1 Envelope", retmax=10)
    record = Entrez.read(handle)
    ids = record["IdList"]

    dengue_sequences = []
    for record_id in ids:
        handle = Entrez.efetch(db="nucleotide", id=record_id, rettype="gb", retmode="text")
        seq_record = SeqIO.read(handle, "genbank")
        dengue_sequences.append(seq_record)

    return dengue_sequences


def write_fasta_file(sequences, output_file="dengue_sequences.fasta"):
    with open(output_file, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")


def clean_up_header(header):
    pattern = r'(\S+).*?\((\w+)\) gene'

    match = re.search(pattern, header)

    if match:
        accession_number = match.group(1)
        gene_name = match.group(2)
        cleaned_header = f">{accession_number}_{gene_name}"
        return cleaned_header
    else:
        return ">Unknown"


def print_multi_fasta(file_path):
    sequences = list(SeqIO.parse(file_path, "fasta"))

    for seq_record in sequences:
        cleaned_header = clean_up_header(seq_record.description)
        print(cleaned_header)
        print(seq_record.seq)


if __name__ == "__main__":
    dengue_sequences = fetch_dengue_sequences()
    write_fasta_file(dengue_sequences)
    print_multi_fasta("dengue_sequences.fasta")
