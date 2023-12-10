from Bio.Seq import Seq
from Bio.SeqUtils import nt_search

dna_seq = Seq(input("Enter the DNA sequence: "))

dna_seq_str = str(dna_seq)


def find_cpg_islands(sequence):
    cpg_positions = nt_search(sequence, "CG")[1:]
    islands = []
    current_island = []

    for pos in cpg_positions:
        if not current_island or pos == current_island[-1] + 1:
            current_island.append(pos)
        else:
            islands.append(current_island)
            current_island = [pos]

    islands.append(current_island)
    return islands


cpg_islands = find_cpg_islands(dna_seq_str)

print("CpG Islands:", cpg_islands)
