from Bio.Data import CodonTable

standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
#standard_table = CodonTable.unambiguous_dna_by_id[1]

# Vertebrate Mitochondrial: ty thể của động vật có xương sống
# mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
# mito_table = CodonTable.unambiguous_dna_by_id[2]

print(standard_table)

#print(mito_table)
print("Stop codons: ", standard_table.stop_codons)
print("Start codons: ", standard_table.start_codons)
print(standard_table)

from Bio import Align
aligner = Align.PairwiseAligner()
aligner = Align.PairwiseAligner(match_score=1.0)
