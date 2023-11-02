import Bio
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import GC

# FASTA parsing example
for gene_record in SeqIO.parse("gene.fna", "fasta"):
    print("Id của chuỗi: ", gene_record.id)
    print(repr(gene_record.seq))
    # print(len(gene_record.seq), "nucleotides")
    mRNA = (
        gene_record.seq[0:86]
        + gene_record.seq[169:286]
        + gene_record.seq[369:819]
        + gene_record.seq[905:]
    )
    print(len(mRNA))

for rna_record in SeqIO.parse("rna.fna", "fasta"):
    print("\nmRNA")
    print("Id của chuỗi: ", rna_record.id)
    print(repr(rna_record.seq))
    # print(len(rna_record.seq), "nucleotides")

# print(mRNA == rna_record.seq)

remainder = len(mRNA) % 3
new_mRNA = mRNA + Seq("N")
protein_seq = new_mRNA[32:1781].translate(to_stop=True)
# print(protein_seq)
# print(len(protein_seq))
# print(len(new_mRNA))


# print("protein_seq_len: ", len(protein_seq))
for protein_record in SeqIO.parse("protein.faa", "fasta"):
    print("\nProtein")
    print("Id của chuỗi: ", protein_record.id)
    print(repr(protein_record.seq))
    print(len(protein_record.seq), "amino acids")

print(protein_seq == protein_record.seq)
# print(mRNA)
# print("\n\n", protein_record.seq)
