# Working with Sequences using fasta files
from Bio import SeqIO
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis


covid = SeqIO.read("covid_sequences.fasta", "fasta")
mers = SeqIO.read("mers_sequence.fasta", "fasta")
sars = SeqIO.read("sars_sequence.fasta", "fasta")
ebola = SeqIO.read("ebola_sequence.fasta", "fasta")
hiv = SeqIO.read("hiv_sequence.fasta", "fasta")

covid_seq = covid.seq
mers_seq = mers.seq
sars_seq = sars.seq
ebola_seq = ebola.seq
hiv_seq = hiv.seq

# Check the length of each sequence
# print("covid_seq ::", len(covid_seq))
# print("mers_seq ::", len(mers_seq))
# print("sars_seq ::", len(sars_seq))
# print("ebola_seq ::", len(ebola_seq))

# Check the length of each sequence
print("GC content of covid_seq ::", gc_fraction(covid_seq))
print("GC content of mers_seq ::", gc_fraction(mers_seq))
print("GC content of sars_seq ::", gc_fraction(sars_seq))
print("GC content of ebola_seq ::", gc_fraction(ebola_seq))
print("GC content of hiv_seq ::", gc_fraction(hiv_seq))

def pad_seq(seq):
    if len(seq) % 3 == 0:
        return seq
    elif len(seq) % 3 == 1:
        return seq + Seq("NN")
    else:
        return seq + Seq("N")
    
covid_proteins = []
mers_proteins = []
sars_proteins = []
ebola_proteins = []
hiv_proteins = []

for i in range(3):
    covid_proteins.append(pad_seq(covid_seq[i:]).translate())
    mers_proteins.append(pad_seq(mers_seq[i:]).translate())
    sars_proteins.append(pad_seq(sars_seq[i:]).translate())
    ebola_proteins.append(pad_seq(ebola_seq[i:]).translate())
    hiv_proteins.append(pad_seq(hiv_seq[i:]).translate())

# print("covid_protein ::", len(covid_protein))
# print("mers_protein ::", len(mers_protein))
# print("sars_protein ::", len(sars_protein))
# print("ebola_protein ::", len(ebola_protein))

# for i in range(3):
#     covid_analysed = ProteinAnalysis(str(covid_proteins[i]))
#     mers_analysed = ProteinAnalysis(str(mers_proteins[i]))
#     sars_analysed = ProteinAnalysis(str(sars_proteins[i]))
#     ebola_analysed = ProteinAnalysis(str(ebola_proteins[i]))
#     hiv_analysed = ProteinAnalysis(str(hiv_proteins[i]))

#     covid_freq = covid_analysed.count_amino_acids()
#     mers_freq = mers_analysed.count_amino_acids()
#     sars_freq = sars_analysed.count_amino_acids()
#     ebola_freq = ebola_analysed.count_amino_acids()
#     hiv_freq = hiv_analysed.count_amino_acids()

    # plt.subplot(2, 3, 1)
    # plt.bar(covid_freq.keys(), covid_freq.values())
    # plt.subplot(2, 3, 2)
    # plt.bar(mers_freq.keys(), mers_freq.values())
    # plt.subplot(2, 3, 3)
    # plt.bar(sars_freq.keys(), sars_freq.values())
    # plt.subplot(2, 3, 4)
    # plt.bar(ebola_freq.keys(), ebola_freq.values())
    # plt.subplot(2, 3, 5)
    # plt.bar(hiv_freq.keys(), hiv_freq.values())
    # plt.show()

cov_n_sars = pairwise2.align.globalxx(
    covid_seq, sars_seq, one_alignment_only=True, score_only=True
)
print(cov_n_sars)
print(cov_n_sars / len(covid_seq) * 100)

cov_n_mers = pairwise2.align.globalxx(
    covid_seq, mers_seq, one_alignment_only=True, score_only=True
)
print(cov_n_mers)
print(cov_n_mers / len(covid_seq) * 100)

cov_n_ebola = pairwise2.align.globalxx(
    covid_seq, ebola_seq, one_alignment_only=True, score_only=True
)
print(cov_n_ebola)
print(cov_n_ebola / len(covid_seq) * 100)

cov_n_hiv = pairwise2.align.globalxx(
    covid_seq, hiv_seq, one_alignment_only=True, score_only=True
)
print(cov_n_hiv)
print(cov_n_hiv / len(covid_seq) * 100)

