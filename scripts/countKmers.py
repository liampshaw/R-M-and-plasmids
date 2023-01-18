from itertools import product

import argparse

# Using Matthew Fahrbach solution from Rosalind problem 
# K-mer counts are enumerated in lexicographical order
# of the k-mers
# e.g. AAAA comes first for k=4


def get_options():
    parser = argparse.ArgumentParser(description='Counts k-mers in a fasta and returns in lexicographic order.',
                                     prog='countKmers')
    parser.add_argument('--fasta', help='fasta file', required=True)
    parser.add_argument('--k', help='value of k', required=True) 
    # To add: subsampling argument (take a random X bp section and count k-mers in that)
    return parser.parse_args()


def get_seq(input_fasta):
    dna = ''
    with open(input_fasta, 'r') as f:
        for line in f.readlines():
            if line.startswith('>'):
                dna += 'Z'
            else:
                dna += line.strip('\n')
    return(dna)


def count_kmers(fasta, k):
    dna = get_seq(fasta)
    table = {};
    for kmer in product('ACGT', repeat=k):
        table[''.join(kmer)] = 0

    for i in range(len(dna) - k + 1):
        kmer = dna[i:i + k]
        if all([base in ['A', 'T', 'C', 'G'] for base in kmer]): # only if kmer is ATCG-based, no other characters allowed
            table[kmer] = table[kmer] + 1
    for kmer in product('ACGT', repeat=k):
        output_string = print(table[''.join(kmer)], end=' ')
    return(output_string)


def main():
    args = get_options()
    input_fasta = args.fasta
    k = int(args.k)
    count_kmers(input_fasta, k)   
#    table = {};
#    for kmer in product('ACGT', repeat=k):
#        table[''.join(kmer)] = 0
#
#    for i in range(len(dna) - k + 1):
#        kmer = dna[i:i + k]
#        if all([base in ['A', 'T', 'C', 'G'] for base in kmer]): # only if kmer is ATCG-based, no other characters allowed
#            table[kmer] = table[kmer] + 1 
##        if not Z' in kmer:
# #           table[kmer] = table[kmer] + 1
#
#    for kmer in product('ACGT', repeat=k):
#        print(table[''.join(kmer)], end=' ')

if __name__=="__main__":
	main()
