# Gets all the RecSeqs from a REBASE fasta file
from Bio import SeqIO
from Bio import Entrez
import sys
import re
import argparse
import itertools as iter
import numpy as np
import pandas as pd

# Set NCBI Entrez email address
Entrez.email = "liam.philip.shaw@gmail.com"

def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='input REBASE fasta', required=True)
    parser.add_argument('--outputprefix', help='output file prefix', required=True)
    parser.add_argument('--verbose', help='whether to print output or not', required=False, action='store_true')
    return parser.parse_args()

ambiguity_codes =   {'A': ['A'],\
                    'G' : ['G'],\
                    'C'	: ['C'],\
                    'T' : ['T'],\
                    'Y'	: ['C', 'T'],\
                    'R'	: ['A','G'],\
                    'W'	: ['A','T'],\
                    'S'	: ['G','C'],\
                    'K'	: ['T','G'],\
                    'M'	: ['C','A'],\
                    'D'	: ['A','G','T'],\
                    'V'	: ['A','C','G'],\
                    'H'	: ['A','C','T'],\
                    'B'	: ['C','G','T'],\
                    'X'	: ['A','C','G','T'],\
                    'N'	: ['A','C','G','T'],\
                    '-' : ['-']}

def targetedBySpecies(target_rs, species, target_rs_dict):
    '''returns 1 if yes, 0 if no)'''
    all_rs = list(set([x for y in target_rs_dict.values() for x in y]))
    if target_rs in all_rs:
        if target_rs in target_rs_dict[species]:
            return(1)
    return(0)

def possibleSequences(dna_sequence):
    '''Takes a DNA sequence (possibly containing ambiguous bases) and returns
    a list of all possible sequences.
    Args:
        dna_sequence (str)
            String of DNA sequence (uppercase)
    Returns:
        possible_strings (list)
            List of strings - all possible DNA sequences
    '''
    # Get all possible bases at each position
    possible_bases = [ambiguity_codes[base] for base in dna_sequence]
    # Go through and add these on
    possible_strings = []
    for bases in possible_bases:
        if len(possible_strings)==0: # if no strings yet, use first base
            possible_strings = list(bases)
        elif len(bases)==1: # if just one possible base at position, add it on
            possible_strings = [x+bases[0] for x in possible_strings]
        else: # if multiple possibilities, add them on one by one
            additions = []
            for base in bases:
                additions.append([x+base for x in possible_strings])
            possible_strings = [x for addition in additions for x in addition]
    return possible_strings

def getRecSeqs(input_file, k=[4,5,6]):
    '''gets all the recseqs for a range of k from an input_file (REBASE format)
    including with their GenBank accession'''
    recseqs = []
    for line in open(input_file, 'r').readlines():
        if line.startswith('>'):
            if 'RecSeq' in line:
                splitline = line.strip('\n').split('\t')
                recseq_item = [x for x in splitline if 'RecSeq' in x][0]
                genbank_item = [x for x in splitline if 'GenBank' in x][0]
                recseqs.append([genbank_item.split(':')[1], recseq_item.split(':')[1]])
    #recseqs_all = []
    #recseqs_set = sorted(list(set(recseqs)))
    recseqs_all = []
    for pair in recseqs: # account for ambiguity
        genbank = pair[0]
        recseq = pair[1]
        recseq_list = recseq.split(',') # some have multiple (not many)
        for recseq in recseq_list:
            recseq = re.sub(' ', '', recseq)
            if len(recseq) in [4,5,6]: # only interested in k=4,5,6
                possible_recseqs = possibleSequences(recseq)
                for possible_recseq in possible_recseqs:
                    recseqs_all.append([genbank, recseq, possible_recseq])
    #recseqs_all = sorted(list(set(recseqs_all))) # convert to set (presence/absence only)
    return(recseqs_all)

def searchEntrez(genbank_id):
    '''searches entrez and returns the species of the genbank id'''
    handle = Entrez.esummary(db="nucleotide", id=genbank_id, retmode="xml")
    try:
        record = Entrez.read(handle)
        species_name = ' '.join(record[0]['Title'].split(' ')[0:2])
        return(species_name)
    except:
        return('Error')


def main():
    '''runs for RE and MT, then combines'''
    args = get_options()
    input_file = str(args.input)
    recseqs_list = getRecSeqs(str(args.input))
    with open(str(args.outputprefix)+'_all_recseqs.csv', 'w') as output_file:
        for i, recseq_pair in enumerate(recseqs_list):
            genbankid = recseq_pair[0]
            species = searchEntrez(genbankid)
            recseqs_list[i].append(species)
            output_file.write(str(args.input)+','+','.join(recseqs_list[i])+'\n')
            if args.verbose:
                print(str(i)+','+','.join(recseqs_list[i]))

    # Also convert to a matrix format for k=4,5,6
    species_names = set([x[3] for x in recseqs_list])
    targeted_rs_dict = {} # Dictionary of species and the recseqs it targets
    for species in species_names:
        targeted_rs_dict[species] = list([x[2] for x in recseqs_list if x[3]==species])

    for k in [4,5,6]:
        possible_kmers = [''.join(x) for x in list(iter.product('AGCT', repeat=k))]
        kmer_species_mat = np.zeros( (len(possible_kmers), len(species_names)) , dtype=int) # matrix of zeroes
        for i, kmer in enumerate(possible_kmers):
            for j, species in enumerate(sorted(targeted_rs_dict.keys())):
                score = targetedBySpecies(kmer, species, targeted_rs_dict)
                kmer_species_mat[i, j] = score # create species matrix
        kmer_species_df = pd.DataFrame(kmer_species_mat, columns = sorted(targeted_rs_dict.keys()), index = possible_kmers)
        kmer_species_df.to_csv(str(args.outputprefix)+'-'+str(k)+'-kmer-db.csv') # write to csv

if __name__ == "__main__":
    main()
