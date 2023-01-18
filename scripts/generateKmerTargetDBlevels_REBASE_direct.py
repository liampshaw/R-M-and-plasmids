# generateKmerTargetDB.py
# Combines the outputs of running rmsFinder on lots of species
# into databases of k-mers (k=4,5,6) that can then be used as
# input for calculateAvoidanceKmerGroupsAcrossDataset.py

import glob # for finding files
import pandas  
import itertools as iter
import pickle

import subprocess
import os
import sys
import pandas as pd
import re
import argparse

from math import floor
import itertools as iter
from numpy import median
from random import sample
import numpy as np
from math import nan


KMER_CATEGORIES = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Palindromic']


alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputdir', help='input directory where species output is', required=True)
    parser.add_argument('--outputdir', help='output directory', required=True)
    return parser.parse_args()

def reverse_complement(seq):
    for k,v in alt_map.items():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.items():
        bases = bases.replace(v,k)
    return bases

def all_palindromes(N):
    '''Returns all palindromes of length N.'''
    palindromes = []
    for seq in iter.product('ATCG', repeat=N):
        seq = ''.join(seq)
        half_seq = seq[0:int(N/2)]
        end_seq = seq[int(N/2):int(N)]
        if half_seq==reverse_complement(end_seq):
            palindromes.append(seq)
    return palindromes

def allCombinationsBases(N):
    '''Returns all combinations of N bases.'''
    seqs = []
    for seq in iter.product('atcg', repeat=N):
        seq = ''.join(seq)
        seqs.append(seq)
    return(seqs)

def targetedBySpecies(target_rs, species, target_rs_dict):
    '''returns 1 if yes, 0 if no)'''
    all_rs = list(set([x for y in target_rs_dict.values() for x in y]))
    if target_rs in all_rs:
        if target_rs in target_rs_dict[species]:
            return(1)
    return(0)

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

def readTargetedRS(files):
    '''returns dict of targeted RS across all files'''
    target_rs_dict = {}
    all_rs = []
    for f in files:
        species = f.split('/')[-2]
        target_rs = []
        for i, line in enumerate(open(f, 'r').readlines()):
            if i!=0:
                split_line = line.strip('\n').split(',')
                target_rs.append(split_line[1]) # rmsFinder output has second entry as target sequence
        target_rs = set(target_rs) # just presence/absence
        nonambiguous_target_rs = [] # check for ambiguous bases, using ambiguity
        for rs in target_rs:
            possible_rs = possibleSequences(rs)
            for nonambiguous_rs in possible_rs:
                nonambiguous_target_rs.append(nonambiguous_rs)
        target_rs = set(nonambiguous_target_rs)
        all_rs.append(list(target_rs))
        target_rs_dict[species] = list(target_rs) # store in dict for each species
    return(target_rs_dict)


# Now convert this to a matrix of targeted RS
def saveKmerDB(word_length, target_rs_dict, output_file):
    '''convert to matrix of targeted RS and save'''
    possible_kmers = [''.join(x) for x in list(iter.product('AGCT', repeat=word_length))]
    palindromes = all_palindromes(word_length)
    kmer_species_mat = np.zeros( (len(possible_kmers), len(target_rs_dict)) , dtype=int) # matrix of zeroes
    for k, kmer in enumerate(possible_kmers):
        for j, species in enumerate(sorted(target_rs_dict.keys())):
            score = targetedBySpecies(kmer, species, target_rs_dict)
            kmer_species_mat[k, j] = score
    kmer_species_df = pd.DataFrame(kmer_species_mat, columns = sorted(target_rs_dict.keys()), index = possible_kmers)
    kmer_species_df.to_csv(output_file)

# Use this dataframe for a given species to return the five groups of k-mers
# within-species targeted RS, within-genus targeted RS, other-genus targeted RS (somewhere in REBASE), not-targeted-but-palindromic RS, other k-mers.
def kmerGroupsForSpecies(species, kmer_species_df, word_length): #, all_rec_seqs_file):
    '''return k-mer groupings for a given species, given the kmer_species_df'''
    possible_kmers = [''.join(x) for x in list(iter.product('AGCT', repeat=word_length))]
    possible_kmer_groups = {k:[] for k in possible_kmers} # start all off without group
    palindromes =  all_palindromes(word_length)     #all_rec_seqs = [line.strip('\n') for line in open(all_rec_seqs_file, 'r')]
    genus = re.sub('_.*', '', species)
    genus_df = kmer_species_df.filter(regex=(genus+'*'))
    for k in possible_kmers:
        if genus_df.loc[k, species]==1:
            possible_kmer_groups[k].append('within-species')
        if sum(genus_df.loc[k])>0:
            possible_kmer_groups[k].append('within-genus')
        elif sum(kmer_species_df.loc[k])>0:
            possible_kmer_groups[k].append('outside-genus-found-in-dataset')
        if k in palindromes:
            possible_kmer_groups[k].append('palindromic')
        if len(possible_kmer_groups[k])==0:
            possible_kmer_groups[k].append('other')
    return possible_kmer_groups

def getSpeciesKmerGroups(species_list, output_dir, word_length):
    '''gets all the species kmer groups and returns them'''
    pan_species_kmer_dict = {}
    for species in species_list:
        species_kmers_file = output_dir+species+'-'+str(word_length)+'mer-categories.txt'
        groups_dict = {}
        with open(species_kmers_file, 'r') as f:
            for line in f.readlines():
                line = line.strip('\n').split(' ')
                groups_dict[line[0]] = line[1:]
        pan_species_kmer_dict[species] = groups_dict
    return(pan_species_kmer_dict)

def getOtherSpeciesForTaxonomicLevel(species, taxonomy_df, level):
    '''gets all other species sharing a particular taxonomic level with another species'''
    level_name = taxonomy_df.loc[species, level]
    return list(taxonomy_df[taxonomy_df[level]==level_name]['Species'])

def classifyKmerForSpecies(kmer, species, kmer_species_df, taxonomy_df):
    '''classifies a kmer for a given species: is it seen in species? then in genus? etc.'''
    taxonomic_levels = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Kingdom']
    # For the k-mer easiest to do all the levels in a loop?
    level_to_return = 'unknown'
    for i, level in enumerate(taxonomic_levels):
        species_to_examine = getOtherSpeciesForTaxonomicLevel(species, taxonomy_df, level)
        #print(level, ':', species_to_examine)
        if any([kmer_species_df.loc[kmer, x] for x in species_to_examine]):
            level_to_return = level
            break
        else:
            pass
    return(level_to_return) 


def main():
    args = get_options()
    input_dir = str(args.inputdir)+'/' # trailing slash

    print(input_dir)
    output_dir = str(args.outputdir)+'/' # trailing slash just in case

    # Values of k to look at
    k_values = [4, 5, 6]




    # Read in taxonomy
    taxonomy_df = pd.read_csv(input_dir+'/taxonomy-filtered.txt', sep=';', header=None)
    taxonomy_df.columns = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    taxonomy_df['Species'] = [re.sub(' ', '_', x) for x in list(taxonomy_df['Species'])]
    taxonomy_df.index = taxonomy_df['Species']

    for k in k_values:
        print(k)
        kmer_species_df = pd.read_csv(output_dir+str(k)+'-kmer-db.csv', index_col=0)
        #kmer_species_df = pd.read_csv(output_dir+str(k)+'-kmer-db.csv', index_col=0)

        answers = {s: {} for s in taxonomy_df.index}

        kmers_of_interest = [kmer for kmer in kmer_species_df.index if sum(kmer_species_df.loc[kmer])!=0]
        palindromes = all_palindromes(k)
        groups = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Palindromic']

        for s in answers.keys():
            print(s)
            for g in groups:
                answers[s][g] = [] 
            for kmer in palindromes:
                answers[s]['Palindromic'].append(kmer)
            for kmer in kmers_of_interest: 
                answer = classifyKmerForSpecies(kmer, s, kmer_species_df, taxonomy_df)
                #print('answer', answer)
                answers[s][answer].append(kmer)
                if answer!="unknown":
                    print(s, kmer, answer)    
        with open(output_dir+str(k)+'-species-target-dict-REBASE-direct.pkl', 'wb') as db_handle:
            pickle.dump(answers, db_handle)

if __name__ == "__main__":
    main()
