# Combines the outputs of running rmsFinder on lots of species
# into databases of k-mers (k=4,5,6)
# And then calculates the avoidance based on those groups

import glob # for finding files
import pandas as pd
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
#import generate_palindromes as gp
from numpy import median
from random import sample
import numpy as np
from math import nan

import generateKmerTargetDBlevels as gk

def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputdir', help='input directory where rmes output is stored', required=True)
    parser.add_argument('--k', help='value of k to use (rmes results)', required=True)
    parser.add_argument('--targetdict', help='pickle file of species k-mer targeting', required=True)
    parser.add_argument('--output', help='output csv file', required=True)
    parser.add_argument('--inclusive', help='whether to treat kmer-hierarchy inclusively', action='store_true', default=False)
    return parser.parse_args()

def summariseRmes(df, groups_dict):
    '''summarise rmes results for a (pre-filtered) dataframe and a kmer grouping'''
    df.index = df['word_rank'] # give ranks as row names for indexing
    median_results = {group: [] for group in gk.KMER_CATEGORIES}
#    groups_dict = {group: [] for group in gk.KMER_CATEGORIES}
#    if inclusive==False:
#        # Make lower dict
#        for group in groups_dict.keys():
#            for kmer in kmers_dict[group]: 
#                groups_dict[group].append(kmer.lower()) # need lower because of rmes results
#    elif inclusive==True:
#        _ = [groups_dict['Palindromic'].append(kmer.lower()) for kmer in kmers_dict['Palindromic']] # need lower due to rmes results
#        cumulative_kmers = kmers_dict['Species']
#        _ = [groups_dict['Species'].append(kmer.lower()) for kmer in cumulative_kmers]
#        for group in ['Genus', 'Family', 'Order', 'Class', 'Phylum', 'Kingdom']:
#            cumulative_kmers += kmers_dict[group]
#            _ = [groups_dict[group].append(kmer.lower()) for kmer in cumulative_kmers] # need lower due to rmes 

    for group, kmers in groups_dict.items():
        filtered_df = df[df['word'].isin(kmers)]
        if len(filtered_df)==0:
            median_rank = nan
            median_score = nan
            n = 0
        else:
            median_rank = median(filtered_df['word_rank'])           
            median_score = median(filtered_df['score'])
            n = len(filtered_df)
            #print(group, n, median_score)
        median_results[group].append(n)
        median_results[group].append(median_rank)
        median_results[group].append(median_score)
    try:
        return([median_results[x] for x in gk.KMER_CATEGORIES]) # Return in the hierarchical order expected
    except:
        return(None)

# Read in file and get median
def summariseRmesSpecies(species, dir, groups_dict, suffix, word_length):
    '''summarise rmes results for a species
    species: name of species
    suffix: e.g. -0-all (.csv included by default)
    '''
    results = pd.read_csv(dir+'/'+species+'/rmes-k-'+str(word_length)+'-'+suffix+'.csv')
    results.index = range(len(results))
    n_words = 4**int(word_length)
    results_array = np.array(results)
    possible_categories = gk.KMER_CATEGORIES
    median_scores = {group: [] for group in possible_categories}
    #print(median_scores)
    for i in range(0, len(results), n_words):
        genome_df = pd.DataFrame(results[i:i+n_words])
        genome_df.columns = results.columns
        genome_name = list(set(genome_df['name']))
        if len(genome_name)!=1:
            print('error')
        else:
            genome_name_short = '.'.join(re.sub('.*\/', '', genome_name[0]).split('.')[:-1])
            genome_subsampling = re.sub('.*-', '', genome_name[0].split('/')[-2])
            genome_part = re.sub('.*\/', '', genome_name[0]).split('.')[-1]
            #print(genome_name)
            genome_summary = summariseRmes(genome_df, groups_dict) 
            
            #print(genome_summary)
            for i, group in enumerate(gk.KMER_CATEGORIES):
                median_scores[group].append([genome_name_short, genome_part, genome_subsampling, genome_summary[i][0], genome_summary[i][1], genome_summary[i][2]])
    return(median_scores)

def main():
    args = get_options()
    k = str(args.k)
    input_dir = str(args.inputdir)
    input_species = list(set([x.split('/')[-2] for x in glob.glob(input_dir+'/*/*csv')]))
    target_dict = str(args.targetdict)
    output_file = str(args.output)
    inclusive = args.inclusive # whether to be inclusive with k-mers or not
    print('Inclusivitiy:', inclusive)
    word_length = k
    print(k)
    with open(target_dict, 'rb') as handle:
        species_kmer_dict = pickle.load(handle)

    # To store all results
    all_results_list = []
    print('Running through the species now...')
    with open(output_file, 'w') as f:
        for species in input_species:
            print(species)
            # Get the group dict of kmer categories
            kmer_dict = species_kmer_dict[species]
            groups_dict = {group: [] for group in gk.KMER_CATEGORIES}
            if inclusive==False:
                for group in groups_dict.keys():
                    for kmer in kmer_dict[group]:
                        groups_dict[group].append(kmer.lower()) # need lower because of rmes results
            elif inclusive==True:
                _ = [groups_dict['Palindromic'].append(kmer.lower()) for kmer in kmer_dict['Palindromic']] # need lower due to rmes results
                cumulative_kmers = kmer_dict['Species']
                _ = [groups_dict['Species'].append(kmer.lower()) for kmer in cumulative_kmers]
                for group in ['Genus', 'Family', 'Order', 'Class', 'Phylum', 'Kingdom']:
                    cumulative_kmers += kmer_dict[group]
                    _ = [groups_dict[group].append(kmer.lower()) for kmer in list(cumulative_kmers)] # need lower due to rmes 
            print(groups_dict)
            # Find all the subsampling files
            subsampling_files = glob.glob(input_dir+'/'+species+'/rmes-k-'+str(k)+'*csv')
            print(subsampling_files)
            subsampling_values = [re.sub('-.*', '', re.sub('.*k-'+str(k)+'-', '', re.sub('.*\\/', '', x))) for x in subsampling_files]
            for subsampling in subsampling_values:
                print(subsampling)
                print('Computing...')
                median_scores = summariseRmesSpecies(species, input_dir, groups_dict, subsampling+'-all', word_length) 
                for kmer_group, scores in median_scores.items():
                    for score in scores:
                        result = [species, kmer_group, score[0], score[1], score[2], score[3], score[4], score[5]]
                        _ = f.write('%s\n' % ','.join([str(x) for x in result])) 
                        print(result)#all_results_list.append([species, kmer_group, score[0], score[1], score[2], score[3], score[4], score[5]])
                   
    #with open(output_file, 'w') as f:
    #    for result in all_results_list:
    #        _ = f.write('%s\n' % ','.join([str(x) for x in result]))

if __name__ == "__main__":
    main()

