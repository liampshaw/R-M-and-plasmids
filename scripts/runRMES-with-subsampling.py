# Takes 3 arguments from commandline:
# 1 - fasta file
# 2 - output table
# 3 - word length desired
import subprocess
import os
import sys
import pandas as pd
import re
import itertools as iter
from Bio import SeqIO
import argparse
import time
import random
from math import floor

# For absolute paths when running from a different folder
_ROOT = os.path.abspath(os.path.dirname(__file__))

# Set the location of rmes on rescomp
rmes_string = '/well/shaw/users/amu125/programs/rmes-master/rmes-build/src/rmes'


def get_options():
    parser = argparse.ArgumentParser(description='Run RMES on a single fasta file.',
                                     prog='runRMES')
    parser.add_argument('--fasta', help='fasta file', required=True)
    parser.add_argument('--output', help='output csv', required=True)
    parser.add_argument('--k', help='k-mer length', required=True)
    parser.add_argument('--bases_subsample', help='subsampling size. default=0 (no subsampling)', required=False, default=0)
    parser.add_argument('--seed', help='seed (for chunk selection)', required=False, default=12345)
    return parser.parse_args()

def makeCompatibleFasta(fasta_file, compatible_fasta_output):
    '''Convert fasta to be compatible with RMES.'''
    seqs = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
   
    with open(compatible_fasta_output, 'w') as f:
        f.write('>1\n')
        for name, record in seqs.items():
                f.write('%sZ' % str(record.seq))
    return
                
            

def makeTmpFile(file_path, suffix, prefix='TMP', randomize=True):
    '''Makes a TMP_ file from a given file descriptor, taking path into account
    so that TMP_ file will be in same directory. Randomize adds a random integer 
    to the tmp filename which is handy for running on same file with different k
    at same time (e.g. on cluster job).'''
    file_path = str(file_path)
    if randomize==True:
        integer_random = str(random.randint(0, 999999))
    if '/' in file_path:
        file_str = re.sub('.*/', '', file_path)
        preamble_str = file_path[:-len(file_str)]
        tmp_path = preamble_str+'TMP_'+integer_random+'_'+file_str+'.'+suffix
    else:
        tmp_path = 'TMP_'+integer_random+'_'+file_path+'.'+suffix
    # check for byte encoding character
    tmp_path = re.sub('\ufeff', '', tmp_path)
    return tmp_path

def runRMES(fasta_file, output_table, word_length):
    '''Runs RMES on a fasta file.''' 
    tmp_fasta = makeTmpFile(fasta_file, suffix='_compatible.fa')
    
    makeCompatibleFasta(fasta_file, tmp_fasta) # Make fasta RMES compatible
    tmp_file = makeTmpFile(fasta_file, 'out')
    #print('Running rmes...')
    rmes_command = [rmes_string,
                    '--gauss', '-s', tmp_fasta,
                    '-l', str(word_length),
                    '--max', '-o', tmp_file]
    rmes_process = subprocess.Popen(rmes_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE) # Runs the process
    rmes_process.wait()
    # Try to see if the expected output file exists - if not, run the command again?
    time.sleep(20) # ait 20 seconds just to be sure
    results = rmesFormat(tmp_file+'.0', int(word_length))
    results.to_csv(output_table, index_label='word')
    print(results)
    # Remove file
    if os.path.exists(tmp_file+'.0'):
        os.remove(tmp_file+'.0')
    else: # run command again if results didn't work
        rmes_process = subprocess.Popen(rmes_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE) # Runs the process
        rmes_process.wait()
        time.sleep(20)
        results = rmesFormat(tmp_file+'.0', int(word_length))
        results.to_csv(output_table, index_label='word')

    #if os.path.exists(tmp_fasta):
    #    os.remove(tmp_fasta)
    if os.path.exists(tmp_file+'.0'):
        os.remove(tmp_file+'.0')
    return

def rmesFormat(raw_output, word_length=5):
    '''Formats the raw output from rmes (because rmes.format doesn't output where a motif has 0 expectation and 0 observations).'''
    line_block = round((4**word_length) / 10) + 1 # how many lines for a block of info, given 10 entries per line, and a gap between entries
    blank_lines = 0
    possible_mers = [''.join(x) for x in list(iter.product('agct', repeat=word_length))]
    results_dict = {k: [] for k in possible_mers}
    start_of_entry_block = 0
    with open(raw_output, 'r') as input:
        for i, line in enumerate(input.readlines()):
            counter = i-7 # correct
            if counter >= 0: # first 7 lines never have useful input
                if line.strip('\n')=='':
                    start_of_entry_block = counter+1 # store this
                    pass
                split_line = line.strip('\n').split(' ')
                if len(split_line)==1 and split_line[0]=='':
                    pass
                else:
                    norm_counter = (counter-start_of_entry_block) % line_block
                    adds =  possible_mers[(norm_counter*10):(norm_counter*10)+10]
                    for n, add in enumerate(adds):
                        results_dict[add].append(split_line[n])
    # sort the results by the scores
    results_df = pd.DataFrame.from_dict(results_dict, orient='index')
    print(results_df)
    results_df.columns = ['count', 'expect', 'sigma2', 'score']
    # convert types
    results_df = results_df.astype(dtype= {"count":"int64",
        "expect":"float64","sigma2":"float64", "score":"float64"})
    # Fill NaN (where expect=0 and count=0) with score=0
    results_df['score'] = results_df['score'].fillna(0)
    # Sort by score
    results_df.sort_values(by='score', inplace=True)
    results_df = results_df.assign(rank=range(1, len(results_df)+1))
    return results_df

def splitFasta(input_fasta, split_length=1000):
    '''Splits a fasta into chunks of length split_length. Assumes there is only one seq in fasta.
    Args:
        input_fasta (str)
            Filename of input fasta
        split_length (int)
            Size of desired chunks. Will return as many chunks of this exact length as possible
            (i.e. some sequence at end from a chunk with length<split_length may be not returned)
    Returns:
        chunks (list)
            List of chunks stored as str
    '''
    seqs = SeqIO.to_dict(SeqIO.parse(input_fasta, 'fasta')) 
    if len(seqs)==1:
        for record in seqs:
            sequence = str(seqs[record].seq)
            chunks = [sequence[(split_length*i):(split_length*(i+1))] for i in range(floor(len(sequence)/split_length))] # floor() is for as many whole chunks as possible
        return(chunks)
    else:
        return


def readRmesTbl(input_tbl):
    '''Reads in RMES tbl output as a data frame.'''
    header_names = ['word', 'count', 'expect', 'sigma2', 'score', 'rank']
    output_df = pd.DataFrame(columns=header_names)
    i = 0
    for line in open(input_tbl, 'r').readlines():
        if i>7: # 9th line has actual input
            input = [re.sub(' ', '', x) for x in line.strip('\n').split('|')]
            if len(input)==6:
                local_df = pd.DataFrame(columns=header_names, data=[input])
                output_df = output_df.append(local_df, ignore_index=True)
        i += 1
    return(output_df)

def main():
    args = get_options()
    fasta = str(args.fasta).encode(errors='ignore').decode('utf-8')
    output = str(args.output).encode(errors='ignore').decode('utf-8')
    bases_subsample = int(str(args.bases_subsample).encode(errors='ignore').decode('utf-8'))
    k = int(args.k)
    seed = int(str(args.seed).encode(errors='ignore').decode('utf-8'))
    tmp_fasta_compatible = makeTmpFile(fasta, suffix='_compatible.fa')
    makeCompatibleFasta(fasta, tmp_fasta_compatible) # Make fasta RMES compatible (in case of multiple sequences in file) 
    if bases_subsample==0: # Define '0' as 'no subsampling' 
        runRMES(tmp_fasta_compatible, output, k)
    else:
        sequence_chunks = splitFasta(tmp_fasta_compatible, split_length=bases_subsample) 
        # print(sequence_chunks)
        if len(sequence_chunks)==0: # if no chunks of length (i.e. bases_subsample>total_sequence_length) then do nothing except create empty output file
            with open(output, 'w') as f:
                f.write('')
            
        else:
            # use a random chunk specified by seed
            chunk_to_use = sequence_chunks[(len(sequence_chunks) % seed) -1]
            tmp_chunk_fasta = makeTmpFile(fasta, suffix='_chunk.fa')
            with open(tmp_chunk_fasta, 'w') as f:
                f.write('>1\n%s' % chunk_to_use) 
            runRMES(tmp_chunk_fasta, output, k)
            if os.path.exists(tmp_chunk_fasta):
                os.remove(tmp_chunk_fasta)
    #if os.path.exists(tmp_fasta_compatible):
    #    os.remove(tmp_fasta_compatible) 

if __name__ == "__main__":
    main()
