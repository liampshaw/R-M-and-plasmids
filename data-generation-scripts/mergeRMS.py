# Merges the rmsFinder output in a folder 
# conda activate panacota_rmes_env
import glob
import pandas as pd
import re
from scipy.stats import wilcoxon
from statistics import median
import argparse
import csv
from io import StringIO 

def get_options():
    parser = argparse.ArgumentParser(description='Merges all RMS output files in a folder matching a pattern.',
                                     prog='mergeRMS')
    parser.add_argument('--dir', help='directory containing RMS output files.', required=True)
    parser.add_argument('--suffix', help='pattern to match to (suffix): RMS, MT, RE', required=True)
    parser.add_argument('--output', help='output file', required=True)
    return parser.parse_args()

def mergeTheFiles(folder, suffix, output):
    '''merges the files in a folder matching a suffix.'''
    rmes_output_files = glob.glob(folder+'/*_'+suffix+'.csv')
    all_data = [] # append data as list
    for f in rmes_output_files: 
        short_name = re.sub('_'+suffix+'.csv', '', f)
        short_name = re.sub('.*\/', '', short_name)
        for i, line in enumerate(open(f, 'r').readlines()):
            if i==0 and line==('\n'):
                pass
            if i!=0:
                data = StringIO(line.strip('\n'))
                reader = csv.reader(data, delimiter=',')
                split_line = list(reader)[0]
                split_line.insert(0, short_name)
                all_data.append(split_line)
    if suffix=='RMS':
        columns_to_use = ['genome','sequence','contig','contig_description','pos_MT','pos_RE','prot_MT','prot_RE','hit_MT','sim_MT','loc_MT','hit_RE','sim_RE','loc_RE']
    else:
        columns_to_use = ['genome','qseqid','sseqid','pident','length','qlen','evalue','coverage_threshold_met','target','contig','position','contig_description','genomic_location','similarity','hit_type','n_REBASE_hits','identicalTarget']

    all_data_df = pd.DataFrame(all_data, columns=columns_to_use)
    all_data_df.to_csv(output, index=False)


def main():
    args = get_options()
    mergeTheFiles(args.dir, args.suffix, args.output)

if __name__ == "__main__":
    main()

