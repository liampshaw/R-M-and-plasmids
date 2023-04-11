# Merges the RMES output in a folder 
# conda activate panacota_rmes_env
import glob
import pandas as pd
import re
from scipy.stats import wilcoxon
from statistics import median
import argparse

def get_options():
    parser = argparse.ArgumentParser(description='Merges all RMES files in a folder matching a pattern.',
                                     prog='mergeRMES')
    parser.add_argument('--dir', help='directory containing RMES output files.', required=True)
    parser.add_argument('--suffix', help='pattern to match to (suffix)', required=True)
    parser.add_argument('--output', help='output file', required=True)
    return parser.parse_args()

def mergeTheFiles(folder, suffix, output, compression=True):
    '''merges the files in a folder matching a suffix.'''
    rmes_output_files = glob.glob(folder+'/*'+suffix)
    all_data = [] # append data as list
    for f in rmes_output_files: 
        short_name = re.sub(suffix, '', f)
        for i, line in enumerate(open(f, 'r').readlines()):
            if i!=0:
                split_line = line.strip('\n').split(',')
                split_line.insert(0, short_name)
                all_data.append(split_line)
    all_data_df = pd.DataFrame(all_data, columns=['name', 'word', 'count', 'expect', 'sigma2', 'score', 'word_rank'])
    # convert dtypes
    all_data_df['word_rank'] = all_data_df['word_rank'].astype('int')
    all_data_df['score'] = all_data_df['score'].astype('float')
    if compression==True:
        all_data_df.to_csv(output+'.gz', index=False, compression='gzip')
    else:
        all_data_df.to_csv(output, index=False)


def main():
    args = get_options()
    mergeTheFiles(args.dir, args.suffix, args.output)

if __name__ == "__main__":
    main()

