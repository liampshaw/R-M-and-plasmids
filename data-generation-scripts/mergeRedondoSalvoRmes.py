# Merges the RMES output for Redondo-Salvo plasmids.
# conda activate rmsFinderEnv
import glob
import pandas as pd
import re
import argparse 

# Read in the Redondo-Salvo database
print("reading in db...")
redondo_salvo_db = pd.read_csv('/well/shaw/users/amu125/data/redondo-salvo/redondo-salvo-2020-DB.csv', index_col=0)
print("done")

def get_options():
    parser = argparse.ArgumentParser(description='Merge output for Redondo-Salvo plasmids.',
                                     prog='mergeRedondoSalvoRmes')
    parser.add_argument('--outdir', help='output directory', required=True)
    parser.add_argument('--k', help='k-mer length', required=True)
    parser.add_argument('--n', help='subsampling', required=True)
    return parser.parse_args()


# Read in RSs for each genera
def mergeRMES(base_dir, out_file_path):
    '''merges all RMES results for plasmids within a given base directory'''
    rmes_output_files = glob.glob(base_dir+'/rmes*.csv')
    all_data = [] # append data as list
    for f in rmes_output_files: 
        plasmid = re.sub('\\.csv', '', re.sub('.*\\/', '', f).split('-')[2])
        for i, line in enumerate(open(f, 'r').readlines()):
            if i!=0:
                split_line = line.strip('\n').split(',')
                split_line.insert(0, plasmid)
                all_data.append(split_line)
    all_data_df = pd.DataFrame(all_data, columns=['plasmid', 'word', 'count', 'expect', 'sigma2', 'score', 'word_rank'])
    # convert dtypes
    all_data_df['word_rank'] = all_data_df['word_rank'].astype('int')
    all_data_df['score'] = all_data_df['score'].astype('float')
    all_data_df.to_csv(out_file_path, index=False)

def main():
    args = get_options()
    k = str(args.k)
    n = str(args.n)
    out_dir = str(args.outdir)
    # Read in the Redondo-Salvo database
    print("reading in db...")
    redondo_salvo_db = pd.read_csv('/well/shaw/users/amu125/data/redondo-salvo/redondo-salvo-2020-DB.csv', index_col=0)
    print("done")
#base_dir = '/well/shaw/users/amu125/projects/restriction-sites/data/plasmids/rmes/k-'+str(k)+'/'
    #out_dir = '/well/shaw/users/amu125/projects/restriction-sites/data/plasmids/rmes/'
#mergeRMES(base_dir, out_dir+'redondo-salvo-k-'+str(k)+'.csv') # do it for the unsubsampled first
    for subsampling_dir in glob.glob('/well/shaw/users/amu125/projects/restriction-sites/data/plasmids/rmes/n-'+n+'/k-'+k+'/'):
        base_dir = subsampling_dir
        #subsampling_n = re.sub('.*n-', '', subsampling_dir) # get the size of subsampling
        mergeRMES(base_dir, out_dir+'redondo-salvo-k-'+k+'-n-'+n+'.csv')
    
if __name__=="__main__":
    main() 
