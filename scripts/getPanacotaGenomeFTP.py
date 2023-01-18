# turns a file containing things like this: /well/shaw/users/amu125/projects/restriction-sites/downloads/Entero_test/prepare_out/tmp_files/GCF_004684365.1_ASM468436v1_genomic.fna_prepare-split5N.fna
# into one containing things like this: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/684/365/GCF_004684365.1_ASM468436v1/GCF_004684365.1_ASM468436v1_genomic.fna.gz 

import re
from math import floor
import sys
import argparse

def get_options():
    parser = argparse.ArgumentParser(description='converts file of panacota genome links for future reference.',
                                     prog='getPanacotaGenomeFTP')
    parser.add_argument('--input', help='input file from panacota: the list of genomes in $prepare_out_dir/LSTINFO-NA-filtered-0.0001_0.06.txt', required=True)
    parser.add_argument('--output', help='output file for ftp links', required=True)
    return parser.parse_args()


def convertPanacotaString(panacota_string):
    '''converts string to ftp link'''
    genome = re.sub('_genomic.*', '', re.sub('.*\\/','', panacota_string))
    general_path = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/'
    genome_id = ''.join([re.sub('\\..*', '', x) for x in genome.split('_')[:2]])
    specific_path = '/'.join([genome_id[3*i:(3*i+3)] for i in range(0, floor(len(genome_id)/3))])
    combined_path = general_path+specific_path+'/'+genome+'/'+genome+'_genomic.fna.gz'
    return combined_path 

# test_string = '/well/shaw/users/amu125/projects/restriction-sites/downloads/Entero_test/prepare_out/tmp_files/GCF_004684365.1_ASM468436v1_genomic.fna_prepare-split5N.fna'
#print(convertPanacotaString(test_string))

def main():
    args = get_options()
    with open(args.output, 'w') as output:
        for i, line in enumerate(open(args.input, 'r').readlines()):
            if i > 0:
                genome_string = line.strip('\n').split('\t')[0]
                output.write('%s\n' % convertPanacotaString(genome_string))

if __name__ == "__main__":
    main()

