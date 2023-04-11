# Stores the statistics on all fasta files in a directory
# (for storing statistics on makeCoreGenomeFastas so that can use
# when looking at subsampling: what proportion of genome component is 
# being subsampled?)

import glob
import sys
import re

input_dir = str(sys.argv[1])
possible_files = glob.glob(input_dir+'*')

for filename in possible_files:
    total_length = 0
    total_CDS = 0
    total_GC = 0
    contig_names = []
    for line in open(filename, 'r').readlines():
        if not line.startswith('>'):
            total_length += len(line.strip('\n'))
            total_GC += line.count('G') + line.count('C') 
        elif line.startswith('>'):
            total_CDS += 1
            contig_name = re.sub('[a-z]', '', line.split('_')[0].split('.')[3])          
            if contig_name not in contig_names:
                contig_names.append(contig_name)  
    short_filename = re.sub('.*\/', '', filename)
    if total_length!=0:
        total_GC_content = float(total_GC)/float(total_length)
    else:
        total_GC_content = 'NA'
    print('.'.join(short_filename.split('.')[0:3]), short_filename.split('.')[3], len(contig_names),total_length, total_CDS, total_GC_content)
