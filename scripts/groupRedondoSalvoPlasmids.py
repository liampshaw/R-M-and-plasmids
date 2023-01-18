#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Group Redondo-Salvo plasmids into PTU fastas
import pandas as pd
from Bio import SeqIO
import re
import os

redondo_salvo_metadata_db = pd.read_csv('/well/shaw/users/amu125/data/redondo-salvo/redondo-salvo-2020-DB.csv', index_col=0)
redondo_salvo_metadata_db = redondo_salvo_metadata_db.assign(species=[re.sub(' ', '_', x) for x in redondo_salvo_metadata_db['TaxSpecies']])

def concatenatePlasmids(seq_dict, plasmid_ids):
    '''Concatenate a list of plasmids given a seq_dict.'''
    concatenated_seq = ''
    for plasmid in plasmid_ids:
        concatenated_seq += seq_dict[plasmid]+'Z'
    return(concatenated_seq)

def makeCompatibleFasta(fasta_file, compatible_fasta_output):
    '''Convert fasta to be compatible with RMES.'''
    seqs = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))

    with open(compatible_fasta_output, 'w') as f:
        f.write('>1\n')
        for name, record in seqs.items():
                f.write('%sZ' % str(record.seq))
    return

for PTU in set(redondo_salvo_metadata_db['PTU']):
    print(PTU)
    PTU_file_write = re.sub('\/', '-', PTU)
    PTU_accessions = redondo_salvo_metadata_db[redondo_salvo_metadata_db['PTU']==PTU].index
    with open('tmp.fa', 'w+') as f:
        for accession in PTU_accessions:
            data = open('/well/shaw/users/amu125/data/redondo-salvo/nt/'+accession+'.fa', 'r').read()
            f.write('%s' % data[:-1]) # write, but strip EOF
    makeCompatibleFasta('tmp.fa', '/well/shaw/users/amu125/data/redondo-salvo/PTU/'+PTU_file_write+'.fa')
    os.remove('tmp.fa')


