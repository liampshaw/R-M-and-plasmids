# Read in core genome file and then produce, for each genome,

# Conda environment: unclear

# Libraries
from Bio import SeqIO
import glob
import re
import argparse
import os

def get_options():
    parser = argparse.ArgumentParser(description='Create core genome nucleotide fasta from PanACoTA output.', prog='makeCoreGenomeFastas')
    parser.add_argument('--core_genome', help='Core genome file (main output of PanACoTA corepers)', required=True)
    parser.add_argument('--gene_dir', help='Directory where genes are stored (usually annotate_output/Genes)', required=True)
    parser.add_argument('--out_dir', help='Output directory for core genome fastas.', required=True)
    return parser.parse_args()

def subsetFasta(input_fasta, seq_names, output_fasta, invert=False):
    '''Subsets a fasta file and creates a new fasta.
    Args:
        input_fasta (str)
            The original fasta file
        seq_names (list)
            The list of sequence names to be subsetted out
        output_fasta (str)
            The output file to write the subset to
        invert (Bool)
            When True, the function returns those sequences *not* in the subset
    Returns:
        None
    '''
    seqs = SeqIO.to_dict(SeqIO.parse(input_fasta, 'fasta'))
    if invert==False:
        subset_seqs = [seqs[record] for record in seq_names]
    if invert==True:
        subset_seqs = [seqs[record] for record in seqs.keys() if record not in seq_names]
    with open(output_fasta, 'w') as output_file:
        for record in subset_seqs:
            output_file.write('>%s\n%s\n' % (record.id, str(record.seq))) 
    return

def read_core_genome(core_genome_file):
    '''Reads in the core genome and storeds as a dict where each genome has list of genes included.
    Args:
        core_genome_file (str)
            Output of PanACoTA corepers
    Returns:
        core_genome_dict (dict)
            Dict of core genome, genomes as keys
    '''
    gene_dict = {}
    core_genome_dict = {}
    with open(core_genome_file, 'r') as f:
        for i, line in enumerate(f.readlines()):
            entries = line.strip('\n').split(' ')
            genes = entries[1:]
            gene_dict[entries[0]] = genes
            if i==0: # for first line, get the names of the genomes
                genome_names = ['.'.join(x.split('.')[:3]) for x in genes]
            for g in genes:
                genome_name = '.'.join(g.split('.')[:3])
                if genome_name in core_genome_dict.keys():
                    core_genome_dict[genome_name].append(g)
                else:
                    core_genome_dict[genome_name] = [g]
    return core_genome_dict


def readGenes():
    files = glob.glob()



def main():
    args = get_options()
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
    core_genome_dict = read_core_genome(args.core_genome)
    for genome in core_genome_dict.keys():
        core_genes = core_genome_dict[genome]
        genome_fasta = args.gene_dir+'/'+genome+'.gen'
        genome_number = re.sub('.*\\.', '', genome)
        all_genes = SeqIO.to_dict(SeqIO.parse(genome_fasta, 'fasta')).keys()
        
        # Three categories of gene: core, accessory (chromosome), and plasmid
        core_chromosome_genes = [gene for gene in core_genes if genome_number+'.0001' in gene]
        accessory_chromosome_genes = [gene for gene in all_genes if gene not in core_genes and genome_number+'.0001' in gene] # Contig number (assuming complete genome) indicates chromosome
        accessory_plasmid_genes = [gene for gene in all_genes if gene not in core_genes and genome_number+'.0001' not in gene] # Contig>0001 indicates plasmid (assuming complete genome)
        # Make core genome fasta
        subsetFasta(input_fasta = args.gene_dir+'/'+genome+'.gen',
                    seq_names = core_chromosome_genes,
                    output_fasta = args.out_dir+'/'+genome+'.core')
        # Make accessory (chromosome) genome fasta
        subsetFasta(input_fasta = args.gene_dir+'/'+genome+'.gen',
                    seq_names = accessory_chromosome_genes,
                    output_fasta = args.out_dir+'/'+genome+'.accessory_chrom')
        # Make accessory (plasmid) genome fasta
        subsetFasta(input_fasta = args.gene_dir+'/'+genome+'.gen',
                    seq_names = accessory_plasmid_genes,
                    output_fasta = args.out_dir+'/'+genome+'.accessory_plas')


if __name__ == "__main__":
    main()
