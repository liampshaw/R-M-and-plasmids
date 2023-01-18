from Bio import SeqIO
import sys

def read_fasta(fasta_file):
    '''Simply uses SeqIO to read in fasta as dict.'''
    return(SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'), lambda rec:rec.id))

seqs = read_fasta(sys.argv[1])
for seq in seqs:
    seq_name = seq.split('.')[0]
    seq_seq = str(seqs[seq].seq)
    with open(seq_name+'.fa', 'w') as f:
        f.write('>%s\n%s\n' % (seq_name, seq_seq))

