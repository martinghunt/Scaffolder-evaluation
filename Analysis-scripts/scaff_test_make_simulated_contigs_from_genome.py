#!/usr/bin/env python3.3

import argparse
from fastaq import *

parser = argparse.ArgumentParser(
    description = 'Used to generate a fake set of contigs from a genome. At regular intervals it puts in a gap and then breaks into contigs',
    usage = '%(prog)s <infile> <gap length> <contig length> <outfile>')
parser.add_argument('infile', help='Name of fasta/q file to be read')
parser.add_argument('gap_length', type=int, help='Length of gaps to be added')
parser.add_argument('contig_length', type=int, help='Length of each contig')
parser.add_argument('outfile', help='Name of output fasta file')
options = parser.parse_args()

seq_reader = sequences.file_reader(options.infile)
f_out = utils.open_file_write(options.outfile)

for seq in seq_reader:
    if len(seq) < 2 * options.contig_length + options.gap_length:
        print('Sequence', seq.id, 'too short (', len(seq), 'bases). Skipping', file=sys.stderr)


    i = 0

    while i + options.contig_length < len(seq):
        contig = sequences.Fasta(seq.id + ':' + str(i+1) + '-' +  str(i+options.contig_length), seq[i:i+options.contig_length])
        print(contig, file=f_out)
        i += options.contig_length + options.gap_length

utils.close(f_out)

