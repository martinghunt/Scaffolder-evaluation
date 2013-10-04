#!/usr/bin/env python3

import argparse
import sys
import copy
import os
from fastaq import *

parser = argparse.ArgumentParser(
    description = 'Used to generate contig test sets and reads.',
    usage = '%(prog)s <infile>')
parser.add_argument('infile', help='Name of config file')
options = parser.parse_args()

# parse the config file
f = utils.open_file_read(options.infile)
contigs = {}
scaffolds = []
ref_seqs = []
trim_bases = -1

for line in f:
    l = line.rstrip().split()
    if l[0] == 'contigs_file':
        all_contigs_infile = l[1]
    elif l[0] == 'outprefix':
        outprefix = l[1]
    elif l[0] == 'scaffold':
        scaffolds.append((l[1], l[2:]))
    elif l[0] == 'ref':
        ref_seqs.append(l[1:])
    elif l[0] == 'trim':
        trim_bases = int(l[1])
    else:
        print('error in config file. Line:', line, file=sys.stderr)
        sys.exit(1)

utils.close(f)
assert trim_bases >= 0

reads_fastq_files = []
contig_seqs = {}
contigs_outfile = outprefix + '.contigs.fa'

# store contigs in memory
tasks.file_to_dict(all_contigs_infile, contig_seqs)

# make trimmed contigs
trimmed_contig_seqs = {}
for id, fa in contig_seqs.items():
    contig = copy.copy(fa)
    contig.trim(trim_bases, trim_bases)
    trimmed_contig_seqs[id] = contig


# write file of contigs to be scaffolded, alphabetical order.
contigs_for_scaffolding = set()

for scaff in scaffolds:
    contigs_for_scaffolding.update(set(scaff[1]))

contigs_for_scaffolding = list(contigs_for_scaffolding)
contigs_for_scaffolding.sort()

f = utils.open_file_write(contigs_outfile)
for c in contigs_for_scaffolding:
    print(trimmed_contig_seqs[c], file=f)
utils.close(f)

# make scaffolds and simulate reads from them
for i in range(len(scaffolds)):
    coverage, contig_list = scaffolds[i]
    scaff_name = 'scaff.' + str(i+1)
    scaff_fname = outprefix + '.' + scaff_name + '.fa'
    seq = sequences.Fasta(scaff_name, ('').join([contig_seqs[c].seq for c in contig_list]))
    f = utils.open_file_write(scaff_fname)
    print(seq, file=f)
    utils.close(f)
    reads_fname = scaff_fname + '.reads.fq'
    reads_fastq_files.append(reads_fname)
    cmd = 'fastaq_to_perfect_reads ' + scaff_fname + ' ' + reads_fname + ' 500 30 ' + coverage + ' 76'
    utils.syscall(cmd)
    os.unlink(scaff_fname)

# cat all the reads files together
reads_fastq = outprefix + '.reads.fq'
fout = utils.open_file_write(reads_fastq)
for fname in reads_fastq_files:
    with open(fname) as infile:
        for line in infile:
            fout.write(line)

    os.unlink(fname)

utils.close(fout)

# make deinterleaved fastq files
tasks.deinterleave(reads_fastq, outprefix + '.reads_1.fq', outprefix + '.reads_2.fq')

# make ref sequences
ref_fasta = outprefix + '.ref.fa'
f = utils.open_file_write(ref_fasta)

for l in ref_seqs:
    seq = sequences.Fasta(''.join(l), ('N'*2*trim_bases).join([trimmed_contig_seqs[c].seq for c in l]))
    print(seq, file=f)

utils.close(f)

# index the reference and generate unique tag for each contig
utils.syscall('samtools faidx ' + ref_fasta)
utils.syscall('scaff_test_make_unique_tags.py ' + contigs_outfile + ' ' + ref_fasta + ' ' + contigs_outfile + '.tag')

