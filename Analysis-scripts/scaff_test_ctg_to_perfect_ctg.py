#!/usr/bin/env python3

import argparse
import os
import copy
import shutil
import sys
from pyfastaq import sequences, utils, intervals, tasks


# check required nucmer programs are in path
progs = ['nucmer', 'delta-filter', 'show-coords']
not_in_path = [p for p in progs if shutil.which(p) is None]
    
if len(not_in_path):
    print('Error! Need these programs to be in your path:', file=sys.stderr)
    print('\n'.join(not_in_path), file=sys.stderr)
    print('Cannot continue', file=sys.stderr)
    sys.exit(1)


def nucmer_file_reader(fname):
    f = utils.open_file_read(fname)
    in_header = True

    for line in f:
        if in_header:
            if line.startswith('['):
                in_header = False
            continue
        yield NucmerHit(line)

    utils.close(f)


class NucmerHit:
    def __init__(self, line):
        # [S1]  [E1]    [S2]    [E2]    [LEN 1] [LEN 2] [% IDY] [LEN R] [LEN Q] [FRM]   [TAGS]
        #1162    25768   24536   4   24607   24533   99.32   640851  24536   1   -1  MAL1    NODE_25757_length_24482_cov_18.920391   [CONTAINS]
        try:
            l = line.rstrip().split('\t')
            self.ref_start = int(l[0])
            self.ref_end = int(l[1])
            self.qry_start = int(l[2])
            self.qry_end = int(l[3])
            self.hit_length_ref = int(l[4])
            self.hit_length_qry = int(l[5])
            self.percent_identity = float(l[6])
            self.ref_length = int(l[7])
            self.qry_length = int(l[8])
            self.frame = int(l[9])
            self.strand = int(l[10])
            self.ref_name = l[11]
            self.qry_name = l[12]

            if len(l) == 14:
                self.tag = l[13][1:-1]
            else:
                self.tag = None
        except:
            print('Error reading this nucmer line:\n' + line, file=sys.stderr)


def update_perfect_contigs(nucmer_hit, ref_fasta, contigs):
    id = nucmer_hit.ref_name + ":" + str(nucmer_hit.ref_start) + '-' + str(nucmer_hit.ref_end)
    contig = sequences.Fasta('x', ref_fasta[nucmer_hit.ref_start-1:nucmer_hit.ref_end])
    contigs[(nucmer_hit.ref_name, nucmer_hit.ref_start, nucmer_hit.ref_end)] = contig


parser = argparse.ArgumentParser(
    description = 'Takes contigs and a reference sequence. Makes a new fasta file of the contigs, but they are now perfect sequences by using the reference instead',
    usage = '%(prog)s [options] <contigs.fa> <reference.fa> <outprefix>')
parser.add_argument('--min_seq_length', type=int, help='Minimum length of contig to output [%(default)s]', default=200)
parser.add_argument('--nucmer_options', help='Options when running nucmer [%(default)s]', default='')
parser.add_argument('contigs_fa', help='Name of contigs fasta file', metavar='contigs.fa')
parser.add_argument('ref_fa', help='Name of reference fasta file', metavar='reference.fa')
parser.add_argument('outprefix', help='Prefix of output files')
options = parser.parse_args()

ref_seqs = {}
tasks.file_to_dict(options.ref_fa, ref_seqs)

nucmer_out_prefix = options.outprefix + '.nucmer'
nucmer_out_delta = nucmer_out_prefix + '.delta'
nucmer_out_filter = nucmer_out_prefix + '.delta-filter'
nucmer_out_coords = nucmer_out_filter + '.coords'

# run nucmer of contigs vs ref
utils.syscall(' '.join(['nucmer', options.nucmer_options, '-p', nucmer_out_prefix, options.ref_fa, options.contigs_fa]))
utils.syscall(' '.join(['delta-filter', '-i 98 -l 180 -q', nucmer_out_delta, '>', nucmer_out_filter]))
utils.syscall(' '.join(['show-coords', '-dTlro', nucmer_out_filter, '>', nucmer_out_coords]))

# load hits into hash. key=ref_name, value=another hash with key=qry_name, value=list of hit positions in that ref seq
nucmer_hits = {}
contigs_to_print = {}

nucmer_reader = nucmer_file_reader(nucmer_out_coords)

for hit in nucmer_reader:
    if hit.ref_name not in nucmer_hits:
        nucmer_hits[hit.ref_name] = {}

    if hit.qry_name not in nucmer_hits[hit.ref_name]:
        nucmer_hits[hit.ref_name][hit.qry_name] = []

    nucmer_hits[hit.ref_name][hit.qry_name].append(intervals.Interval(min(hit.ref_start, hit.ref_end), max(hit.ref_start, hit.ref_end)))

# merge all the overalpping hits for each list of hits corresponding to one contig
for ref_name, d in nucmer_hits.items():
    for qry_name, hits in d.items():
        intervals.merge_overlapping_in_list(hits)

        for hit in hits:
            if hit.end - hit.start + 1 >= options.min_seq_length:
                if ref_name not in contigs_to_print:
                    contigs_to_print[ref_name] = []

                contigs_to_print[ref_name].append(copy.copy(hit))

# remove any contigs that are completely contained in another contig
for ref, l in contigs_to_print.items():
    intervals.remove_contained_in_list(l)

# print the final perfect contigs
f_out = utils.open_file_write(options.outprefix + '.fa')
counter = 1
last_id = None
for ref_name in sorted(contigs_to_print):
    counter = 1

    for interval in contigs_to_print[ref_name]:
        id = ':'.join([str(x) for x in [ref_name, counter, interval.start, interval.end]])
        print(sequences.Fasta(id, ref_seqs[ref_name][interval.start - 1: interval.end]), file=f_out)
        counter += 1

utils.close(f_out)

