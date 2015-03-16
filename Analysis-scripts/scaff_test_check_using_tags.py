#!/usr/bin/env python3

import argparse
import copy
import os
import sys
from operator import attrgetter

from pyfastaq import *
import pysam

parser = argparse.ArgumentParser(
    description = 'Checks the accuracy of scaffolding, using a file of unique tags for each contig that went into the scaffolder',
    usage = '%(prog)s [options] <insert size> <tags_files_prefix> <scaffolds.fa> <reference.fa.fai> <outprefix>')
parser.add_argument('-c','--circular', action='append', help='Use this if a sequence is circular. This option can be used multiple times (once for each sequence that is cirulcar)', metavar='SEQUENCE_NAME')
parser.add_argument('--all_circular', action='store_true', help='Use this if all sequences are circular (overrides --circular)')
parser.add_argument('insert_size', type=float, help='To count as within the correct distance, this is the cutoff of difference between the expected distance and observed distance of two tags. So set it to the insert size to allow the distance error to be one insert size')
parser.add_argument('tags_files_prefix', help='Prefix of files made by scaffold_test_make_unique_tags.py')
parser.add_argument('scaffolds_fa', help='Name of scaffolds fasta file', metavar='scaffolds.fa')
parser.add_argument('reference_fai', help='Name of reference .fai file', metavar='reference.fa.fai')
parser.add_argument('outprefix', help='Prefix of output files')
options = parser.parse_args()

WRONG_ORIENTATION = 1
WRONG_CHR = 2
WRONG_DIST = 4
WRONG_ORDER = 8
NA = 16
DUP = 32

class Tag:
    def __init__(self, id, chr, position, length):
        self.id = id
        self.chr = chr
        self.position = int(position)
        self.length = int(length)
        self.ordered_index = -1 # will be putting these into lists, ordered
                                # by position. This will be the index in that
                                # list after sorting it.

    def __lt__(self, other):
        assert self.chr == other.chr  # don't intend to compare tags from different chromosomes
        return self.position < other.position

    def __str__(self):
        return '\t'.join([str(x) for x in [self.id, self.chr, self.position, self.length, self.ordered_index]])

def get_tag_from_sam(tag, tag_list):
    for (key, val) in tag_list:
       if key == tag:
           return val

    return None

reference_lengths = {}
tasks.lengths_from_fai(options.reference_fai, reference_lengths)

tags = {} # id -> tag
tags_by_chr = {}
tags_tsv_file = options.tags_files_prefix + '.tags.tsv'
tags_fa_file = options.tags_files_prefix + '.uniquely-tagged.tags.fa'

if options.circular:
    circular_seqs = set(options.circular)
else:
    circular_seqs = set()

# load tags from file
f = utils.open_file_read(tags_tsv_file)
for line in f:
    a = line.rstrip().split('\t')
    assert ' ' not in a[-1]
    (chr, pos, strand) = a[-1].split(':::')
    assert strand == '+'
    tag  = Tag(a[0], chr, pos, len(a[3]))
    assert tag.id not in tags
    tags[tag.id] = tag

    if tag.chr not in tags_by_chr:
        tags_by_chr[tag.chr] = []

    tags_by_chr[tag.chr].append(tag)

utils.close(f)

print('Got', len(tags), 'tags', file=sys.stderr)

# sort the tags into reference order for each chromosome
for chr, l in tags_by_chr.items():
    l.sort()
    for i in range(len(l)):
        l[i].ordered_index = i

# map the tags to the scaffolds
samfile = options.outprefix + '.map_tags.sam'
bamfile = options.outprefix + '.map_tags.bam'
sorted_bamfile = options.outprefix + '.map_tags.sorted.bam'
utils.syscall('bowtie2-build ' + options.scaffolds_fa + ' ' + options.scaffolds_fa)
utils.syscall('bowtie2 -f -x ' + options.scaffolds_fa + ' -U ' + tags_fa_file + ' -S ' + samfile)
utils.syscall('samtools view -T ' + options.scaffolds_fa + ' -bS ' + samfile + ' > ' + bamfile)
os.unlink(samfile)
utils.syscall('samtools sort ' + bamfile + ' ' + sorted_bamfile[0:-4])
#os.unlink(bamfile)

# Load the hits into memory
previous_sam = None
previous_tag = None
flag_counts = {k:0 for k in [0,1,2,4,5,8,12,16]}
tags_from_bam = set()
tag_distances = []
f_log = utils.open_file_write(options.outprefix + '.log')
f_tags = utils.open_file_write(options.outprefix + '.tags.gz')
skipped_tags = 0
sam_reader = pysam.Samfile(sorted_bamfile, "rb")
sam_writer = pysam.Samfile(options.outprefix + '.tag_pairs.bam', 'wb', header=sam_reader.header)

for current_sam in sam_reader.fetch(until_eof=True):
    if not current_sam.is_unmapped:
        tags_from_bam.add(current_sam.qname)
        as_tag = get_tag_from_sam('AS', current_sam.tags)
        xs_tag = get_tag_from_sam('XS', current_sam.tags)

        if as_tag != 0:
            print('Nonzero alignemnt score', current_sam, file=f_log)

        if xs_tag is not None and xs_tag >= as_tag:
            print('Non-unique best hit', current_sam, file=f_log)
    else:
        print('Unmapped', current_sam, file=f_log)


    if previous_sam is None:
        previous_sam = current_sam
        previous_tag = copy.copy(tags[current_sam.qname])
        continue

    flag = 0
    current_tag = copy.copy(tags[current_sam.qname])

    if current_sam.qname == previous_sam.qname: # should never happen, because of how mapper was run
        print('Warning: tag hit twice', current_sam, sep='\n\t', file=sys.stderr)
        flag += DUP
    elif (not current_sam.is_unmapped) and (not previous_sam.is_unmapped) and previous_sam.tid == current_sam.tid:
        if current_tag.chr != previous_tag.chr:  # if tags from different chromosomes
            flag += WRONG_CHR
        else:  # tags are from the same chromosome
            same_strand = previous_sam.is_reverse == current_sam.is_reverse
            mapped_distance = abs(current_sam.pos - previous_sam.pos)
            expected_distance = abs(current_tag.position - previous_tag.position)
            circular = False

            if options.all_circular or current_tag.chr in circular_seqs:
                previous_dist_to_right_end = reference_lengths[current_tag.chr] - previous_tag.position - previous_tag.length
                current_dist_to_right_end = reference_lengths[current_tag.chr] - current_tag.position - current_tag.length

                if previous_tag.position < current_tag.position:
                    expected_circular_distance = previous_tag.position + current_dist_to_right_end
                else:
                    expected_circular_distance = current_tag.position + previous_dist_to_right_end

                if expected_circular_distance < 2 * options.insert_size:
                    if expected_circular_distance < expected_distance:
                        circular = True

                    expected_distance = min(expected_circular_distance, expected_distance)

            tag_distance_error = abs(expected_distance - mapped_distance)

            if same_strand:
                if not previous_sam.is_reverse:
                    tag_order_ok = (not circular and previous_tag.ordered_index < current_tag.ordered_index) \
                                    or (circular and previous_tag.ordered_index > current_tag.ordered_index)
                elif previous_sam.is_reverse:
                    tag_order_ok = (not circular and previous_tag.ordered_index > current_tag.ordered_index) \
                                    or (circular and previous_tag.ordered_index < current_tag.ordered_index)

                if tag_order_ok:
                    tag_distances.append((expected_distance, mapped_distance))

                    if (not circular and abs(previous_tag.ordered_index - current_tag.ordered_index) != 1) \
                         or (circular and sorted([previous_tag.ordered_index, current_tag.ordered_index]) == [0, len(tags_by_chr[current_tag.chr])-1]):
                        skipped_tags += 1
                else:
                    flag += WRONG_ORDER
            else:
                flag += WRONG_ORIENTATION

            if tag_distance_error > options.insert_size:
                flag += WRONG_DIST
    else:
        flag += NA

    flag_counts[flag] = flag_counts.get(flag, 0) + 1
    print(flag, file=f_tags)
    sam_writer.write(previous_sam)
    sam_writer.write(current_sam)

    previous_sam = current_sam
    previous_tag = copy.copy(tags[current_sam.qname])

utils.close(f_tags)
sam_writer.close()

assert set(tags.keys()).issuperset(tags_from_bam)
lost_tags = set(tags.keys()).difference(tags_from_bam)
for t in lost_tags:
    print('LOST', t, file=f_log)

for n in sorted(flag_counts.keys()):
    print(n, flag_counts[n], sep='\t', file=f_log)

print('lost', len(lost_tags), sep='\t', file=f_log)
print('skipped', skipped_tags, sep='\t', file=f_log)
utils.close(f_log)

f = utils.open_file_write(options.outprefix + '.tag_distances')
print('#expected\tobserved', file=f)
for t in tag_distances:
   print(t[0], t[1], sep='\t', file=f)
utils.close(f)

f = utils.open_file_write(options.outprefix + '.tag_distances.R')
print('a=read.csv(file="' + options.outprefix + '.tag_distances' + r'''", header=T, sep="\t")''', file=f)
print('pdf("' + options.outprefix + '.tag_distances.R.pdf")', file=f)
print('plot(a, xlim=c(0,10000), ylim=c(0,10000))', file=f)
print('dev.off()', file=f)
utils.close(f)
utils.syscall('R CMD BATCH ' + options.outprefix + '.tag_distances.R')
