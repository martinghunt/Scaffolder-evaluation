#!/usr/bin/env python3

import argparse
import os
from pyfastaq import utils
import pysam

parser = argparse.ArgumentParser(
    description = 'Works out the layout of the contigs within scaffolds, using the file *.tags_and_sam.gz file made by the script scaffold_test_check_using_tags.py',
    usage = '%(prog)s [options] <inprefix> <outprefix>')
parser.add_argument('inprefix', help='Prefix of input files. Use the outprefix when scaff_test_check_using_tags.py was run')
parser.add_argument('outprefix', help='Prefix of output files')
options = parser.parse_args()

# load flags into memory
f = utils.open_file_read(options.inprefix + '.tags.gz')
flags = f.readlines()
utils.close(f)
flags = [int(x) for x in flags]

# load sam records into memory
sam_reader = pysam.Samfile(options.inprefix + '.tag_pairs.bam', 'rb')
lines = []
for sam in sam_reader:
    lines.append(sam)

nodes = {}

# loop over flag pairs, making graph nodes and adjacency lists
for i in range(0, len(lines), 2):
    flag = flags[int(i/2)]
    sam1 = lines[i]
    sam2 = lines[i+1]

    if sam1.qname not in nodes:
        nodes[sam1.qname] = set()


    if sam2.qname not in nodes:
        nodes[sam2.qname] = set()

    if sam1.tid == sam2.tid:

        if sam1.pos < sam2.pos:
            if (not sam1.is_reverse) and (not sam2.is_reverse):
                nodes[sam1.qname].add(sam2.qname)
            elif sam1.is_reverse and sam2.is_reverse:
                nodes[sam2.qname].add(sam1.qname)
            else:
                nodes[sam1.qname].add(sam2.qname)
                nodes[sam2.qname].add(sam1.qname)

        else:
            if (not sam1.is_reverse) and (not sam2.is_reverse):
                nodes[sam2.qname].add(sam1.qname)
            elif sam1.is_reverse and sam2.is_reverse:
                nodes[sam1.qname].add(sam2.qname)
            else:
                nodes[sam1.qname].add(sam2.qname)
                nodes[sam2.qname].add(sam1.qname)

# make pdf of the graph using graphviz
cmd = 'echo "digraph G {'
first  = True

for node, l in sorted(nodes.items()):
    if first:
        first = False
    else:
        cmd += ';'

    if len(l):
        cmd += ';'.join([node + '->' + x for x in l])
    else:
        cmd += node

cmd += '}"  | dot -Tpdf > ' + options.outprefix + '.pdf'

make_graph = options.outprefix + '.make_graph.sh'
f = utils.open_file_write(make_graph)
print(cmd, file=f)
utils.close(f)
utils.syscall('bash ' + make_graph)
