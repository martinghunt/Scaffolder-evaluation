Scaffolder Evaluation
=====================

Wrapper scripts to run genome assembly scaffolding tools and
scripts to analyse the output for accuracy.

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Wrapper scripts
===============

The `Wrapper-scripts/` directory contains a script for each of the tools
listed below. Note that none of these tools are included in the code. Each
one must be installed separately and the scripts assume that the relevant
programs are in your $PATH.

The scaffolding tools with links to their publications and where to download
the code are as follows.

 * ABySS [code] [ABySS code], [paper] [ABySS paper]
 * Bambus2 (part of the AMOS package) [code] [Bambus2 code], [paper] [Bambus2 paper]
 * GRASS [code] [GRASS code], [paper] [GRASS paper]
 * MIP [code] [MIP code], [paper] [MIP paper]
 * Opera [code] [Opera code], [paper] [Opera paper]
 * SCARPA [code] [SCARPA code], [paper] [SCARPA paper]
 * SGA [code] [SGA code], [paper] [SGA paper]
 * SOAPdenovo2 [code] [SOAPdenovo2 code], [paper] [SOAPdenovo2 paper]
 * SOPRA [code] [SOPRA code], [paper] [SOPRA paper]
 * SSPACE [code] [SSPACE code], [paper] [SSPACE paper]

Run each script without any options to get the usage instructions.

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Analysis scripts
================

The `Analysis-scripts/` directory contains scripts to analyse the accuracy
of a genome scaffolder.

Installation
------------

Prerequisites:

 * [Python3] [python]
 * The Python package [Fastaq] [Fastaq]
 * The Python package [Pysam] [Pysam]
 * [MUMmer] [MUMmer]
 * [Bowtie2] [Bowtie2]
 * [SAMtools] [SAMtools]
 * [Graphviz] [Graphviz]
 * [R] [R]
 * [blastall], from NCBI's legacy BLAST executables

Take a copy of the scripts in the `Analysis-scripts/` directory and make sure
they are all in your $PATH to run the analysis.


Protocol to analyse scaffolds
-----------------------------

The starting point is the following data:

 * `assembly_contigs.fasta` - a FASTA file of the contigs that are to be scaffolded
 * `reference.fasta` - a FASTA file of the reference genome
 * `reads_1.fastq`, `reads_2.fastq` - FASTQ files of the reads that will be used to scaffold the contigs.

First, make artificial contigs from the assembly contigs.

    scaff_test_ctg_to_perfect_ctg.py assembly_contigs.fa reference.fa artificial_contigs

Now make sequence tags and only keep the artificial contigs that could be successfully tagged.

    scaff_test_make_unique_tags.py artificial_contigs.fa reference.fa artificial_contigs.tag

Run your scaffolding tool of choice on the tagged artificial contigs
`artificial_contigs.tag.uniquely-tagged.contigs.fa`, to make a file of scaffolds called `scaffolds.fa`.
Then check the output of the scaffolder:

    samtools faidx reference.fa
    scaff_test_check_using_tags.py <insert size> artificial_contigs.tag scaffolds.fa reference.fa.fai check_scaffolds

where `<insert size>` should be set to the insert size of your paired reads.
This makes numerous files called `check_scaffolds.*`. The important one is

    check_scaffolds.log

which has the counts of tag types at its end. The tag types are bitwise flags
(in a similar style to the samtools flag). The meanings are:

0 -- correct pair of tags.

1 -- tags originate from same reference sequence, but their orientation in the scaffolds is incorrect.

2 -- tags originate from different reference sequences.

4 -- tags originate from the same reference sequence but are the wrong distance apart.

8 -- tags originate from the same reference sequence but are not in the correct order.

For example, a flag of 5 means that 4+1 happened, i.e. a pair of tags that
originated from the same reference sequence, but their orientation and order
were incorrect.  Similarly, 12=8+4 means that two tags were from the same
reference sequence, but were the wrong distance apart and in the wrong order.
Finally, "lost" means that the tag was not found in the output of the
scaffolder and "skipped" means a pair of correct tags, but they skipped
over another tag.

Test cases
----------

Make test case contigs and reads for scaffolding using

    scaff_test_make_test_cases.sh <seed>

where `<seed>` is one of 1,2,...,5 to make one of the 5 repetitions of the data.
The important files for test case n are:

 * `test.n.ref.fa` - the 'reference' genome
 * `test.n.contigs.fa.tag.uniquely-tagged.contigs.fa` - the tagged contigs for scaffolding
 * `test.n.reads_1.fq`, `test.n.reads_2.fq` - the reads to use for scaffolding

Run the scaffolder to make a file of scaffolds called `scaffolds.fa`.  The
minimum, average and maximum insert sizes are 350, 500 and 650. Use 0.4
for the 'standard deviation' of insert size with SSPACE.  Check the
accuracy of the scaffolds with the following.

    scaff_test_check_using_tags.py 500 test.n.contigs.fa.tag scaffolds.fa test.n.ref.fa.fai check_scaffolds
    scaff_test_check_using_tags_get_contig_layout.py check_scaffolds check_scaffolds.graph

The second script made a using Graphviz, so you can easily see what the scaffolder
did.

Note that if you use SGA on these test cases, it will break unless you
specify `-b N` when running `sga-astat.py`, where `N` is the number of contigs
to be scaffolded. In the wrapper script, this can be set using `-astat_ops
"-b N"`.

Simulated data
--------------

Generate simulated contigs from the S. aureus genome from the [GAGE] dataset,
excluding its plasmid sequences:

    wget http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/genome.fasta
    samtools faidx genome.fasta
    samtools faidx genome.fasta NC_010079 > genome_no_plasmid.fa

Then generate 3kb and 10kb contigs from this sequence:

    scaff_test_make_simulated_contigs_from_genome.py genome_no_plasmid.fa 300 10000 genome_no_plasmid.fake_contigs.10kb_300gap.fa
    scaff_test_make_simulated_contigs_from_genome.py genome_no_plasmid.fa 50 3000 genome_no_plasmid.fake_contigs.3kb_50gap.fa

Generate reads with scripts from the Fastaq package:

    fastaq_to_perfect_reads genome_no_plasmid.fa fake_reads.3kb.fq 3000 200 20 76
    fastaq_to_perfect_reads genome_no_plasmid.fa fake_reads.500bp.fq 500 30 20 76
    fastaq_deinterleave fake_reads.3kb.fq fake_reads.3kb_1.fq fake_reads.3kb_2.fq
    fastaq_deinterleave fake_reads.500bp.fq fake_reads.500bp_2.fq fake_reads.500bp_2.fq

We now have our reference sequence, contigs and paired reads. The analysis
proceeds as above in the "Protocol to analyse scaffolds" section.


  [ABySS code]: http://www.bcgsc.ca/platform/bioinfo/software/abyss
  [ABySS paper]: http://genome.cshlp.org/content/19/6/1117
  [Bambus2 code]: http://sourceforge.net/projects/amos/
  [Bambus2 paper]: http://bioinformatics.oxfordjournals.org/content/27/21/2964.long
  [GRASS code]: https://code.google.com/p/tud-scaffolding/
  [GRASS paper]: http://bioinformatics.oxfordjournals.org/content/28/11/1429
  [MIP code]: http://www.cs.helsinki.fi/u/lmsalmel/mip-scaffolder/
  [MIP paper]: http://bioinformatics.oxfordjournals.org/content/27/23/3259
  [Opera code]: http://sourceforge.net/projects/operasf/files/version%201.0/
  [Opera paper]: http://online.liebertpub.com/doi/abs/10.1089/cmb.2011.0170
  [SCARPA code]: http://compbio.cs.toronto.edu/hapsembler/scarpa.html
  [SCARPA paper]: http://bioinformatics.oxfordjournals.org/content/29/4/428
  [SGA code]: https://github.com/jts/sga
  [SGA paper]: http://genome.cshlp.org/content/22/3/549
  [SOAPdenovo2 code]: http://soap.genomics.org.cn/soapdenovo.html
  [SOAPdenovo2 paper]: http://www.gigasciencejournal.com/content/1/1/18
  [SOPRA code]: http://www.physics.rutgers.edu/~anirvans/SOPRA/
  [SOPRA paper]: http://www.biomedcentral.com/1471-2105/11/345/
  [SSPACE code]: http://www.baseclear.com/landingpages/basetools-a-wide-range-of-bioinformatics-solutions/
  [SSPACE paper]: http://bioinformatics.oxfordjournals.org/content/27/4/578
  [Python]: http://www.python.org/
  [Fastaq]: https://github.com/sanger-pathogens/Fastaq
  [Pysam]: http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/samtools/api.html
  [MUMmer]: http://mummer.sourceforge.net/
  [Bowtie2]: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
  [blastall]: http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
  [Graphviz]: http://www.graphviz.org/
  [SAMtools]: http://samtools.sourceforge.net/
  [R]: http://www.r-project.org/
  [GAGE]: http://gage.cbcb.umd.edu/
