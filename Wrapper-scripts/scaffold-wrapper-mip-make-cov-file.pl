#!/usr/bin/env perl

use strict;
use warnings;

if ($#ARGV != 2) {
    print STDERR "usage: $0 <in.bam> <in.fasta.fai> <out.mip.coverage>

Given a BAM file, makes a 'coverage' file required as input to the
MIP scaffolder. Assumes samtools is in your path.
";
    exit(1);
}

my $bam_file = $ARGV[0];
my $fai_file = $ARGV[1];
my $outfile = $ARGV[2];

my %total_bases_mapped; # contig name -> total bases mapped
my %sequence_lengths;   # contig name -> length
my @sequence_order;     # sequences in order in the fai file

# get sequence lengths into a hash
open F, "$fai_file" or die "Error opening file '$fai_file'";
while (<F>) {
    chomp;
    my ($name, $length) = split;
    $sequence_lengths{$name} = $length;
    push @sequence_order, $name
}
close F;

# run mpileup and add up the total bases mapped
open F, "samtools mpileup -d 1000 -A $bam_file |" or die "Error opening file '$bam_file'";
while (<F>) {
    chomp;
    my ($name, undef, undef, $depth) = split;
    $total_bases_mapped{$name} += $depth;
}
close F or die $!;

# write the output file
open F, ">$outfile" or die "Error opening file '$outfile'";
foreach my $i (0..$#sequence_order) {
    my $seq = $sequence_order[$i];
    my $cov = $total_bases_mapped{$seq} ? ($total_bases_mapped{$seq} / $sequence_lengths{$seq}) : 0;
    print F ($i + 1) . "\t$seq\t$sequence_lengths{$seq}\t$cov\n";
}
close F;
