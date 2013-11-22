#!/usr/bin/env perl

use strict;
use warnings;
use POSIX;
use Getopt::Long;
use File::Basename;
use File::Spec;
use File::Spec::Link;
use File::Which;

my %options;
my $options_ok = GetOptions(\%options,
    'noclean'
);

if (!($#ARGV == 3 or $#ARGV == 5) or !($options_ok)) {
    print STDERR "usage: $0 [options] <contigs.fa> <output_directory> <reads_1_fwd.fq> <reads_1_rev.fq> [<reads_2_fwd.fq> <reads_2_rev.fq>]

where extra in [square brackets] can be used with a second library.

Notes:

1. This script assumes ABySS scripts are in your path. Specifically:
   abyss-map, AdjList, abyss-scaffold, PathConsensus, MergeContigs

2. This script was tested on Abyss version 1.3.6

Options:

-noclean
    Use this to not clean up files.
    Default is to delete most files when finished
";

    exit(1);
}

my $contigs = File::Spec->rel2abs($ARGV[0]);
my $outdir = $ARGV[1];
my @reads_fwd = (File::Spec->rel2abs($ARGV[2]));
my @reads_rev = (File::Spec->rel2abs($ARGV[3]));

if ($#ARGV == 5){
    push @reads_fwd, File::Spec->rel2abs($ARGV[4]);
    push @reads_rev, File::Spec->rel2abs($ARGV[5]);
}


mkdir $outdir or die $!;
chdir $outdir or die $!;

open F, ">run.sh" or die $!;

print F "set -e
fastaq_to_fasta -l0 $contigs - | awk '\$1~/^>/ {i++; \$1=\">\"i} {print}' > contigs.fa
";

my $dot_files = "";
my $kmer = 1000000;

foreach my $i (0..$#reads_fwd) {
    my $readlength = get_read_length($reads_fwd[$i]);
    $kmer = $readlength < $kmer ? $readlength : $kmer;
}

$kmer = floor($kmer * 0.8);
$kmer = 64 < $kmer ? 64 : $kmer; # 64 is max kmer size (if abyss compiled with the deafults)

foreach my $i (0..$#reads_fwd) {
    print F "abyss-map -l$kmer $reads_fwd[$i] $reads_rev[$i] contigs.fa | abyss-fixmate -h lib.$i.hist | sort -snk3 -k4 | DistanceEst --dot -k$kmer -s200 -n10 -o lib.$i.dist.dot lib.$i.hist\n";
    $dot_files .= " lib.$i.dist.dot";
}

print F "abyss-scaffold -k$kmer -s200 -n10 -g x-6.path.dot contigs.fa $dot_files > x-6.path
PathConsensus -k$kmer -p0.9 -s x-7.fa -g x-7.adj -o x-7.path contigs.fa contigs.fa x-6.path
cat contigs.fa x-7.fa | MergeContigs -k$kmer -o scaffolds.fa - x-7.adj x-7.path
";

unless ($options{noclean}) {
    print F "rm lib.* contigs.fa x-*\n";
}

close F;

exec "bash run.sh" or die $!;

sub get_read_length {
    my $filename = shift;
    open FIN, $filename or die $!;
    <FIN>;
    my $seq = <FIN>;
    chomp $seq;
    close FIN;
    return length $seq;
}
