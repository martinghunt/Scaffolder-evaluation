#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Spec;

my %options = (
    'astat_ops' => '',
);
my $options_ok = GetOptions(\%options,
    'astat_ops=s',
    'noclean',
);

unless ($#ARGV >= 3 or !($options_ok)) {
    print STDERR "usage: $0 [options] <contigs.fa> <in.raw.bam> <in.sorted.bam> <outdir> [in.raw.bam2, in.raw.bam3, ...]

where the optional [in.raw.bam2,....] is list of second, third ..etc mapped libraries.
All libraries must be innies (i.e. FR (i.e. --> <--))
The first library given should be the highest coverage.

Script assumes that these are in your path:
sga, sga-bam2de.pl, sga-astat.py

Options:

-astat_ops \"STRING\"
    Options to pass to the sga-astat.py script. Put them in quotes, so they
    make it through to the call in one piece.
    (e.g. -astat_ops \"-b 10 -g 1000000\")

-noclean
    Don't clean up files. Default is to delete most files when finished
";
    exit(1);
}

my @raw_bams;
my @de_files;

my $contigs = File::Spec->rel2abs($ARGV[0]);
push @raw_bams, File::Spec->rel2abs($ARGV[1]);
my $sorted_bam = File::Spec->rel2abs($ARGV[2]);
my $outdir = $ARGV[3];

foreach (@ARGV[4..$#ARGV]) {
    push @raw_bams, File::Spec->rel2abs($_);
}

my $run_script = "run.sh";

mkdir $outdir or die $!;
chdir $outdir or die $!;

# write and run bash script to run all the stages
open F, ">$run_script" or die $!;

print F "set -e\n";

for my $i (0..$#raw_bams) {
    my $de_file = "bam2de.out.$i";
    push @de_files, "$de_file.de";
    print F "sga-bam2de.pl --mina 1 --prefix $de_file $raw_bams[$i]\n";
}

my $pe_string = "--pe " . join(' --pe ', @de_files);
print F "sga-astat.py $options{astat_ops} $sorted_bam > contigs.astat
sga scaffold -a contigs.astat $pe_string -o sga.scaffold.out $contigs
sga scaffold2fasta --write-unplaced -o scaffolds.fa --use-overlap -f $contigs sga.scaffold.out
";

unless ($options{noclean}) {
    print F "rm bam2de.out.* contigs.astat\n";
}

close F;

exec "bash $run_script" or die $!;
