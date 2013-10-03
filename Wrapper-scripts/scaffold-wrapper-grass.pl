#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Spec;

my %options = (
    'bin' => ''
);

my $options_ok = GetOptions(\%options,
    'bin=s',
);


if ($#ARGV != 6 or !($options_ok)) {
    print STDERR "usage: $0 [options] <contigs.fa> <reads_1.fq> <reads_2.fq> <insert size> <insert std> <expected read coverage> <outdir>

Options:

-bin /path/to/directory/
    Script assumes that dataLinker and scaffoldOptimizer are in your path
    unless you use this option to give the directory in which both these
    programs exist.
";
    exit(1);
}


my $contigs = File::Spec->rel2abs($ARGV[0]);
my $reads1 = File::Spec->rel2abs($ARGV[1]);
my $reads2 = File::Spec->rel2abs($ARGV[2]);
my $insert_mean = $ARGV[3];
my $insert_sd = $ARGV[4];
my $exp_cov = $ARGV[5];
my $outdir = $ARGV[6];

my $dataLinker = $options{bin} ? File::Spec->catfile($options{bin}, 'dataLinker') : 'dataLinker';
my $scaffoldOptimizer = $options{bin} ? File::Spec->catfile($options{bin}, 'scaffoldOptimizer') : 'scaffoldOptimizer';

my $run_script = "run.sh";

mkdir $outdir or die $!;
chdir $outdir or die $!;

# write and run bash script to run all the stages
open F, ">$run_script" or die $!;

print F "set -e
$dataLinker -bwathreads 1 -illumina $reads1 $reads2 $insert_mean $insert_sd -readcoverage read.coverage -output datalinker.opt $contigs
$scaffoldOptimizer -repeat-coverage $exp_cov read.coverage -output scaffolds.fa datalinker.opt
";

close F;

exec "bash $run_script" or die $!;
