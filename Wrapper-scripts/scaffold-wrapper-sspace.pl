#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use File::Spec::Link;
use File::Which;

my %options = (
    'sspace_ops' => '',
    'sspace' => ''
);

my $options_ok = GetOptions(\%options,
    'sspace_ops=s',
    'noclean',
    'sspace=s',
);

if (!($#ARGV == 5 or $#ARGV == 9) or !($options_ok)) {
    print STDERR "usage: $0 [options] <contigs.fa> <reads_1.fq> <reads_2.fq> <output dir> <insert size> <insert sd> [<reads_1.fq> <reads_2.fq> >insert size> <insert sd>]

where extra in [square brackets] can be used with a second library.

Options:

-sspace /path/to/SSPACE_Basic_v2.0.pl
    Script assumes that SSPACE_Basic_v2.0.pl is in your path unless
    you use this option to specify the full path to the SSPACE script

-sspace_ops \"STRING\"
    Put SSPACE options in quotes. Unless this option is used,
    SSPACE's defaults will be used.

-noclean
    Don't clean up files. Default is to delete most files when finished
";
    exit(1);
}


my @reads_1;
my @reads_2;
my @insert_sizes;
my @insert_sds;


my $contigs = File::Spec->rel2abs($ARGV[0]);
push @reads_1, File::Spec->rel2abs($ARGV[1]);
push @reads_2, File::Spec->rel2abs($ARGV[2]);
my $outdir = $ARGV[3];
push @insert_sizes, $ARGV[4];
push @insert_sds, $ARGV[5];


# can't run SSPACE if it's a symbolic link. Need to use the actual file, otherwise
# it crashes when trying to use some modules that come with SSPACE.
my $sspace_exe = $options{sspace} ? $options{sspace} : which('SSPACE_Basic_v2.0.pl');
$sspace_exe = File::Spec::Link->resolve($sspace_exe);

my $run_script = 'run.sh';

if ($#ARGV == 9) {
    push @reads_1, File::Spec->rel2abs($ARGV[6]);
    push @reads_2, File::Spec->rel2abs($ARGV[7]);
    push @insert_sizes, $ARGV[8];
    push @insert_sds, $ARGV[9];
}


mkdir $outdir or die $!;
chdir $outdir or die $!;

open F, ">lib" or die $!;
for my $i (0..$#reads_1) {
    print F "LIB $reads_1[$i] $reads_2[$i] $insert_sizes[$i] $insert_sds[$i] FR\n";
}
close F;

open F, ">$run_script" or die $!;
print F "set -e
perl $sspace_exe $options{sspace_ops} -l lib -s $contigs
ln -s standard_output.final.scaffolds.fasta scaffolds.fa\n";

unless ($options{noclean}) {
    print F "rm -fr reads/ bowtieoutput/ intermediate_results/ pairinfo/\n";
}

close F;

exec "bash $run_script" or die $!;
