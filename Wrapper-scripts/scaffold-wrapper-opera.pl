#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use File::Spec;
use File::Spec::Link;
use File::Which;
use Getopt::Long;

my %options = ('bin' => '');
my $options_ok = GetOptions(\%options,
    'bin=s',
    'noclean'
);

unless ($#ARGV == 5 or $#ARGV == 7 or !($options_ok)) {
    print STDERR "usage: $0 [options] <contigs.fa> <reads_1.fq> <reads_2.fq> <outdir> <bwa|bowtie> <kmer> [<reads_1.fq> <reads_2.fq>]

where the optional [<reads_1.fq> <reads_2.fq>] can be used with a second library.
kmer is supposed to be set to be the same as that used by the assembler when
making the contigs. If this is not known/relevant set it to zero.

options:

-bin /path/to/directory/
    Script assumes that opera is in your path, and also that the script
    preprocess_reads.pl that comes with OPERA is in the same directory as
    opera. You can use this option to give the directory in which these scripts
    exist.

-noclean
    Don't clean up files. Default is to delete most files when finished
";
    exit(1);
}

my @libs;
push @libs, {};

my $contigs = File::Spec->rel2abs($ARGV[0]);
$libs[0]{reads_1} = File::Spec->rel2abs($ARGV[1]);
$libs[0]{reads_2} = File::Spec->rel2abs($ARGV[2]);
$libs[0]{id} = 1;
my $outdir = $ARGV[3];
my $mapper = $ARGV[4];
my $kmer = $ARGV[5];

my $map_file = "preprocess.map";
my $config_file = "preprocess.config";
my $run_script = "run.sh";

my $opera_bin = $options{bin};

unless ($opera_bin) {
    my $opera = $options{bin} ? File::Spec->catfile($options{bin}, 'opera') : which('opera');
    $opera = File::Spec::Link->resolve($opera);
    (undef, $opera_bin) = fileparse($opera);
}

if ($#ARGV == 7){
    push @libs, {};
    $libs[$#libs]{reads_1} = File::Spec->rel2abs($ARGV[6]);
    $libs[$#libs]{reads_2} = File::Spec->rel2abs($ARGV[7]);
    $libs[$#libs]{id} = scalar (@libs);
}

mkdir $outdir or die $!;
chdir $outdir or die $!;

if ($kmer != 0) {
    # write config file
    open F,  ">config_file" or die $!;

    print F "# Output folder for final results
output_folder=OUT

# Contig file
contig_file=$contigs

# value of kmer used to produce contig file
kmer=$kmer

";

    for my $h (@libs) {
       print F "[LIB]
# Mapped read locations
map_file=preprocess.lib.$h->{id}.map

";
    }

    close F;
}

# write and run bash script to run all the stages
open F, ">$run_script" or die $!;

print F "set -e
ln -s $contigs contigs.fa\n";

my $mapping_files_list = "";
my @mapping_files;

for my $h (@libs) {
    print F "perl $opera_bin/preprocess_reads.pl contigs.fa $h->{reads_1} $h->{reads_2} preprocess.lib.$h->{id}.map $mapper\n";
    $mapping_files_list .= ",preprocess.lib.$h->{id}.map";
    push @mapping_files, "preprocess.lib.$h->{id}.map";
}

$mapping_files_list = substr($mapping_files_list, 1);

if ($kmer != 0) {
    print F "$opera_bin/opera config_file\n";
}
else {
    print F "$opera_bin/opera contigs.fa $mapping_files_list OUT\n";
}

unless ($options{noclean}) {
    print F "rm OUT/pairedEndReads_preprocess contigs.fa.*\n";
    print F "rm " . join(' ', @mapping_files) . "\n";

}

print F "ln -s OUT/scaffoldSeq.fasta scaffolds.fa\n";

close F;

exec "bash $run_script" or die $!;
