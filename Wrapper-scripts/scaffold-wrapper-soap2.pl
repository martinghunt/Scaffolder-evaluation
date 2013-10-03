#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Spec;

my %options;

my $options_ok = GetOptions(\%options,
    'noclean',
);

unless ($#ARGV == 4 or $#ARGV == 7 or !($options_ok)) {
    print STDERR "usage: $0 [options] <contigs.fa> <reads_1.fq> <reads_2.fq> <outdir> <insert size> [<reads_1.fq> <reads_2.fq> <insert size>]

Use optional extra in square brackets to use a second library.
Note: assumes all reads are the same length for each library.
Also assumes finalFusion (from the prepare module) and SOAPdenovo-63mer
are in your path.

Options:

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
my $outdir = $ARGV[3];
$libs[0]{insert} = $ARGV[4];
$libs[0]{id} = 1;

if ($#ARGV == 7){
    push @libs, {};
    $libs[$#libs]{reads_1} = File::Spec->rel2abs($ARGV[5]);
    $libs[$#libs]{reads_2} = File::Spec->rel2abs($ARGV[6]);
    $libs[$#libs]{insert} = $ARGV[7];
    $libs[$#libs]{id} = scalar (@libs);
}

my $soap2_prepare = 'finalFusion';
my $soap2 = 'SOAPdenovo-63mer';
my $config_file = 'config_file';
my $run_script = 'run.sh';
my $outprefix = 'out';


# get the max read length
my $max_read_length = 0;

for my $h (@libs) {
    my $f = $h->{reads_1};
    open F, "$f" or die "Error opening '$f'";
    <F>;
    my $seq = <F>;
    chomp $seq;
    my $read_length = length($seq);
    close F;

    if ($read_length > $max_read_length) {
        $max_read_length = $read_length;
    }
}

if ($max_read_length == 0) {
    print STDERR "Error getting max read length from fastq file(s)\n";
    exit(1);
}

mkdir $outdir or die $!;
chdir $outdir or die $!;

# write config file
open F,  ">$config_file" or die $!;

print F "
#maximal read length
max_rd_len=$max_read_length

";


for my $h (@libs) {
   print F "[LIB]
#average insert size
avg_ins=$h->{insert}
#if sequence needs to be reversed. 0=innies --> <--. 1=outies <-- -->
reverse_seq=0
#use for scaffolding only
asm_flags=2
#in which order the reads are used while scaffolding
rank=$h->{id}
#fastq files
q1=$h->{reads_1}
q2=$h->{reads_2}

";
}

close F;

# write and run bash script to run all the stages
open F, ">$run_script" or die $!;

print F "set -e
echo \"RUN PREPARE\"
$soap2_prepare -p 1 -D -c $contigs -g $outprefix
echo \"RUN MAP\"
$soap2 map -p 1 -s $config_file -g $outprefix
echo \"RUN SCAFF\"
touch $outprefix.Arc
$soap2 scaff -p 1 -g $outprefix
ln -s $outprefix.scafSeq scaffolds.fa
";

unless ($options{noclean}) {
    print F "ls out.* | grep -v '$outprefix.scafS' | xargs rm\n";
}

close F;

exec "bash $run_script" or die $!;

