#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Spec;
use File::Spec::Link;
use File::Which;
use Getopt::Long;

my %options = (
    'mapper_ops' => '',
    'bin' => ''
);

my $options_ok = GetOptions(\%options,
    'bin=s',
    'mapper_ops=s',
    'noclean'
);

if (!($#ARGV == 4 or $#ARGV == 6) or !($options_ok)) {
    print STDERR "usage: $0 [options] <contigs.fa> <reads_interleaved.fq> <output dir> <insert size> <mapper> [<reads_interleaved2.fq> <insert size2>]

where:
 - mapper must be one of bowtie,bowtie2,bwa. It is assumed that the mapper
   is in your path.
 - extra in [square brackets] can be used for a second library.

Options:

-bin /path/to/directory/
    Script assumes that scarpa is in your path, and also that scarpa_parser
    and scarpa_process that come with SCARPA are in the same directory as
    scarpa (after following any symlinks).
    You can use this option to give the directory in which these exist, in
    case scarpa is not in your path

-mapper_ops \"STRING\"
    Put mapper options in quotes to pass them into the mapping call.
    If not given, defaults used for bowtie2 and bwa, -m 1 is used for
    bowtie.

-noclean
    Don't clean up files. Default is to delete most files when finished
";
    exit(1);
}


my $contigs = File::Spec->rel2abs($ARGV[0]);
my $reads = File::Spec->rel2abs($ARGV[1]);
my $outdir = $ARGV[2];
my $insert_size = $ARGV[3];
my $mapper = $ARGV[4];
my $reads2 = '';
my $insert_size2 = 0;


if ($#ARGV == 6) {
    $reads2 = File::Spec->rel2abs($ARGV[5]);
    $insert_size2 = $ARGV[6];
}

my $scarpa_bin = $options{bin};

unless ($scarpa_bin) {
    my $scarpa = $options{bin} ? File::Spec->catfile($options{bin}, 'scarpa') : which('scarpa');
    $scarpa = File::Spec::Link->resolve($scarpa);
    (undef, $scarpa_bin) = fileparse($scarpa);
}

my $scarpa = File::Spec->catfile($scarpa_bin,'scarpa');
my $scarpa_process =  File::Spec->catfile($scarpa_bin, 'scarpa_process');
my $scarpa_parser =  File::Spec->catfile($scarpa_bin, 'scarpa_parser');
my $bowtie = 'bowtie';
my $bowtie2 = 'bowtie2';
my $bwa = 'bwa';

my $index_cmd;
my $mapping_cmd;

my $scarpa_contigs = 'contigs.fa.scarpa.fa';
my $scarpa_reads = 'reads.fq.scarpa.fq';
my $scarpa_info = 'contigs.fa.scarpa.info';
my $scarpa_map_file = 'parser.map_file';

if ($mapper eq 'bowtie') {
    $index_cmd = "$bowtie-build $scarpa_contigs $scarpa_contigs";
    if ($options{mapper_ops}  eq '') {
        $options{mapper_ops} = '-m 1';
    }
    $mapping_cmd = "$bowtie $options{mapper_ops} --sam $scarpa_contigs $scarpa_reads";
}
elsif ($mapper eq 'bowtie2') {
    $index_cmd = "$bowtie2-build $scarpa_contigs $scarpa_contigs";
    $mapping_cmd = "$bowtie2 $options{mapper_ops} -x $scarpa_contigs -U $scarpa_reads";
}
elsif ($mapper eq 'bwa') {
    $index_cmd = "$bwa index $scarpa_contigs $scarpa_contigs";
    $mapping_cmd = "$bwa aln $options{mapper_ops} $scarpa_contigs $scarpa_reads > tmp.bwa.sai\n" .
                   "$bwa samse -n 1 $scarpa_contigs tmp.bwa.sai $scarpa_reads";
}
else {
    print STDERR "Error! Mapper must be one of: bowtie,bowtie2,bwa\n";
    exit(1);
}

$mapping_cmd .= " | $scarpa_parser > $scarpa_map_file";

my $run_script = 'run.sh';

mkdir $outdir or die $!;
chdir $outdir or die $!;

open F, ">$run_script" or die $!;
print F "set -e
ln -s $contigs contigs.fa
";

if ($reads2) {
    print F "ln -s $reads reads1.fq
ln -s $reads2 reads2.fq
$scarpa_process -c contigs.fa -f reads1.fq -i $insert_size -f reads2.fq -i $insert_size2
cat reads1.fq.scarpa.fq reads2.fq.scarpa.fq > $scarpa_reads
";
}
else {
    print F "ln -s $reads reads.fq
$scarpa_process -c contigs.fa -f reads.fq -i $insert_size
";
}

print F "
$index_cmd
$mapping_cmd
$scarpa -c $scarpa_contigs -i $scarpa_map_file -l $scarpa_info -o scaffolds.fa
";

unless ($options{noclean}) {
    print F "rm $scarpa_reads contigs.fa.* scaffolds.fa.scafftmp $scarpa_map_file\n";
}
close F;

exec "bash $run_script" or die $!;

