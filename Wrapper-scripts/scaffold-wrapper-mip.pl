#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Spec;
use File::Spec::Link;
use File::Which;

my %options = (
    'bin' => '',
    'mapper_ops' => '',
);
my $options_ok = GetOptions(\%options,
    'bin=s',
    'mapper_ops=s',
    'noclean'
);

if (!($#ARGV == 9 or $#ARGV == 14) or !($options_ok)) {
    print STDERR "usage: $0 [options] <contigs.fa> <bwa|bowtie|bowtie2> <reads_1.fq> <reads_2.fq> <genome size upper bound> <min insert> <insert length> <max insert> <contig_cov_file> <output_directory> [<reads_1.fq> <reads_2.fq> <min insert> <insert length> <max insert>]

where extra in [square brackets] can be used with a second library.
Mapper must be one of bowtie,bowtie2,bwa.
To generate contig_cov_file use scaffold-wrapper-mip-make-cov-file.pl

Options:

-bin /path/to/directory/
    Script assumes that mip-scaffolder.pl is in your path, and also
    that the scripts filter-mappings.sh and merge.sh that come with MIP
    are in the same directory as mip-scaffolder.pl.
    You can use this option to give the directory in which all these
    scripts exist

-mapper_ops \"STRING\"
    Put mapper options in quotes to pass them into the mapping call.
    If not given, defaults used for bowtie2 and bwa, -m 1 is used for
    bowtie.

-noclean
    Don't clean up files. Default is to delete most files when finished
";

    exit(1);
}

my %lib1_data;
my %lib2_data;

my $contigs = File::Spec->rel2abs($ARGV[0]);
my $mapper = $ARGV[1];
$lib1_data{reads_1} = File::Spec->rel2abs($ARGV[2]);
$lib1_data{reads_2} = File::Spec->rel2abs($ARGV[3]);
my $genome_max = $ARGV[4];
$lib1_data{insert_min} = $ARGV[5];
$lib1_data{insert} = $ARGV[6];
$lib1_data{insert_max} = $ARGV[7];
$lib1_data{id} = 1;
my $contig_coverage = File::Spec->rel2abs($ARGV[8]);
my $outdir = $ARGV[9];

my $mip_scripts_dir = $options{bin};

unless ($mip_scripts_dir) {
    my $mip_scaffolder_pl = $options{bin} ? File::Spec->catfile($options{bin}, 'mip-scaffolder.pl') : which('mip-scaffolder.pl');
    $mip_scaffolder_pl = File::Spec::Link->resolve($mip_scaffolder_pl);
    (undef, $mip_scripts_dir) = fileparse($mip_scaffolder_pl);
}

if ($#ARGV == 14) {
    $lib2_data{reads_1} = File::Spec->rel2abs($ARGV[10]);
    $lib2_data{reads_2} = File::Spec->rel2abs($ARGV[11]);
    $lib2_data{insert_min} = $ARGV[12];
    $lib2_data{insert} = $ARGV[13];
    $lib2_data{insert_max} = $ARGV[14];
    $lib2_data{id} = 2;
}

my $bowtie = 'bowtie';
my $bowtie2 = 'bowtie2';
my $bwa = 'bwa';
$mapper eq 'bwa' or $mapper eq 'bowtie' or $mapper eq 'bowtie2' or die "Mapper must be bwa, bowtie or bowtie2";

mkdir $outdir or die $!;
chdir $outdir or die $!;

my $stage_config_string = "[STAGE]
# Maximum biconnected component size. (optional)
#maximum_biconnected_component=50
# Maximum allowed degree in scaffolding graph. (optional)
#maximum_degree=50
# Maximum coverage for nonrepetitive contig. (optional)
#maximum_coverage=20
# The maximum overlap between contigs that is allowed without checking for
# sequence similarity. By default this is set based on the variablility in
# insert size lengths of each library. (optional)
#maximum_overlap=100
# The minimum support for an edge. (optional)
#minimum_support=2
# Should edges with negative estimated distance be checked for sequence
# similarity or removed automatically? (optional)
#check_negative_edges=1
# The maximum allowed error level when checking for sequence similarity
# (optional)
#alignment_error=0.1
";

# write config file
open F, ">config_file" or die $!;

print F "# Upper bound for genome length (required)
genome_length=$genome_max

#parameter specifications for the first stage
$stage_config_string

# library specification for the first stage
";

print F lib_data_to_library_config_string(\%lib1_data);

if (%lib2_data) {
    print F "
#parameter specifications for the second stage
$stage_config_string

# library specification for the second stage
";

    print F lib_data_to_library_config_string(\%lib2_data);
}

close F;

# write bash script to run all the stages
my $index_cmd;

if ($mapper eq 'bwa'){
    $index_cmd = "$bwa index contigs.fa contigs.fa";
}
elsif ($mapper eq 'bowtie'){
    $index_cmd = "$bowtie-build contigs.fa contigs.fa";
    if ($options{mapper_ops} eq '') {
        $options{mapper_ops} = '-m 1';
    }
}
elsif ($mapper eq 'bowtie2'){
    $index_cmd = "$bowtie2-build contigs.fa contigs.fa";
}


open F, ">run.sh" or die $!;

print F "set -e
ln -s $contigs contigs.fa
ln -s $contig_coverage contig_coverage
$index_cmd\n";
print F lib_data_to_mapping_string(\%lib1_data, "contigs.fa", $mip_scripts_dir, $mapper, $options{mapper_ops});

if (%lib2_data) {
    print F lib_data_to_mapping_string(\%lib2_data, "contigs.fa", $mip_scripts_dir, $mapper, $options{mapper_ops});
}
print F "$mip_scripts_dir/mip-scaffolder.pl config_file contigs.fa contig_coverage MIP-out/\n";
print F "ln -s MIP-out/scaffolds.fasta scaffolds.fa\n";

unless ($options{noclean}) {
    print F "rm lib.* MIP-out/lib.* MIP-out/*.current contigs.fa.*\n";
}

close F;

exec "bash run.sh" or die $!;


sub lib_data_to_library_config_string {
    my $data = shift; # ref to hash of data

    return "[LIBRARY]
# File in SAM format containing mappings for the mate pair reads
# to the contigs
mappings=lib.$data->{id}.filtered-mappings
# Orientation of the mate pairs (in current version must be SOLID)
orientation=SOLID
# Insert length
insert_length=$data->{insert}
# Minimum insert length
min_insert_length=$data->{insert_min}
# Maximum insert length
max_insert_length=$data->{insert_max}
";
}


sub lib_data_to_mapping_string {
    my $data = shift; # ref to hash of data
    my $contigs_fname = shift;
    my $scripts_dir = shift;
    my $mapper = shift;
    my $mapper_ops = shift;

    #my $s = "awk '{print substr(\$1,1,length(\$1)-2)\"_R3\"; for (i=0;i<3;i++){getline;print}}' $data->{reads_1} > lib.$data->{id}.reads_R3.fq
#fastn_revcomp $data->{reads_2} | awk '{print substr(\$1,1,length(\$1)-2)\"_F3\"; for (i=0;i<3;i++){getline;print}}' > lib.$data->{id}.reads_F3.fq\n";
    my $s = "scaffold-wrapper-mip-convert-fastqs.pl $data->{reads_1} $data->{reads_2} lib.$data->{id}.reads\n";

    if ($mapper eq 'bwa') {
        $s .= "$bwa aln $mapper_ops $contigs_fname lib.$data->{id}.reads_R3.fastq > lib.$data->{id}.bowtie.reads_R3.sai
$bwa samse -n 1 $contigs_fname lib.$data->{id}.bowtie.reads_R3.sai lib.$data->{id}.reads_R3.fastq > lib.$data->{id}.bowtie.reads_R3.sam
$bwa aln $mapper_ops $contigs_fname lib.$data->{id}.reads_F3.fastq > lib.$data->{id}.bowtie.reads_F3.sai
$bwa samse -n 1 $contigs_fname lib.$data->{id}.bowtie.reads_F3.sai lib.$data->{id}.reads_F3.fastq > lib.$data->{id}.bowtie.reads_F3.sam\n";
    }
    elsif ($mapper eq 'bowtie') {
        $s .= "$bowtie $mapper_ops --sam $contigs_fname lib.$data->{id}.reads_R3.fastq > lib.$data->{id}.bowtie.reads_R3.sam
$bowtie $mapper_ops --sam $contigs_fname lib.$data->{id}.reads_F3.fastq > lib.$data->{id}.bowtie.reads_F3.sam\n";
    }
    elsif ($mapper eq 'bowtie2') {
        $s .= "$bowtie2 $mapper_ops $contigs_fname -U lib.$data->{id}.reads_R3.fastq -S lib.$data->{id}.bowtie.reads_R3.sam
$bowtie2 $mapper_ops $contigs_fname -U lib.$data->{id}.reads_F3.fastq -S lib.$data->{id}.bowtie.reads_F3.sam\n";
    }

    $s .= "
$scripts_dir/merge.sh lib.$data->{id}.bowtie.reads_F3.sam lib.$data->{id}.bowtie.reads_R3.sam lib.$data->{id}.merged-mappings
$scripts_dir/filter-mappings.sh lib.$data->{id}.merged-mappings.sorted1 lib.$data->{id}.merged-mappings.sorted2 lib.$data->{id}.filtered-mappings

";
}
