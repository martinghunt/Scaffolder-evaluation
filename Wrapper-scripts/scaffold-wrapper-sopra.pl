#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use File::Spec;
use File::Spec::Link;
use File::Which;

my %options = ('mapper_ops' => '', 's_scaf_ops' => '', bin => '');
my $options_ok = GetOptions(\%options,
    'bin=s',
    'mapper_ops=s',
    's_scaf_ops=s',
    'noclean');

if (!($#ARGV == 6 or $#ARGV == 9) or !($options_ok)) {
    print STDERR "usage: $0 [options] <contigs.fa> <bwa|bowtie|bowtie2> <reads_1.fq> <reads_2.fq> <reads_per_split> <insert size> <outdir> [<reads_1.fq> <reads_2.fq> <insert size>]

Assumes reads are in FR orientation.
Use optional extras to use a second library

options:
-bin /path/to/directory/
    Script assumes that the SOPRA script s_scaf_v1.4.6.pl is in your path,
    and also that all other SOPRA scripts are in the same directory as
    s_scaf_v1.4.6.pl. You can use this option to give the directory in which
    these scripts exist.

-mapper_ops \"STRING\"
    Put mapper options in quotes to pass them into the mapping call.
    If not given, defaults used for bowtie2 and bwa, -m 1 -v 0 is used for
    bowtie.

-s_scaf_ops \"STRING\"
    Options to pass to the SOPRA scaffolding script s_scaf_vx.x.x.pl.

-noclean
    Don't clean up files. Default is to delete most files when finished
";
    exit(1);
}

my @libs;
push @libs, {};


my $contigs = File::Spec->rel2abs($ARGV[0]);
my $mapper = $ARGV[1];
$libs[0]{reads1} = File::Spec->rel2abs($ARGV[2]);
$libs[0]{reads2} = File::Spec->rel2abs($ARGV[3]);
$libs[0]{id} = 1;
my $reads_per_split = $ARGV[4];
$libs[0]{insert} = $ARGV[5];
my $outdir = $ARGV[6];

if ($#ARGV == 9) {
    push @libs, {};
    $libs[$#libs]{reads1} = File::Spec->rel2abs($ARGV[7]);
    $libs[$#libs]{reads2} = File::Spec->rel2abs($ARGV[8]);
    $libs[$#libs]{insert} = $ARGV[9];
    $libs[$#libs]{id} = scalar (@libs);
}


my $sopra_bin = $options{bin};

unless ($sopra_bin) {
    my $sopra_scaf = $options{bin} ? File::Spec->catfile($options{bin}, 's_scaf_v1.4.6.pl') : which('s_scaf_v1.4.6.pl');
    $sopra_scaf = File::Spec::Link->resolve($sopra_scaf);
    (undef, $sopra_bin) = fileparse($sopra_scaf);
}

my $run_script = "run.sh";
my $bowtie = 'bowtie';
my $bowtie2 = 'bowtie2';
my $bwa = 'bwa';

$mapper eq 'bwa' or $mapper eq 'bowtie' or $mapper eq 'bowtie2' or die "Mapper must be bwa, bowtie or bowtie2";


mkdir $outdir or die $!;
chdir $outdir or die $!;


my $split_reads_files = "";
my $sam_parsed_files = "";
my $split_reads_dir = "Split_reads";
mkdir $split_reads_dir or die $!;
my $last_file_index = 1;

# split the reads
for my $h (@libs) {
    my $file_count = 0;
    my $read_count = 0;
    open FWD, $h->{reads1} or die $!;
    open REV, $h->{reads2} or die $!;
    open OUT, ">$split_reads_dir/r.$last_file_index.fa" or die $!;


    while (my $fwd_line = <FWD>) {
        my $rev_line = <REV>;

        if (substr($fwd_line, 1, -3) ne substr($rev_line, 1, -3)) {
            die "Read name mismatch $fwd_line, $rev_line";
        }

        $fwd_line =~ s/@/>/;
        $rev_line =~ s/@/>/;

        if ($read_count >= $reads_per_split) {
            close OUT;
            $read_count = 0;
            open OUT, ">$split_reads_dir/r." . ($last_file_index + $file_count + 1) . ".fa";
            $file_count++;
        }

        my $read1 = <FWD>;
        my $read2 = <REV>;

        print OUT "$fwd_line$read1$rev_line$read2";
        $read_count+=2;

        <FWD>; <FWD>;
        <REV>; <REV>;
    }

    close OUT;
    close FWD;
    close REV;

    for my $i (($last_file_index)..($last_file_index + $file_count)) {
        $split_reads_files .= " $split_reads_dir/r.$i.fa";
        $sam_parsed_files .= " -parsed Split_map/$i.sam_parsed_p$i -d $h->{insert}";
    }

    $last_file_index += $file_count;
}

# make indexing and mapping commands. Depends on the mapper being used
my $index_cmd;
my $map_cmd;

if ($mapper eq 'bwa'){
    $index_cmd = "$bwa index OUT/contigs_sopra.fasta OUT/contigs_sopra.fasta";
    $map_cmd = "$bwa $options{mapper_ops} aln OUT/contigs_sopra.fasta OUT/r.\$i\\\_sopra.fasta > Split_map/\$i.aln.sai\n    " .
               "$bwa samse OUT/contigs_sopra.fasta Split_map/\$i.aln.sai OUT/r.\$i\\\_sopra.fasta > Split_map/\$i.sam";
}
elsif ($mapper eq 'bowtie'){
    $index_cmd = "$bowtie-build OUT/contigs_sopra.fasta OUT/contigs_sopra.fasta";
    if ($options{mapper_ops} eq '') {
        $options{mapper_ops} = '-m 1 -v 0';
    }
    $map_cmd = "$bowtie $options{mapper_ops} -f --sam OUT/contigs_sopra.fasta OUT/r.\$i\\\_sopra.fasta > Split_map/\$i.sam";
}
elsif ($mapper eq 'bowtie2'){
    $index_cmd = "$bowtie2-build OUT/contigs_sopra.fasta OUT/contigs_sopra.fasta";
    $map_cmd = "$bowtie2 $options{mapper_ops} -x OUT/contigs_sopra.fasta -f -U OUT/r.\$i\\\_sopra.fasta -S Split_map/\$i.sam";
}

# write and run bash script to run all the stages
open F, ">$run_script" or die $!;

print F "set -e
echo \"RUN s_prep_contigAseq_v1.4.6.pl\"
perl $sopra_bin/s_prep_contigAseq_v1.4.6.pl -contig $contigs -mate $split_reads_files -a OUT

echo \"RUN Bowtie build on the contigs\"
$index_cmd
";

unless ($options{noclean}) {
    print F "rm -rf Split_reads/\n";
}


print F "
echo \"RUN Map and run s_parse_sam_v1.4.6.pl  for each split file\"
mkdir Split_map
for i in {1..$last_file_index}
do
    $map_cmd
    perl $sopra_bin/s_parse_sam_v1.4.6.pl -sam Split_map/\$i.sam -a OUT -p \$i

done

echo \"RUN 03.s_read_parsed_sam_v1.4.6.pl\"
perl $sopra_bin/s_read_parsed_sam_v1.4.6.pl $sam_parsed_files -a OUT -pt $last_file_index

echo \"RUN s_scaf_v1.4.6.pl\"
perl $sopra_bin/s_scaf_v1.4.6.pl $options{s_scaf_ops} -o OUT/orientdistinfo_c5 -a OUT
";


if ($options{noclean}) {
    print F "ln -s OUT/scaffolds_*.fasta scaffolds.fa\n";
}
else {
    print F "mv OUT/scaffolds_*.fasta scaffolds.fa\n";
    print F "rm -fr Split_map/ OUT/\n";
}

close F;

exec "bash $run_script" or die $!;
