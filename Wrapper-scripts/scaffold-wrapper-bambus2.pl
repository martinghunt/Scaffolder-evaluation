#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Spec;
use File::Spec::Link;
use File::Which;

my %options;
my $options_ok = GetOptions(\%options,
    'noclean'
);

if (!($#ARGV == 4 or $#ARGV == 7) or !($options_ok)) {
    print STDERR "usage: $0 [options] <contigs.fa> <output_directory> <1.bam> <1.insert_size> <1.insert_stdev> [<2.bam> <2.insert_size> <2.insert_stdev>]

where extra in [square brackets] can be used with a second library.

Important Notes:

1. This script has been tested with AMOS version 3.1.0.
   The bin/goBambus.py script in the AMOS download has bugs and needs replacing.
   This wrapper script has been tested with this version of goBambus2.py:
   http://sourceforge.net/p/amos/code/ci/38c22fccba9eb66b1d40fa3aa62b267a5da7852b/tree/src/Pipeline/goBambus2.py
   You need to replace the original AMOS bin/goBambus2.py script with the newer one.

2. This script also assumes that abyss-samtoafg from ABySS is in your path.

Options:

-noclean
    Use this to not clean up files.
    Default is to delete most files when finished
";

    exit(1);
}

my $contigs = File::Spec->rel2abs($ARGV[0]);
my $outdir = $ARGV[1];
my @bam_files;
my @inserts;
my @stdevs;
my %contig_name_to_id;

for (my $i=2; $i<=$#ARGV; $i+=3) {
    push @bam_files, File::Spec->rel2abs($ARGV[$i]);
    push @inserts, $ARGV[$i+1];
    push @stdevs, $ARGV[$i+2];
}


my $bank_transact = which('bank-transact');
$bank_transact or die "Error finding AMOS bin";
$bank_transact = File::Spec::Link->resolve($bank_transact);
my (undef, $amos_bin) = fileparse($bank_transact);
my $go_bambus = File::Spec->catfile($amos_bin, 'goBambus2.py');
-e $go_bambus or die "Error finding goBambus2.py script";

mkdir $outdir or die "Error mkdir $outdir";
chdir $outdir or die "Error chdir $outdir";

system_call("fastaq_to_fasta -l 0 $contigs tmp.contigs.fa");


open F, 'tmp.contigs.fa' or die "Error opening tmp.contigs.fa";
while (<F>) {
    if (/^>/) {
        my ($name) = split;
        $name = substr($name, 1);
        if ($contig_name_to_id{$name}) {
            die "Got name $name twice. Cannot continue";
        }
        $contig_name_to_id{$name} = 1 + scalar (keys %contig_name_to_id);
    }
}
close F;


if ($#bam_files) {
    my ($fragments, $contig_names1) = split_afg("samtools view $bam_files[0] | abyss-samtoafg -e 1 -i 1 -m $inserts[0] -s $stdevs[0] tmp.contigs.fa - |", 'tmp.1.afg.split', 0, \%contig_name_to_id);
    my (undef, $contig_names2) = split_afg("samtools view $bam_files[1] | abyss-samtoafg -e 2 -i 2 -m $inserts[1] -s $stdevs[1] tmp.contigs.fa - |", 'tmp.2.afg.split', 2 * $fragments, \%contig_name_to_id);
    system_call("cat tmp.1.afg.split.LIB tmp.2.afg.split.LIB > tmp.afg");
    system_call("cat tmp.1.afg.split.RED.FRG tmp.2.afg.split.RED.FRG >> tmp.afg");
    my %contig_names;
    foreach (keys %{$contig_names1}) {
        $contig_names{$_} = 1;
    }
    foreach (keys %{$contig_names2}) {
        $contig_names{$_} = 1;
    }

    for my $id (keys %contig_names) {
        if ($contig_names1->{$id}) {
            system_call("cat tmp.1.afg.split.CTG.$id tmp.1.afg.split.TLE.$id >> tmp.afg");
            if ($contig_names2->{$id}) {
                system_call("cat tmp.2.afg.split.TLE.$id >> tmp.afg");
            }
        }
        else {
            system_call("cat tmp.2.afg.split.CTG.$id tmp.2.afg.split.TLE.$id >> tmp.afg");
        }
        system_call("echo } >> tmp.afg");
    }

}
else {
    system_call("samtools view $bam_files[0] | abyss-samtoafg -e 1 -i 1 -m $inserts[0] -s $stdevs[0] tmp.contigs.fa - > tmp.afg");
}


system_call("bank-transact -cb tmp.bnk -m tmp.afg");
system_call("python $go_bambus tmp.bnk scaff");
system_call(q~awk '!/^#/ && $6 != "fragment" {print $6}' tmp.scaff.linear.agp > tmp.used_contigs.ids~);
system_call("fastaq_filter -v --ids_file tmp.used_contigs.ids tmp.contigs.fa tmp.unused_contigs.fa");
system_call("cat scaff.scaffold.linear.fasta tmp.unused_contigs.fa > scaffolds.fa");

unless ($options{noclean}) {
    system_call("rm -fr tmp.* myreps scaff.*\n");
}


sub split_afg {
    my $infile = shift;
    my $outprefix = shift;
    my $offset = shift;
    my $contig_name_to_id = shift;
    my $fragments = 0;
    my %contig_names;

    open FIN, $infile or die $!;

    my $change_ids = 0;
    my %contig_TLE;

    my $line = <FIN>;

    if ($line eq "{LIB\n") {
        open FOUT, ">$outprefix.LIB" or die "Error opening $outprefix.LIB";
        while ($line ne "{RED\n") {
            print FOUT $line;
            $line = <FIN>;
        }
        close FOUT;
    }
    else {
        die "Error getting library info from second library";
    }

    if ($line eq "{RED\n") {
        open FOUT, ">$outprefix.RED.FRG" or die "Error opening $outprefix.RED.FRG";
        while ($line ne "{CTG\n") {
            if ($line =~ /^(iid:|frg:)(\d+)\n$/) {
                $line = $1 . ($2 + $offset) . "\n";
                if ($1 eq "frg:") {
                    $fragments++;
                }
            }
            print FOUT $line;
            $line = <FIN>;
        }
        close FOUT;
    }
    else {
        die "Error getting sequence and fragment info from second library";
    }

    while ($line eq "{CTG\n") {
        $line = <FIN>;
        if ($line =~ /^iid:\d+\n$/) {
            my $iid_line = $line;
            $line = <FIN>;
            my $contig_name;
            if ($line =~ /^eid:(\S+)\n$/) {
                $contig_name = $1;
            }
            else {
                die "Error getting eid from contig";
            }

            $contig_name_to_id->{$contig_name} or die "Got contig '$contig_name' in BAM but not in contigs file";
            $contig_names{$contig_name} = 1;
            open FOUT, ">$outprefix.CTG.$contig_name" or die "Error opening $outprefix.CTG.$contig_name";
            print FOUT "{CTG\n";
            print FOUT "iid:" . $contig_name_to_id->{$contig_name} . "\n";
            print FOUT $line;
            my @lines;

            while($line ne "{TLE\n") {
                $line = <FIN>;
                push @lines, $line;
            }
            pop @lines;
            foreach (@lines) {
                print FOUT $_;
            }

            open FOUT, ">$outprefix.TLE.$contig_name" or die "Error opening $outprefix.TLE.$contig_name";
            while ($line eq "{TLE\n") {
                for (0..4) {
                    print FOUT $line;
                    $line = <FIN>;
                }
            }
            close FOUT;
            $line eq "}\n" or die $!;
        }
        else{
            die "Error getting contig iid around this line: $line";
        }

        if (eof(FIN)) {
            last;
        }
        $line = <FIN>;
    }

    close FIN or die $!;
    return $fragments, \%contig_names;
}


sub system_call {
    my $cmd  = shift;
    if (system($cmd)) {
        print STDERR "Error in system call:\n$cmd\n";
        exit(1);
    }
}
