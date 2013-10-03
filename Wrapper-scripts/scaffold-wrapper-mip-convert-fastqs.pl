#!/usr/bin/env perl

use strict;
use warnings;

if ($#ARGV != 2) {
    die qq/usage:$0 <in_1.fastq> <in_2.fastq> <out_prefix>

Converts forward and reverse reads fastq files into format for
the MIP scaffolder. Assumes reads are innies.
/;
}

my $in_1 = $ARGV[0];
my $in_2 = $ARGV[1];
my $out_1 = $ARGV[2] . "_R3.fastq";
my $out_2 = $ARGV[2] . "_F3.fastq";

# Here's what MIP needs:
# The fwd reads must become the reverse reads, with /1 at
# the end of each name becoming _R3
#
# The rev reads become the fwd reads, they must be
# reverse complemented and /2 in each name becomes _F3

open IN, $in_1 or die $!;
open OUT, ">$out_1" or die $!;

while (my $line = <IN>) {
    $line =~ /^@/ or die "Expected line starting with '\@'. Got this:\n$line";
    chomp $line;

    print OUT substr($line, 0, -2) . "_R3\n";

    for my $i (0..2) {
         $line = <IN>;
         print OUT $line;
    }
}

close IN;
close OUT;

open IN, $in_2 or die $!;
open OUT, ">$out_2" or die $!;

while (my $line = <IN>) {
    $line =~ /^@/ or die "Expected line starting with '\@'. Got this:\n$_";
    chomp $line;
    print OUT substr($line, 0, -2) . "_F3\n";
    $line = <IN>;
    chomp $line;
    $line =~ tr/ACGTacgt/TGCAtgca/;
    $line = reverse($line);
    print OUT "$line\n";
    $line = <IN>;
    print OUT "+\n";
    $line = <IN>;
    chomp $line;
    print OUT reverse($line) . "\n";
}

close IN;
close OUT;
