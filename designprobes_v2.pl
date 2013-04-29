#!/usr/bin/perl -w
use strict;

# Tool for designing probes based on a sliding window
# Cuts sequence every $probesize base pairs
# Input: sequence file, probe size, tiling density, output file
# Katrina Galkina

# Usage
# designprobes_v2.pl [filename] [probesize] [tilingdensity] [outfile]

my $numArgs = $#ARGV + 1;
my $filename_in = $ARGV[0];
my $probesize = $ARGV[1];
my $tilingdensity = $ARGV[2];
my $filename_out = $ARGV[3];

open (SEQ, $filename_in) || die "cannot open $filename_in.\n";
open (OUT, "> ".$filename_out) || die "cannot open $filename_out.\n";

my @bases;

# Input file should not contain spaces or carriage returns or a header (resembling FASTA format, but without header)
while (<SEQ>) {

    chomp;

    my $len = length($_);

#Probe generation will start with the first base pair of your input sequence 
#To change this, set $i to the  starting base pair of your choice (instead of the default 1)
#Set $tilingdensity to 1 to move one base pair at a time. 
#To move more than 1 base pair along the sequence, set $tilingdensity to that number

    for( my $i = 1; $i <= $len-($probesize-1); $i+=$tilingdensity) {
        print OUT substr($_, $i, $probesize)."\n";  #create a substring of sequence starting with i (default is 1) ending in probesize
    }

}

close (SEQ);
close (OUT);

