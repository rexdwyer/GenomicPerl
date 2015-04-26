#!/usr/bin/perl -I /home/dwyer/book/perl -I .

use strict;
use PAMBuilder;

my $A = collectCounts();
printMatrix($A,"Raw Counts");
my $A = addPseudocounts($A);
printMatrix($A,"Raw Counts plus Pseudocounts");
my $P = freqsFromCounts($A);
printMatrix($P,"Frequencies");
my $M = ratesFromFreqs($P);
my $PAMM = rescaleRates($M, 0.01);
for (my $i=1; $i<=256; $i+=$i) {
    print "Mutation rate of PAM$i is ", averageMutationRate($PAMM), "\n";
    printMatrix(lodsFromFreqs(freqsFromRates($PAMM)),"PAM$i (half-bits)");
    $PAMM = multiplyRates($PAMM,$PAMM);
}



