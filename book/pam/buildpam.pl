#!/usr/bin/perl -I /home/dwyer/book/perl -I .

use strict;
use Util;
use PAMBuilder;

my $A = collectCounts();
printMatrix($A,"Raw Frequencies");
my $A = addPseudocounts($A,1);
printMatrix($A,"Raw Frequencies plus Pseudocounts");
my $P = countsToProbs($A);
printMatrix($P,"Probabilities");

my $L = probsToLods($P);
printMatrix($L,"LOD Scores (half-bits)");

my $M = probsToRates($P);
my $PAMM1 = rescaleRates($M, 0.01);
print "Mutation rate of PAM1 is ", averageMutationRate($PAMM1), "\n";
printMatrix(probsToLods(ratesToProbs($PAMM1)),"PAM1 (half-bits)");
my $PAMM2 = multiplyRates($PAMM1,$PAMM1);
my $PAMM4 = multiplyRates($PAMM2,$PAMM2);
my $PAMM8 = multiplyRates($PAMM4,$PAMM4);
my $PAMM16 = multiplyRates($PAMM8,$PAMM8);
my $PAMM32 = multiplyRates($PAMM16,$PAMM16);
my $PAMM64 = multiplyRates($PAMM32,$PAMM32);
my $PAMM128= multiplyRates($PAMM64,$PAMM64);
my $PAMM256 = multiplyRates($PAMM128,$PAMM128);
printMatrix(probsToLods(ratesToProbs($PAMM32)),"PAM32 (half-bits)");
print "Mutation rate of PAM32 is ", averageMutationRate($PAMM32), "\n";
printMatrix($PAMM64,"PAM64 mutation rates");
printMatrix(ratesToProbs($PAMM64),"PAM64 probabilities");
printMatrix(probsToLods(ratesToProbs($PAMM64)),"PAM64 (half-bits)");
printMatrix(probsToLods(ratesToProbs($PAMM128)),"PAM128 (half-bits)");
printMatrix(probsToLods(ratesToProbs($PAMM256)),"PAM256 (half-bits)");
print "Mutation rate of PAM256 is ", averageMutationRate($PAMM256), "\n";
