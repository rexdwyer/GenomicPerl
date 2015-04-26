#!/usr/bin/perl -I . -I /home/dwyer/book/chap6

use strict;
use Prosite;
use PrositeSuffixTree;
use SeqReader;

my $seqFile = new SeqReader("turkey.fasta");
my $allSeqs;
while (my $protein = $seqFile->readSeq()) {
    $allSeqs .= "<$protein>"
}
$seqFile->close();

my $ptree = new PrositeSuffixTree($allSeqs."=");

my $prosite = new Prosite("prosite.dat");
while (my $motif = $prosite->readMotif()) {
    $ptree->listMotifs($motif);
}
$prosite->close();


