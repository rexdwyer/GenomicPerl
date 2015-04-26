#!/usr/bin/perl -I /home/dwyer/book/perl

use strict;
use SeqReader;
use Util;
# An implementation of the Needleman-Wunsch algorithm
# for global alignment of DNA sequences.


my @M;    # scoring matrix
my $g = -2;   # gap penalty

sub p {
    my ($sAa,$tAa) = @_;
    return (($sAa eq $tAa)?1:-1);
}

sub similarity {
    my($s,$t) = @_;
    undef @M;
    my ($i,$j);
    my ($sLen,$tLen) = (length $s, length $t);

    for (my $i=$sLen; $i>=0; $i--) { $M[$i][length $t] = undef; }
    my @js = (1..$tLen);
    foreach $i (0..$sLen) { $M[$i][0] = $g * $i; }
    foreach $j (@js) { $M[0][$j] = $g * $j; }
    my @ss = ("_", (split //, $s));
    my @ts = ("_", (split //, $t));
    foreach $i (1..$sLen) {
	my $sI = $ss[$i];
	my $i1=$i-1; 
	foreach $j (@js) {
	    my $j1=$j-1;
	    if ($sI eq $ts[$j]) {
		$M[$i][$j] = $M[$i1][$j1] + 1;
	    } else {
		my $best = $M[$i1][$j1] - 1;
		my $try = $M[$i1][$j] + $g;
		$best = $try if $try > $best;
		$try = $M[$i][$j1] + $g;
		$best = $try if $try > $best;
		$M[$i][$j] = $best;
	    }
	}
    }
    return ( $M[$sLen][$tLen] );
}

sub getAlignment {
    my($s,$t) = @_;
    my $i = length($s);
    my $j = length($t);
    return ( "-"x$j, $t) if ($i==0);
    return ( $s, "-"x$i) if ($j==0);
    
    my $sLast = substr($s,-1);
    my $tLast = substr($t,-1);
    if ($M[$i][$j] == $M[$i-1][$j-1] + p($sLast,$tLast)) {
	my ($sa, $ta) = getAlignment(substr($s,0,-1), substr($t,0,-1));
	return ($sa . $sLast , $ta . $tLast);
    } elsif ($M[$i][$j] == $M[$i-1][$j] + $g) {
	my ($sa, $ta) = getAlignment(substr($s,0,-1), $t);
	return ($sa . $sLast , $ta . "-");
    } else {
	my ($sa, $ta) = getAlignment($s, substr($t,0,-1));
	return ($sa . "-" , $ta . $tLast);
    }
}

#MAIN
$|=1;
my $turkeyFile = new SeqReader "./turkey.fasta";
while (my @turkeySeq = $turkeyFile->readSeq()) {
    my $turkeyProtein = @turkeySeq[0];
    my $tLen = length $turkeyProtein;
    my $catFile = new SeqReader "./cat.fasta";
    my @best = (-100000,undef,undef);
    while (my @catSeq = $catFile->readSeq()) {
	my $catProtein = @catSeq[0];
	my $cLen = length $catProtein;
	my $score = "*";
	if ($best[0] < ($g*abs($cLen-$tLen) + min($cLen,$tLen))) {
	    $score = similarity($turkeyProtein,$catProtein);
	    @best = ($score, $turkeySeq[1], $catSeq[1], 
		     getAlignment($turkeyProtein,$catProtein))
		if $score >$best[0];
	}
	print STDERR "$score ";
    }
    $catFile->close();
    warn "\n$best[1]\n    $best[2]\n";
    map((print "$_\n"), @best);
}
$turkeyFile->close();
