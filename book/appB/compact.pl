#!/usr/bin/perl -I /home/dwyer/book/perl

use strict;
use Util;

sub compactSimilarity {
    my ($s,$t) = @_;
    my ($sLen, $tLen) = (length($s),length($t));
    my @a;
    $#a = $tLen;  ## fast allocate.
    foreach my $j (0..$tLen) { $a[$j] = -2*$j; }
    foreach my $i (1..$sLen) {
	unshift @a, -2*$i;
	foreach my $j (1..$tLen) {
	    ## Loop Invariant:
	    ## a[k] = M[i-1][k-1]  for k>=j; a[k] = M[i][k] for k<j;
	    my $m =  (substr($s,$i-1,1) eq substr($t,$j-1,1)) ? +1 : -1;
	    $a[$j] = max($a[$j]+$m, $a[$j-1]-2, $a[$j+1]-2);
	}
	pop @a;
    }
    return(@a);
}

sub compactAlignment {
    my($s,$t) = @_;
    my ($sLen, $tLen) = (length($s),length($t));
    if ($sLen == 1) {
	if ($t =~ /^(.*)$s(.*)$/) {
	    my $ss = "-"x length($1) . $s . "-"x length($2);
	    my $score = 3 - 2*$tLen;
	    return "$ss\n$t\n$score";
	} else {
	    my $score = 1 - 2*$tLen;
	    return "-"x($tLen-1) . "$s\n$t\n$score";
	}
    } else {
	my $mid = int($sLen/2);
	my $sBeg = substr($s,0,$mid);
	my $sEnd = substr($s,$mid, $sLen-$mid);
        my @aBeg = compactSimilarity($sBeg, $t);
        my @aEnd = reverse(compactSimilarity((scalar reverse $sEnd), 
					  (scalar reverse $t)));
	my ($kMax, $kMaxVal) = (0, $aBeg[0]+$aEnd[0]);
        foreach my $k (1..$tLen) { 
	    ($kMax, $kMaxVal) = ($k, $aBeg[$k]+$aEnd[$k])
		if ($aBeg[$k]+$aEnd[$k]) > $kMaxVal;
	}
	my ($ssBeg,$ttBeg,$scoreBeg) = 
	    split("\n", compactAlignment($sBeg, substr($t,0,$kMax)));
	my ($ssEnd,$ttEnd,$scoreEnd) = 
	    split("\n", compactAlignment($sEnd, substr($t,$kMax,$tLen-$kMax)));
        my $score = $scoreBeg + $scoreEnd;
	return "$ssBeg$ssEnd\n$ttBeg$ttEnd\n$score";
    }
}

#MAIN
{
    while (my $s = <DATA>) {
	chomp($s);
	my $t = <DATA>; chomp($t);
	print "\n", compactAlignment($s,$t), "\n";
    }
}
