#!/usr/bin/perl -I /home/dwyer/book/perl

## An implementation of the recursive, space-saving algorithm
## for global alignment, S&M p. 59 & 16, simplified.

use strict;
use Util;

my ($g,$h) = (-2,-3);

sub maxIx {
    my($m,@l) = @_;
    my $mi=0; 
    my $i=0;
    foreach my $x (@l) { 
	$i++;
	if ($x > $m) {
	    $m=$x;
	    $mi=$i;
	}
    }
    return $mi;
}

sub compactScore {
    my ($s,$t) = @_;
    my ($i,$j,@a,@b,@c,$olda, $oldb, $oldc, $tempa, $tempb, $tempc);
    my $tlen = length($t);
    
    ##initialize row  [0,*]
    undef (@a,@b,@c); 
    ($#a,$#b,$#c) = ((0 x $tlen) x 3);
    foreach $j (0..$tlen) { 
	($a[$j],$b[$j],$c[$j]) = (-(1<<30), $h+$g*$j, -(1<<30));
    }
    $a[0] = 0;
    
    
    ##fill in row [$i,*]
    foreach $i (1..length($s)) {
        ($olda, $oldb, $oldc) = ($a[0], $b[0], $c[0]);
	($a[0], $b[0], $c[0]) = (-(1<<30), -(1<<30), $h + $g*$i);
	foreach $j (1..$tlen) {
            ($tempa,$tempb,$tempc)  = ($a[$j], $b[$j], $c[$j]);
            my $m = (substr($s,$i-1,1) eq substr($t,$j-1,1))?1:-1;
	    
	    $a[$j] = $m + &max($olda,$oldb,$oldc);
	    $b[$j] = $g + &max($a[$j-1]+$g, $b[$j-1], $c[$j-1]+$g);
	    $c[$j] = $g + &max($tempa+$h, $tempb+$h, $tempc);
	    ($olda,$oldb,$oldc) = ($tempa,$tempb,$tempc);
	}
    }
    return(\@a,\@b,\@c);
}

sub revString() {
    return join('',reverse(split(//,$_[0])));
}

sub recursiveAlign {
    my($s,$t) = @_;
#	print ">recursiveAlign($s,$t)\n";
    my $slen = length($s);
    my $tlen = length($t);
    if ($slen == 0) {
	return( ("-" x $tlen) . "\n$t\n" . ($g*$tlen) );
    } elsif ($tlen == 0) {
	return( "$s\n" . ("-" x $slen) . "\n" . ($g*$slen) );
    } elsif ($slen == 1) {
	my $k = index $t, $s;
	my $ss = ("-" x  &max($k,0)) . $s;
	$ss .=  ("-" x ($tlen - length($ss)));
	my $score = (($tlen-1) * $g) + (($k>=0) ? 1 : -1);
	my $result = "$ss\n$t\n$score";
#		print "<recursiveAlign($s,$t) = \n$result\n";
	return $result;
    } else {
	my $mid = int($slen/2);
	my $sBeg = substr($s,0,$mid);
	my $sEnd = substr($s,$mid);
	my @refs = &compactScore($sBeg, $t);
        my @aBeg = @{$refs[0]};
        my @bBeg = @{$refs[1]};
        my @cBeg = @{$refs[2]};
	my @refs = &compactScore(&revString($sEnd),&revString($t));
        my @aEnd = reverse(@{$refs[0]});
        my @bEnd = reverse(@{$refs[1]});
        my @cEnd = reverse(@{$refs[2]});
#	print join(',',@aBeg), "A\n";
#	print join(',',@bBeg), "B\n";
#	print join(',',@cBeg), "C\n";
#	print join(',',@aEnd), "A\n";
#	print join(',',@bEnd), "B\n";
#	print join(',',@cEnd), "C\n";
	
	my ($k,@a,@ix);
        foreach $k (0..$tlen) { 
	    my @options = 
		($aBeg[$k]+$aEnd[$k],
		 $aBeg[$k]+$bEnd[$k],
		 $aBeg[$k]+$cEnd[$k],
		 $bBeg[$k]+$aEnd[$k],
		 $bBeg[$k]+$bEnd[$k]-$h,
		 $bBeg[$k]+$cEnd[$k],
		 $cBeg[$k]+$aEnd[$k],
		 $cBeg[$k]+$bEnd[$k],
		 $cBeg[$k]+$cEnd[$k]-$h);
	    $ix[$k] = &maxIx(@options);
	    $a[$k] = $options[$ix[$k]];
	}
#	print join(',',@a), "AA\n";
#	print join(',',@ix), "IX\n";
	my $k = &maxIx(@a);
	my $ix = $ix[$k];
	my $score = $a[$k];
	my $begIx = int($ix / 3);
	my $endIx = int($ix % 3);
	my ($ssBeg,$ttBeg,$scoreBeg);
	my ($ssEnd,$ttEnd,$scoreEnd);
#	print "$score,$ix,$begIx,$endIx\n";
	
	if ($begIx == 0) {
	    ($ssBeg,$ttBeg,$scoreBeg) = 
		split("\n", &recursiveAlign($sBeg, substr($t,0,$k)));
	} elsif ($begIx == 1) {
	    ($ssBeg,$ttBeg,$scoreBeg) = 
		split("\n", &recursiveAlign($sBeg, substr($t,0,$k-1)));
	    $ssBeg .= "-";
	    $ttBeg .= substr($t,$k-1,1);
	} elsif ($begIx == 2) {
	    ($ssBeg,$ttBeg,$scoreBeg) = 
		split("\n",&recursiveAlign(substr($sBeg,0,-1),
					    substr($t,0,$k)));
	    $ssBeg .= substr($sBeg,-1);  ## last letter
	    $ttBeg .= "-"
	    } else { print "CASE ERROR!\n"; }
	
	if ($endIx == 0) {
	    ($ssEnd,$ttEnd,$scoreEnd) = 
		split("\n", &recursiveAlign($sEnd, substr($t,$k)));
	} elsif ($endIx == 1) {
	    ($ssEnd,$ttEnd,$scoreEnd) = 
		split("\n", &recursiveAlign($sEnd, substr($t,$k+1)));
	    $ssEnd = "-" . $ssEnd;
	    $ttEnd = substr($t,$k,1) . $ttEnd;
	} elsif ($endIx == 2) {
	    ($ssEnd,$ttEnd,$scoreEnd) = 
		split("\n",&recursiveAlign(substr($sEnd,1),substr($t,$k)));
	    $ssEnd = substr($sEnd,0,1) . $ssEnd;
	    $ttEnd = "-" . $ttEnd;
	} else { print "CASE ERROR!\n"; }
	
	my $result = "$ssBeg$ssEnd\n$ttBeg$ttEnd\n$score";
#	print "<recursiveAlign($s,$t) = \n$result\n";
	return $result;
    }
}

##MAIN
{
    while (my $s = <DATA>) {
	chomp($s);
	my $t = <DATA>; chomp($t);
	print "\n", &recursiveAlign($s,$t), "\n";
    }
}

__END__
AB
ABC
ABC
AB
AC
ABC
ABC
AC
RECONSTRUCTION
UNCONDITIONALLY
THANKSGIVING
CHANGELING
ACURSEDFIENDBRINGSGRIEFANDPAIN
ABLESSEDFRIENDBRINGSRELIEFAGAIN
