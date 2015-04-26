#!/usr/bin/perl -I /home/dwyer/book/perl
use strict;
use Util;

## An implementation of the Needleman-Wunsch algorithm
## for global alignment of DNA sequences.


my @M;    ## scoring matrix
my $g = -2;   ## gap penalty

sub p {  ## match / mismatch reward/penalty
    my ($aa1, $aa2) = @_;
    return ($aa1 eq $aa2)?1:-1;
}

sub similarity {
    my($s,$t) = @_;
    foreach my $i (0..length($s)) { 
	foreach my $j (0..length($t)) { 
	    $M[$i][$j] = -999999;
	}
    }
    $M[0][0]=0;
    foreach my $layer (0..(length($s)+length($t))) { 
	my $i = min($layer,length($s)); 
	my $j = $layer-$i;
	while ($j<=length($t) && $i>=0) {
	    $M[$i+1][$j] = max($M[$i+1][$j], $M[$i][$j]+$g);
	    $M[$i][$j+1] = max($M[$i][$j+1], $M[$i][$j]+$g);
	    my $p =  p(substr($s,$i,1),substr($t,$j,1));
	    $M[$i+1][$j+1] = max($M[$i+1][$j+1], $M[$i][$j]+$p);
	    $i--; $j++;
	}
    }
    return ( $M[length($s)][length($t)] );
}

sub getAlignment {
    my ($s,$t) = @_;
    my ($i,$j) = (length($s), length($t));
    return ( "-"x$j, $t) if ($i==0);
    return ( $s, "-"x$i) if ($j==0);
    my ($sLast,$tLast) = (substr($s,-1),substr($t,-1));
    
    if ($M[$i][$j] == $M[$i-1][$j-1] + p($sLast,$tLast)) { ## Case 1
        ## last letters are paired in the best alignment
	my ($sa, $ta) = getAlignment(substr($s,0,-1), substr($t,0,-1));
	return ($sa . $sLast , $ta . $tLast );
    } elsif ($M[$i][$j] == $M[$i-1][$j] + $g) { ## Case 2
        ## last letter of the first string is paired with a gap
	my ($sa, $ta) = getAlignment(substr($s,0,-1), $t);
	return ($sa . $sLast , $ta . "-");
    } else { ## Case 3: last letter of the 2nd string is paired with a gap
	my ($sa, $ta) = getAlignment($s, substr($t,0,-1));
	return ($sa . "-" , $ta . $tLast );
    }
}

##MAIN
{
    my $s = <STDIN>; chomp($s);
    my $t = <STDIN>; chomp($t);
    print "Similarity score: ", similarity($s,$t), "\n";
    print "Alignment: \n";
    foreach my $x (getAlignment($s,$t)) {
	print $x,"\n";
    }
}
