#!/usr/local/bin/perl
use strict;

my @c;
my @s;
my %bonds = (GU=>1,UG=>1,AU=>2,UA=>2,CG=>3,GC=>3);

sub max {
    my ($x,$y) = @_;
    if ($x>$y) { return $x; } else { return $y};
}

sub foldRna {
    my ($s) = @_;
    my $slen = length $s;
    @s = ('X', split(//, $s));
    
    for (my $len = 5; $len <= $slen; $len++) { 
	for (my $i=1; $i<=$slen-$len+1; $i++) {
	    my $j = $i+$len-1;
	    $c[$i][$j] = max($c[$i+1][$j],
			     $bonds{$s[$i].$s[$j]}+$c[$i+1][$j-1]);
	    
	    for (my $k=$i+1; $k<$j; $k++) {
		$c[$i][$j] = max($c[$i][$j], 
                                 $c[$i][$k]+$c[$k+1][$j]);
	    }
	}
    }
}

sub traceBack {
    my ($i,$j) = @_;
    my $cij = $c[$i][$j];

    return ("."  x ($j-$i+1)) if ($cij==0);
    return "." . traceBack($i+1,$j) 
        if ($cij == $c[$i+1][$j]);
    return "(" . traceBack($i+1,$j-1) . ")" 
	if ($cij == $bonds{$s[$i].$s[$j]}+ $c[$i+1][$j-1]);
    for (my $k = $i+1; $k < $j; $k++) {
	return traceBack($i,$k) . traceBack($k+1,$j)
	    if ($cij == ($c[$i][$k]+$c[$k+1][$j]));
    }
}

$|=1;
my $basestring = <STDIN>;
chomp($basestring);
foldRna($basestring);
print "$basestring\n", traceBack(1, length $basestring), "\n";
