#!/usr/bin/perl
use strict; 

my $node = -2;
my %parent;
$parent{-2}=-1;
    my %path;

my $size = 10;
generateRandomTree(0,$size-1,-1);
for (my $i=0; $i<$size; $i++) {
    undef %path;
    my $j=$i;
    my $cnt=1;
    while ($j != -1) {
	$path{$j} = $cnt++;
	$j = $parent{$j};
    }
    for (my $j=$i+1; $j<$size; $j++) {
	my $jcnt = find($j);

	for (my $k=$j+1; $k<$size; $k++) {

	    my $kcnt = find($k);

	    for (my $m=$k+1; $m<$size; $m++) {
		my $mcnt = find($m);

		my ($mincnt, $mininx) = ($jcnt,$j);
		($mincnt, $mininx) = ($kcnt,$k) if $kcnt<$jcnt;
		($mincnt, $mininx) = ($mcnt,$m) if $mcnt<$kcnt;
		
		print "$i,$j;$k,$m\n" if ($mininx == $j);
		print "$i,$k;$j,$m\n" if ($mininx == $k);
		print "$i,$m;$j,$k\n" if ($mininx == $m);
	    }
	}		
    }
}

sub find {
    my ($n) = @_;
    while (!defined($path{$n})) { $n = $parent{$n}; };
    return $path{$n};
}

sub generateRandomTree {
    my ($lo,$hi,$par) = @_;
    if ($lo==$hi) {
	print "parent{$lo} = $par\n";
	return( $parent{$lo} = $par );
    }
    my $me = $node--;
    $parent{$me} = $par;
    print "parent{$me} = $par\n";
    my $split = $lo + int(rand($hi-$lo));
    generateRandomTree($lo,$split,$me);
    generateRandomTree($split+1,$hi,$me);
}
