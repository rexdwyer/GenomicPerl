#!/usr/bin/perl
use strict; 

my @c;
my @x;

sub objectiveValue {
    my $val;
    for (my $i=0; $i<@c; $i++) {
	for (my $j=0; $j<$i; $j++) {
	    next if $c[$i][$j] == 0;
	    my $sum;
	    for (my $k=0; $k<$#c; $k++) {
		$sum += $x[$i][$k] * $x[$j][$k];
	    }
	    $val += $c[$i][$j] * $sum;
	}
    }
    return $val;
}

sub printX {
    print "\n";
    for (my $i=0; $i<@c; $i++) {
	for (my $j=0; $j<$#c; $j++) {
	    printf "%9.6f ", $x[$i][$j];
	}
	print "\n";
    }
    print "\n";
}
	    
sub readC {
    while (0 && <STDIN>) {
#	my ($l1,$l2,$r1,$r2,$confid) = m/(\d*),(\d*);(\d*),(\d*):(\d*)/;
	my ($l1,$l2,$r1,$r2) = m/(\d*),(\d*);(\d*),(\d*)/;
	my $confid = 1;
	$c[$l1][$l2] += 2*$confid;
	$c[$r1][$r2] += 2*$confid;
	$c[$l1][$r1] -= $confid;
	$c[$l1][$r2] -= $confid;
	$c[$l2][$r1] -= $confid;
	$c[$l2][$r2] -= $confid;
    }

    my $max;
    for (my $i=0; $i<@c; $i++) {
	for (my $j=0; $j<$i; $j++) {
	    $c[$i][$j] += $c[$j][$i];
	    $c[$j][$i] = $c[$i][$j];
	    $max = $c[$i][$j] if $c[$i][$j]>$max;
	    printf "%5d", 	$c[$i][$j];
	}
	print "\n";
    }
    print "max=$max\n";
    $c[4][4]=1;

    $max *=100;
    $max = -1;

    for (my $i=0; $i<@c; $i++) {
	for (my $j=0; $j<$i; $j++) {
	    $c[$i][$j] = $max;
	    $c[$j][$i] = $max;
	    printf "%5d", $c[$i][$j];
	}
	print "\n";
    }
}

sub initX {
    for (my $i=0; $i<$#c; $i++) {
	$x[$i][$i] = 1;
    }
    $x[$#c][$#c-1] = 1;
}


readC();
initX();
printX();

my $objectiveDelta = 1;
my $objective = objectiveValue();

while ($objectiveDelta > 1E-4) {
    for (my $i = 0; $i<@c; $i++) {
	my @Ci;
	my $sumSq;
	for (my $k=0; $k<$#c; $k++) {
	    my $Cik = 0;
	    for (my $j=0; $j<@c; $j++) {
		$Cik += $c[$i][$j] * $x[$j][$k];
	    }
	    $Ci[$k] = $Cik;
	    $sumSq += $Cik*$Cik;
	}
	my $normalization = sqrt($sumSq);
	for (my $k=0; $k<$#c; $k++) {
	    $x[$i][$k] = $Ci[$k] / $normalization;
	}
	print "Objective Value = ", objectiveValue(), "\n";
    }
    (my $oldObjective,$objective) = ($objective,objectiveValue());
    $objectiveDelta = $objective - $oldObjective;
    printX();
    print "Objective Value = $objective\n";

    for (my $i=0; $i<@c; $i++) {
	for (my $j=0; $j<@c; $j++) {
	    my $sum;
	    for (my $k=0; $k<@c; $k++) {
		$sum += $x[$i][$k] * $x[$j][$k];
	    }
	    printf "%8.6f ", $sum;
	}
	print "\n";
    }
}

	
