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
		my $diff = $x[$i][$k] - $x[$j][$k];
		$sum += $diff*$diff;
	    }
	    $val += $c[$i][$j] * sqrt($sum);
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
    $max = 1;

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
    for (my $i=0; $i<@c; $i++) { $x[$i][0]=-1; }
    $x[0][0]=1;
    return;

    for (my $i=0; $i<$#c; $i++) {
	$x[$i][$i] = 1;
	$x[$#c][$i] = - 1.0 / sqrt($#c);
    }
}


readC();
initX();
printX();



my $objective = objectiveValue();
my $delta = 0.5;

while ($delta > 1E-4) {
    print "************ Objective Value = $objective; Delta=$delta\n";
    my $change=1;
    while ($change--) {
	print "Change=$change\n";
	for (my $i = 1; $i<@c; $i++) {
	    my @Ci;
	    my $sumSq;
	    my $best = undef; 
	    my $bestdelta;
	    my $bestObjective = $objective;
	    my @oxi;
	    my $sumSq;
	    for (my $j=0; $j<$#c; $j++) { 
		$oxi[$j] = $x[$i][$j];
		$sumSq += $oxi[$j]*$oxi[$j];
	    }
	    
	    for (my $k=0; $k<$#c; $k++) {  ## move vector i in direction k.

		foreach my $sdelta ($delta,-$delta) {
		    
		    for (my $j=0; $j<$#c; $j++) { $x[$i][$j] = $oxi[$j]; }
		    $x[$i][$k] += $delta;
		    my $newSumSq 
			= $sumSq - $oxi[$k]*$oxi[$k] + $x[$i][$k]*$x[$i][$k];
		    my $len = sqrt($newSumSq);
		    next if $len==0;
		    for (my $j=0; $j<$#c; $j++) { $x[$i][$j] /= $len; }
		    my $newObjective = objectiveValue();
		    ($best,$bestdelta,$bestObjective,$change) = 
			($k,$sdelta,$newObjective,1)
			if ($newObjective > $bestObjective);
		}
	    }
	    if (defined($best)) {
		print "$i,$best,$bestdelta\n";
		for (my $j=0; $j<$#c; $j++) { $x[$i][$j] = $oxi[$j]; }
		$x[$i][$best] += $bestdelta;
		my $newSumSq = $sumSq - $oxi[$best]*$oxi[$best] 
		    + $x[$i][$best]*$x[$i][$best];
		my $len = sqrt($newSumSq);
		for (my $j=0; $j<$#c; $j++) { $x[$i][$j] /= $len; }
		$objective = objectiveValue();		
	    }
	    print "$i: Objective Value = $objective\n";
	    printX();
	}
	print "Objective Value = $objective; Delta=$delta\n";
	printX();
    }
    $delta /= 2.0;
}
    
    
