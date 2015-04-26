#!/usr/bin/perl -I ../perl
use strict; 
use Util;

my @c;
my @x;
my $n;   ## number of species
my $d;   ## number of dimensions = n-1
my @dist;

sub innerProd {  ## inner product of two vectors.
    my ($u,$v) = @_;
    my $limit = min($#$u,$#$v);
    my $sum;
    for (my $j=0; $j<=$limit; $j++) { 
	$sum += $$u[$j]*$$v[$j];
    }
    return $sum;
}

sub distance {  ## Euclidean distance between two vectors.
    my ($u,$v) = @_;
    my $limit = max($#$u,$#$v);
    my $sumSq;
    for (my $j=0; $j<=$limit; $j++) { 
	my $diff = $$u[$j]-$$v[$j];
	$sumSq += $diff*$diff;
    }
    return sqrt($sumSq);
}

sub distanceSq {  ## Euclidean distance between two vectors.
    my ($u,$v) = @_;
    my $limit = max($#$u,$#$v);
    my $sumSq;
    for (my $j=0; $j<=$limit; $j++) { 
	my $diff = $$u[$j]-$$v[$j];
	$sumSq += $diff*$diff;
    }
    return $sumSq;
}

sub reducedRowEchelon {
    my ($M) = @_;
    my $rows = @$M;
    my $cols = 1+max(map($#{$_||[]}, @$M));
    my $maxI = min($rows,$cols);
    for (my $i=0; $i<$maxI; $i++) {  ## for each species !=i
	my $Mii = $$M[$i][$i];
	(warn "reducedRowEchelon: M[$i][$i] is 0\n", next) unless $Mii;
	for (my $j=0; $j<$rows; $j++) {  ## for each species !=i
	    next if $i==$j;
	    my $factor = $$M[$j][$i] / $Mii;
	    $$M[$j][$i] = 0;
	    for (my $k=$i+1; $k<$cols; $k++) {  ## for each dimension
		$$M[$j][$k] -= $factor * $$M[$i][$k];
	    }
	}
    }
}

sub objectiveValue {
    my ($c,$x) = @_;
    my $val;
    for (my $i=0; $i<$n; $i++) {
	for (my $j=0; $j<$i; $j++) {
	    next if $$c[$i][$j] == 0;
	    $val += $$c[$i][$j] * distanceSq($$x[$i],$$x[$j]);
	}
    }
    return $val;
}

sub printMatrix {
    my ($M, $format) = @_;
    my $rows = @$M;
    my $cols = 1+max(map($#{$_||[]}, @$M));
    $format = $format || "%12.3f ";
    print "\n";
    for (my $i=0; $i<$rows; $i++) {
	for (my $j=0; $j<$cols; $j++) {
	    printf $format, $$M[$i][$j];
	}
	print "\n";
    }
    print "\n";
}
	    
sub printLowerTriangularMatrix {
    my ($M, $format) = @_;
    my $rows = @$M;
    $format = $format || "%12.3f ";
    print "\n";
    for (my $i=0; $i<$rows; $i++) {
	for (my $j=0; $j<=$i; $j++) {
	    printf $format, $$M[$i][$j];
	}
	print "\n";
    }
    print "\n";
}
	    
sub readC {
    $d=0;
    while (<STDIN>) {
	my ($l1,$l2,$r1,$r2,$confid) = m/(\d*),(\d*);(\d*),(\d*):(\d*)/;
#	my ($l1,$l2,$r1,$r2) = m/(\d*),(\d*);(\d*),(\d*)/;
#	my $confid = 1;
	# Next two pairs should be close; penalize (subtract from objective) if far apart.
	$c[$l1][$l2] -= 2*$confid;
	$c[$r1][$r2] -= 2*$confid;
	# Next four pairs should be far; reward (add to objective) if far apart.
	$c[$l1][$r1] += $confid;
	$c[$l1][$r2] += $confid;
	$c[$l2][$r1] += $confid;
	$c[$l2][$r2] += $confid;
	$d = max($d,$l1,$l2,$r1,$r2);

    }
    $n = $d+1;

    for (my $i=0; $i<$n; $i++) {
	for (my $j=0; $j<=$i; $j++) {
	    $c[$i][$j] += $c[$j][$i];
	    $c[$j][$i] = $c[$i][$j];
	    printf "%5d", ($c[$i][$j] || 0.0);
	}
	print "\n";
    }
}

sub translate {
    my ($M,$v) = @_;  ## index of vector to be translated to origin.
    my $rows = @$M;
    my $cols = 1+max(map($#{$_||[]}, @$M));
    for (my $k=0; $k<$cols; $k++) {  ## for each dimension
	my $vk = $$v[$k];
	for (my $j=0; $j<$rows; $j++) {  ## for each species
	    $$M[$j][$k] -= $vk;
	}
    }
}

sub findHyperplaneEqn {
    my ($x, $i) = @_;
    my @M;
    my $jj=0;
    for (my $j=0; $j<$n; $j++) {  ## for each species
	next if $j==$i;           ## except the one being optimized
	for (my $k=0; $k<$d; $k++) {  ## for each dimension
	    $M[$jj][$k] = $$x[$j][$k];
	}
	$M[$jj][$d] = 1;
	$jj++;
    }
    reducedRowEchelon(\@M);

    my $prod=1;
    for (my $k=0; $k<$n-1; $k++) { $prod *= $M[$k][$k]; }
    my @a;
    for (my $j=0; $j<$d; $j++) { 
	$a[$j] = $prod * $M[$j][$d] / $M[$j][$j];
    }
    my $norm = distance(\@a,[]);
    for (my $j=0; $j<$d; $j++) { $a[$j] /= $norm; }
    return @a;
}

sub checkDistanceInvariant {
    my ($x) = @_;
    my $inv;
    for (my $i=0; $i<$n; $i++) {  ## for each species pair
	for (my $j=0; $j<$i; $j++) {
	    $inv += distanceSq($$x[$i], $$x[$j]);
	}
    }
    return $inv;
}

sub findMinPositive {
    my ($M, $mask) = @_;
    my $rows = @$M;
    my $cols = 1+max(map($#{$_||[]}, @$M));
    my $minEntry = 1E99;
    my @minPair;
    for (my $i=0; $i<$rows; $i++) {
	next if ($$mask[$i] > -1);
	for (my $j=0; $j<$cols; $j++) {
	    next if ($$mask[$j] > -1);
	    ($minEntry,@minPair) = ($$M[$i][$j],$i,$j)
		if ($$M[$i][$j] > 0 && $$M[$i][$j] < $minEntry);
	}
    }
    return @minPair;
}
    
sub constructMidpointTree {
    print "\nConstructing Tree by Midpoint Method\n";
    my @parent = ((-1) x $n);
    my @descendents = map([$_], (0..$n-1));
    for (my $avail = $n; $avail>1; $avail--) {
	my ($i,$j) = findMinPositive(\@dist, \@parent);
	push @parent, -1;
	my $nu = $#parent;
	$parent[$i] = $parent[$j] = $nu;
	$descendents[$nu] = [(@{$descendents[$i]}, @{$descendents[$j]})];
	print "descendents of $nu: {", join(',', @{$descendents[$nu]}), "}\n";
	for (my $k=0; $k<$d; $k++) {
	    $x[$nu][$k] = 0.5*($x[$i][$k]+$x[$j][$k]);
	}
	for (my $j=0; $j<$nu; $j++) {
	    $dist[$nu][$j] = distance($x[$nu],$x[$j]);
	}
    }
#   printLowerTriangularMatrix(\@dist, "%4d");
#   printMatrix(\@x);
}


sub centroid {
    my @cen;
    for (my $k=0; $k<$d; $k++) {
	my $coord;
	foreach (@_) { $coord += $$_[$k]; }
	$cen[$k] = $coord / @_;
    }
    return \@cen;
}

sub constructCentroidTree {
    print "\nConstructing Tree by Centroid Method\n";
    my @parent = ((-1) x $n);
    my @descendents = map([$_], (0..$n-1));
    for (my $avail = $n; $avail>1; $avail--) {
	my ($i,$j) = findMinPositive(\@dist, \@parent);
	push @parent, -1;
	my $nu = $#parent;
	$parent[$i] = $parent[$j] = $nu;
	$descendents[$nu] = [(@{$descendents[$i]}, @{$descendents[$j]})];
	print "descendents of $nu: {", join(',', @{$descendents[$nu]}), "}\n";
	my $weight = 1/@{$descendents[$nu]};
	$x[$nu] = centroid(map( $x[$_], @{$descendents[$nu]}) );
	for (my $j=0; $j<$nu; $j++) {
	    $dist[$nu][$j] = distance($x[$nu],$x[$j]);
	}
    }
#   printLowerTriangularMatrix(\@dist, "%4d");
#   printMatrix(\@x);
}

sub constructCircumcenterTree {
    print "\nConstructing Tree by Circumcenter Method\n";
    my @parent = ((-1) x $n);
    my @descendents = map([$_], (0..$n-1));
    for (my $avail = $n; $avail>1; $avail--) {
	my ($i,$j) = findMinPositive(\@dist, \@parent);
	push @parent, -1;
	my $nu = $#parent;
	$parent[$i] = $parent[$j] = $nu;
	$descendents[$nu] = [(@{$descendents[$i]}, @{$descendents[$j]})];
	print "descendents of $nu: {", join(',', @{$descendents[$nu]}), "}\n";
	my $weight = 1/@{$descendents[$nu]};
	$x[$nu] = circumcenter(map( $x[$_], @{$descendents[$nu]}) );
	for (my $j=0; $j<$nu; $j++) {
	    $dist[$nu][$j] = distance($x[$nu],$x[$j]);
	}
    }
#   printLowerTriangularMatrix(\@dist, "%9d");
#   printMatrix(\@x);
}

sub circumcenter {
    my ($O,@x) = @_;
    translate(\@x,$O);
    my @lambdaEqns;
    for (my $i=0; $i<@x; $i++) {
	$lambdaEqns[$i][@x] = innerProd($x[$i],$x[$i]);
	for (my $l=0; $l<@x; $l++) {
	    $lambdaEqns[$i][$l] = 2*innerProd($x[$i],$x[$l]);
	}
    }
    reducedRowEchelon(\@lambdaEqns);
    my @cen;
    for (my $i=0; $i<@x; $i++) {
	my $lambdai = $lambdaEqns[$i][@x] / $lambdaEqns[$i][$i];
	for (my $j=0; $j<$d; $j++) {
	    $cen[$j] += $lambdai * $x[$i][$j];
	}
    }

    ## Translate back to original coordinate system.
    my @minusO = map(-$_, @$O);
    translate(\@x,\@minusO);
    translate([\@cen],\@minusO);
    return \@cen;
}


readC();
$x[$d][$d-1] = 0;    ## establish dimension
for (my $i=0; $i<$d; $i++) { $x[$i][$i] = 1; }
printMatrix(\@x);

my ($spread,$itnum);
my $objective = objectiveValue(\@c,\@x);
do {
    for (my $i=0; $i<$n; $i++) {  ## for each species
	$itnum++;
	printf "iteration %4d    invariant=%6.3f", $itnum, checkDistanceInvariant(\@x);
	## x[i][0], ..., x[i][d-1] are being optimized.
	## First, translate: put vector to optimize at origin.

	## Collect coefficients of quadratic equation for $lambda.
	my ($A,$B,$C) = (0,0,0);
	for (my $j1=0; $j1<$n; $j1++) {
	    next if $j1==$i;
	    for (my $j2=0; $j2<$n; $j2++) {
		next if $j2==$i;
		my $ip = innerProd($x[$j1],$x[$j2]);
		$A -= $ip;
		$C += $c[$i][$j1]*$c[$i][$j2]*$ip;
	    }
	    $A += ($n-1)* (innerProd($x[$j1],$x[$j1]) - distanceSq($x[$i],$x[$j1]))
	}
	my $lambda1 = sqrt(-$C / $A);
	my $lambda2 = -$lambda1;
	
	my (@nu1,@nu2);
	$#nu1 = $#nu2 = $d-1;
	for (my $k=0; $k<$d; $k++) {
	    my ($nu1k, $nu2k);
	    for (my $j=0; $j<$n; $j++) {
		next if $j==$i;
		$nu1k += ($c[$i][$j]+$lambda1)*$x[$j][$k];
		$nu2k += ($c[$i][$j]+$lambda2)*$x[$j][$k];
	    }
	    $nu1[$k] = $nu1k / (($n-1)*$lambda1);
	    $nu2[$k] = $nu2k / (($n-1)*$lambda2);
	}
	my ($objective1, $objective2);
	for (my $j=0; $j<$n; $j++) {
	    next if $j==$i;
	    $objective1 += $c[$i][$j]*distanceSq(\@nu1, $x[$j]);
	    $objective2 += $c[$i][$j]*distanceSq(\@nu2, $x[$j]);
	}
	$x[$i] = ($objective1>$objective2) ? \@nu1 : \@nu2;

	$objective = objectiveValue(\@c,\@x);
	printf "    Objective Value = %10.2f\n", $objective;
#	printMatrix(\@x);
    }
    printMatrix(\@x);

    my ($maxratio, $minratio) = (0,1E99);
    for (my $i=0; $i<$n; $i++) {
	for (my $j=0; $j<$i; $j++) {
	    my $distij = distance($x[$i],$x[$j]);
	    my $ratio = $distij/($dist[$i][$j] || 1);
	    $maxratio = max($maxratio, $ratio);
	    $minratio = min($minratio, $ratio);
	    $dist[$i][$j] = $distij;
	}
    }
    printLowerTriangularMatrix(\@dist);
    $spread = $maxratio/$minratio;
    printf "Objective Value = %f spread = %7.4f\n", $objective, $spread;
} while ($itnum<3*$n || $spread>1.001);

my $maxdist = max(map(max(@{$_||[]}), @dist));
printLowerTriangularMatrix(\@x);
printLowerTriangularMatrix(\@dist);
constructMidpointTree();
constructCentroidTree();
constructCircumcenterTree();



