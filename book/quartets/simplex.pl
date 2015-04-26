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
#	my ($l1,$l2,$r1,$r2,$confid) = m/(\d*),(\d*);(\d*),(\d*):(\d*)/;
	my ($l1,$l2,$r1,$r2) = m/(\d*),(\d*);(\d*),(\d*)/;
	my $confid = 1;
	# Next two pairs should be close; penalize (subtract from objective) if far apart.
	$c[$l1][$l2] -= 3*$confid;
	$c[$r1][$r2] -= 3*$confid;
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

sub findVolume {
    my ($x) = @_;
    my @M;
    for (my $j=0; $j<$n; $j++) {  ## for each species
	for (my $k=0; $k<$d; $k++) {  ## for each dimension
	    $M[$j][$k] = $$x[$j][$k];
	}
	$M[$j][$d] = 1;
    }
    reducedRowEchelon(\@M);
    my $vol = 1;
    for (my $j=0; $j<$n; $j++) {  $vol *= $M[$j][$j]; }
    if (abs($vol) > 1.1) {
	print "Volume = $vol\n";
	printMatrix($x);
	die "Volume error";
    }
    return $vol;
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
	    $dist[$nu][$j] = int(0.5 + distance($x[$nu],$x[$j]));
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
	    $dist[$nu][$j] = int(0.5 + distance($x[$nu],$x[$j]));
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
	    $dist[$nu][$j] = int(0.5 + distance($x[$nu],$x[$j]));
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

my ($spread,$itnum);
my $objective = objectiveValue(\@c,\@x);
do {
    for (my $i=0; $i<$n; $i++) {  ## for each species
	$itnum++;
#	printf "iteration %4d    volume=%6.3f", $itnum, findVolume(\@x);
	## x[i][0], ..., x[i][d-1] are being optimized.
	## First, translate: put vector to optimize at origin.
	translate(\@x, $x[$i]);
	## Now find unit normal to hyperplane defined by the other points.
	my @a = findHyperplaneEqn(\@x, $i);

	## Now compute some sums we need.
	my $lambda;
	for (my $k=0; $k<$d; $k++) {
	    for (my $j=0; $j<$n; $j++) {
		$lambda += $c[$i][$j]*$x[$j][$k]*$a[$k];
	    }
	}
	my $sum_ci; 
	for (my $j=0; $j<$n; $j++) { $sum_ci += $c[$i][$j]; }
#	printf "   lambda=%6.2f", $lambda;

	for (my $k=0; $k<$d; $k++) {
	    my $sum_cx_ik;
	    for (my $j=0; $j<$n; $j++) {
		$sum_cx_ik += $c[$i][$j] * $x[$j][$k];
	    }
	    ## We don't move all the way to the optimal point;
	    ## just 1/2 of the distance from the old point.
	    $x[$i][$k] = 0.5 * 
		($x[$i][$k] + (($sum_cx_ik-$a[$k]*$lambda) / $sum_ci));
	}

	$objective = objectiveValue(\@c,\@x);
#	printf "    Objective Value = %10.2f\n", $objective;
#	printMatrix(\@x);
    }
#   printMatrix(\@x);

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
#   printLowerTriangularMatrix(\@dist);
    $spread = $maxratio/$minratio;
    printf "Objective Value = %f spread = %7.4f\n", $objective, $spread;
} while ($itnum<3*$n || $spread>1.001);

my $maxdist = max(map(max(@{$_||[]}), @dist));
printf "\n\nScaling to maxdist=%8.2f:\n", $maxdist;
foreach (@dist) { foreach (@$_) { $_ = int(0.5 + (999*$_/$maxdist));}}
foreach (@x) { foreach (@$_) { $_ *= (999/$maxdist);}}
printLowerTriangularMatrix(\@x);
printLowerTriangularMatrix(\@dist, "%4d");
constructMidpointTree();
constructCentroidTree();
constructCircumcenterTree();



