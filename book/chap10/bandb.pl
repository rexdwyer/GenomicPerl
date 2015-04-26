#!/usr/bin/perl -I /home/dwyer/book/perl
## Branch & Bound multiple alignment

use strict;
use Util;

my $g=-2; ## gap penalty

sub planMultipleAlignment {
    my ($relevant, @ss) = @_;
    my %M;     ## holds best alignment scores.
    my %how;  ## holds best predecessor of entries of %M.

    my $numDirections = (1<<@ss)-1;
    my @pendingLayers = (); foreach (@ss) { push @pendingLayers, []; }
    my $origin = join(',' , ((0) x (scalar @ss)));
    my @goal = map { length($_) } @ss;
    my $goalLayer=0; foreach (@goal) { $goalLayer += $_; }
    my $currentLayer = [$origin];
    
    $M{$origin} = 0;
    my $count = 0;
    for (my $layer=0; $layer<=$goalLayer; $layer++) {
	print "LAYER $layer:";
	while (my $vec = pop @$currentLayer) {  ### something in current layer.
	    my @vec = split(',', $vec);
	    my $score = $M{$vec};
	    
	    my $inRange = 1;
	    foreach (0..$#vec) {
		($inRange=0, last) if $vec[$_] > $goal[$_];
	    }
	    if ($inRange && &$relevant($score, \@vec)) {
#		print "Relevant: $vec $M{$vec}\n";
		foreach my $direction (1..$numDirections) {
		    my @succ = @vec;
		    my @column = (("-") x @succ);
		    my $layerDelta = 0;
		    for (my $i=0; $direction; $i++, $direction>>=1) {
			if ($direction & 1) {
			    $succ[$i]++;
			    $column[$i] = substr($ss[$i],$succ[$i]-1,1);
			    $layerDelta++;
			}
		    }
		    my $succ = join(',',@succ);
		    my $nuscore = $score + scoreColumn(@column);
		    if (!defined($M{$succ})) {
			($M{$succ},$how{$succ}) = ($nuscore, $vec);
#			print "push $succ onto list $layerDelta\n";
			push @{$pendingLayers[$layerDelta-1]}, $succ;
		    } elsif ($nuscore >= $M{$succ}) {
			($M{$succ},$how{$succ}) = ($nuscore, $vec);
		    }
		}
	    } else { 
#		print "Not relevant: $vec $M{$vec} $how{$vec}\n";
		delete $M{$vec};
	        delete $how{$vec};
	    }
	}
	push @pendingLayers, [];
	$currentLayer = shift @pendingLayers;
    }
    return ($M{join(',',@goal)}, \%how);
}

sub scoreColumn {
    my @col = @_;
    my ($gaps,$aas,$score) = (0,0,0);
    while (@col) {
	my $aa = shift @col;
	($gaps++, next) if $aa eq "-";
	$aas++;
	foreach my $aa1 (@col) {
	    next if $aa1 eq "-";
#	    $score += $blosum62{$aa,$aa1};
	    $score += ($aa eq $aa1) ? +1 : -1;
	}
    }
    return $score + ($g * $gaps * $aas);
}

sub prependColumnToAlignment {
    my ($A, @col) = @_;
    foreach (@$A) {$_ = (shift @col).$_};
}


sub reconstructAlignment {
    my ($how,@ss) = @_;
    my @result = (("") x @ss);
    my $origin = join(',', map {"0"} @ss);
    my $current = join(',', map {length($_)} @ss);
    my @gaps = (("-")x@ss);
    while ($current ne $origin) {
	my @current = split(',',$current);
	my $previous = $$how{$current};
	my @previous = split(',',$previous);
	my @column = @gaps;
	for (my $i=0; $i<@ss; $i++) {
	    if ($current[$i] != $previous[$i]) {
		$column[$i] = substr($ss[$i],$previous[$i],1);
	    }
	}
	prependColumnToAlignment(\@result, @column);
	$current = $previous;
    }
    return \@result;
}


sub tunnelAlign {
    my ($width,@ss) = @_;
    my @goal = map { length($_) } @ss;
    my $goalSq = 0;
    foreach (@goal) { $goalSq += $_*$_; }

    sub tunnelRelevant {
	my ($score,$vecref) = @_;
	my ($i,$veclen2, $proj, $sslen2) = (0,0,0,0);
	foreach $i (0..@ss-1) {
	    $veclen2 += $$vecref[$i]*$$vecref[$i];
	    $proj += $$vecref[$i]*$goal[$i];
	}
	return (($veclen2 - ($proj*$proj / $goalSq)) < $width*$width);
    }

    my ($score, $how) = planMultipleAlignment(\&tunnelRelevant, @ss);
    return ($score,reconstructAlignment($how,@ss));
}

sub computeSuffixPairSimilarities {
    my(@ss) = @_;
    my @c;
    foreach my $p (0..$#ss) {
	my $s = $ss[$p];
	my $m = length $s;
	foreach my $q ($p+1..$#ss) {
	    my $t = $ss[$q];
	    my $n = length $t;
	    $c[$p][$q][$m][$n] = 0;
	    for (my $i=$m; $i>=0; $i--) { $c[$p][$q][$i][$n] = $g * ($m-$i); }
	    for (my $j=$n; $j>=0; $j--) { $c[$p][$q][$m][$j] = $g * ($n-$j); }
	    for (my $i=$m-1; $i>=0; $i--) {
		for (my $j=$n-1; $j>=0; $j--) {
		    my $match = scoreColumn(substr($s,$i,1),substr($t,$j,1));
		    $c[$p][$q][$i][$j] = max($c[$p][$q][$i+1][$j] + $g,
					     $c[$p][$q][$i][$j+1] + $g,
					     $c[$p][$q][$i+1][$j+1] + $match);
		}
	    }
	}
    }
    return \@c;
}

sub branchBoundAlign {
    my ($lowerBound, @ss) = @_;

    my $C = computeSuffixPairSimilarities(@ss);

    sub branchBoundRelevant {
	my ($score,$vecref) = @_;
	my $K = @ss;

	my $upperBound = $score;
	foreach my $p (0..$#ss) {
	    foreach my $q ($p+1..$#ss) {
		$upperBound += $$C[$p][$q][$$vecref[$p]][$$vecref[$q]];
	    }
	}
#	print join(',',@$vecref), " return($lowerBound <= $upperBound)\n";
	return ($lowerBound <= $upperBound);
    }

    my ($score, $how) = planMultipleAlignment(\&branchBoundRelevant, @ss);
    return ($score,reconstructAlignment($how,@ss));
}

sub scoreMultipleAlignment {
    my ($alignment) = @_;
    my $score;
    foreach my $i (0..length($$alignment[0])-1) {
	$score += scoreColumn(map {substr($_,$i,1)} @$alignment);
    }
    return $score;
}

sub printAlignment {
    my ($title, $alignment) = @_;
    print "\n***** $title\n";
    foreach (@$alignment) { print "$_\n"; }
    my $score = scoreMultipleAlignment($alignment);
    print "score $score\n\n";
}    

##MAIN
my @ss;
while (my $s = <DATA>) {
    chomp($s);
    push @ss, $s if $s;  ## tolerate blank lines
}

my ($tunnelScore, $tunnelAlignment) = tunnelAlign(2,@ss);
printAlignment("Tunnel Strategy", $tunnelAlignment);

my ($lcScore, $lcAlignment) = branchBoundAlign($tunnelScore, @ss);
printAlignment("Branch & Bound Strategy", $lcAlignment);
__END__
GVLTDVQVALVKSSFEEFNANIPKNTHRFFTLVLEIAPGAKDLFSFLKGSS
SPLTADEASLVQSSWKAVSHNEVEILAAVFAAYPDIQNKFSQFAGK
VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLS
VHLSGGEKSAVTNLWGKVNINELGGEALGRLLVVYPWTQRFFEAFGDLS
VLSAADKTNVKGVFSKIGGHAEEYGAETLERMFIAYPQTKTYFPHFDLS
