#!/usr/bin/perl -I /home/dwyer/book/perl
## Carillo-Lipman multiple alignment

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
    my $goal = join(',' , map { length($_) } @ss );
    my $currentLayer = [$origin];
    
    $M{$origin} = 0;
    my $layer = 0;
    my $count = 0;
    while (@$currentLayer) {  ### something in current layer.
	print "LAYER $layer:";
	my $relevantCount = 0;
	while (my $vec = pop @$currentLayer) {  ### something in current layer.
#	    print "$vec $M{$vec}\n";
	    my @vec = split(',', $vec);
	    my $score = $M{$vec};
	    if (inRange(\@vec, \@ss) && &$relevant($score, \@vec)) {
		$relevantCount++;
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
		delete $M{$vec};
	        delete $how{$vec};
	    }
	}
	print " $relevantCount relevant entries\n";
#	print "summary of pending layers: ";
#	map { print scalar @$_, " "; } @pendingLayers;
#	print "\n";
        $layer++;
	push @pendingLayers, [];
	$currentLayer = shift @pendingLayers;
    }
    return ($M{$goal}, \%how);
}

sub scoreColumn {
    my @col = @_;
    my $col = join(',',@col);
    my $score = 0;
    while (@col>1) {
	my $aa = shift @col;
	foreach my $aa1 (@col) {
#           $score += $blosum62{$aa,$aa1};
	    if ($aa eq $aa1) {
		$score += ($aa eq "-") ? 0 : 1;
	    } else {
		$score += ($aa eq "-" || $aa1 eq "-") ? $g : -1;
	    }
	}
    }
    return $score;
}

sub inRange {
    my ($vecref, $ssref) = @_;
    foreach (0..@$vecref-1) {
	return 0 if ($$vecref[$_] > length($$ssref[$_]));
    }
    1;
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

    sub tunnelRelevant {
	my ($score,$vecref) = @_;
	my ($i,$veclen2, $proj, $sslen2) = (0,0,0,0);
	foreach $i (0..@ss-1) {
	    $veclen2 += $$vecref[$i]*$$vecref[$i];
	    $proj += $$vecref[$i]*length($ss[$i]);
	    $sslen2 += length($ss[$i])*length($ss[$i]);
	}
	return (($veclen2 - ($proj*$proj / $sslen2)) < $width*$width);
    }

    my ($score, $how) = planMultipleAlignment(\&tunnelRelevant, @ss);
    print "$score\n";
    return ($score,reconstructAlignment($how,@ss));
}

sub nw {  ## Needleman-Wunsch for two strings.
    my($s,$t) = @_;
    my @a;
    my $m = length $s;
    my $n = length $t;
    foreach my $i (0..$m) { $a[$i][0] = $g * $i; }
    foreach my $j (0..$n) { $a[0][$j] = $g * $j; }
    foreach my $i (1..$m) {
	foreach my $j (1..$n) {
            my $match = (substr($s,$i-1,1) eq substr($t,$j-1,1))?1:-1;
	    $a[$i][$j] = max($a[$i-1][$j] + $g,
			     $a[$i][$j-1] + $g,
			     $a[$i-1][$j-1] + $match);
#           print "a[$i][$j] = $a[$i][$j]\n";
	}
    }
    return \@a;
}

sub fillC {
    my(@ss) = @_;
    my $K = @ss;
    my @c;
    my ($nwa1, $nwa2);
    my ($i,$j,$p,$q);
    foreach $i (1..$K-1) {
	my $s = $ss[$i-1];
	my $m = length $s;
	foreach $j ($i+1..$K) {
	    my $t = $ss[$j-1];
	    my $n = length $t;
	    $nwa1 = nw($s,$t);
	    $nwa2 = nw(scalar reverse $s, scalar reverse $t);
            foreach $p (0..$m) {
		foreach $q (0..$n) {
		    $c[$i][$j][$p][$q] = $$nwa1[$p][$q] + $$nwa2[$m-$p][$n-$q];
		}
	    }
	}
    }
    return \@c;
}

sub lipmanCarilloAlign {
    my ($lowerBound, @ss) = @_;

    my $C = fillC(@ss);

    sub carilloRelevant {
	my ($score,$vecref) = @_;
	my $K = @ss;
	my $upperBound = 0;
	foreach my $i (1..$K-1) {
	    foreach my $j ($i+1..$K) {
		$upperBound += $$C[$i][$j][$$vecref[$i-1]][$$vecref[$j-1]];
	    }
	}
	return ($lowerBound <= $upperBound);
    }

    my ($score, $how) = planMultipleAlignment(\&carilloRelevant, @ss);
    print "$score\n";
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
    push @ss, $s;
    print "$s\n";
}

$|=1;  ## don't buffer output

my ($tunnelScore, $tunnelAlignment) = tunnelAlign(2,@ss);
printAlignment("Tunnel Strategy", $tunnelAlignment);

my ($lcScore, $lcAlignment) = lipmanCarilloAlign($tunnelScore, @ss);
printAlignment("Lipman & Carillo's Strategy", $lcAlignment);
__END__
GVLTDVQVALVKSSFEEFNANIPKNTHRFFTLVLEIAPGAKDLFSFLKGSS
SPLTADEASLVQSSWKAVSHNEVEILAAVFAAYPDIQNKFSQFAGK
VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLS
VHLSGGEKSAVTNLWGKVNINELGGEALGRLLVVYPWTQRFFEAFGDLS
VLSAADKTNVKGVFSKIGGHAEEYGAETLERMFIAYPQTKTYFPHFDLS
