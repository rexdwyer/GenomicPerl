#!/usr/bin/perl -I /home/dwyer/book/perl

use strict;
use Util;

my @M;    ## scoring matrix
my @how;    
my $g = -11;   ## gap penalty
my %blosum;

sub readblosum {
    open BLOSUM, "<blosum62";
    my $header = <BLOSUM>;
    chomp $header;
    my @bases = split /\s+/, $header;
    shift @bases;
    while (<BLOSUM>) {
	my @scores = split;
	my $aa = shift @scores;
        foreach (@bases) { 
	    $blosum{$_,$aa} = shift @scores;
	}
    }
}
	
sub p {
    my ($aa1, $aa2) = @_;
    return $blosum{$aa1,$aa2};
}


sub sim2 {
    my($s,$t) = @_;
    foreach my $i (0..length($s)) { $M[$i][0] = $g * $i; }
    foreach my $j (0..length($t)) { $M[0][$j] = $g * $j; }
    foreach my $i (1..length($s)) {
	foreach my $j (1..length($t)) {
	    my $p =  p(substr($s,$i-1,1),substr($t,$j-1,1));
	    $M[$i][$j] = 
                max($M[$i-1][$j] + $g,
                    $M[$i][$j-1] + $g,
                    $M[$i-1][$j-1] + $p);
	}
    }
    return ( $M[length($s)][length($t)] );
}

sub getAlignment2 {
    my ($s,$t) = @_;
    my ($i,$j) = (length($s), length($t));
    return ( "-"x$j, $t) if ($i==0);
    return ( $s, "-"x$i) if ($j==0);
    my ($sLast,$tLast) = (substr($s,-1),substr($t,-1));
    
    if ($M[$i][$j] == $M[$i-1][$j-1] + p($sLast,$tLast)) { ## Case 1
	print p($sLast,$tLast), " + ";
        ## last letters are paired in the best alignment
	my ($sa, $ta) = getAlignment2(substr($s,0,-1), substr($t,0,-1));
	return ($sa . $sLast , $ta . $tLast );
    } elsif ($M[$i][$j] == $M[$i-1][$j] + $g) { ## Case 2
	print "$g + ";
        ## last letter of the first string is paired with a gap
	my ($sa, $ta) = getAlignment2(substr($s,0,-1), $t);
	return ($sa . $sLast , $ta . "-");
    } else { ## Case 3: last letter of the 2nd string is paired with a gap
	print "$g + ";
	my ($sa, $ta) = getAlignment2($s, substr($t,0,-1));
	return ($sa . "-" , $ta . $tLast );
    }
}



sub similarity {
    my($s1,$s2,$s3) = @_;

    ### Fill in edges of cube.
    foreach my $i1 (0..length($s1)) { $M[$i1][0][0] = $g * $i1 * 2; }
    foreach my $i2 (0..length($s2)) { $M[0][$i2][0] = $g * $i2 * 2; }
    foreach my $i3 (0..length($s3)) { $M[0][0][$i3] = $g * $i3 * 2; }

    ### Fill in sides of cube.
    ## Side 1
    foreach my $i1 (1..length($s1)) {
	my $aa1 = substr($s1,$i1-1,1);
	foreach my $i2 (1..length($s2)) {
	    my $aa2 = substr($s2,$i2-1,1);
	    $M[$i1][$i2][0] = max($M[$i1-1][$i2][0] + $g + $g,
				  $M[$i1][$i2-1][0] + $g + $g,
				  $M[$i1-1][$i2-1][0] + $g + p($aa1,$aa2));
	}
    }

    ## Side 2
    foreach my $i1 (1..length($s1)) {
	my $aa1 = substr($s1,$i1-1,1);
	foreach my $i3 (1..length($s3)) {
	    my $aa3 = substr($s3,$i3-1,1);
	    $M[$i1][0][$i3] = max($M[$i1-1][0][$i3] + $g + $g,
				  $M[$i1][0][$i3-1] + $g + $g,
				  $M[$i1-1][0][$i3-1] + $g + p($aa1,$aa3));
	}
    }
    
    ## Side 3
    foreach my $i2 (1..length($s2)) {
	my $aa2 = substr($s2,$i2-1,1);
	foreach my $i3 (1..length($s3)) {
	    my $aa3 = substr($s3,$i3-1,1);
	    $M[$i2][0][$i3] = max($M[0][$i2-1][$i3] + $g + $g,
				  $M[0][$i2][$i3-1] + $g + $g,
				  $M[0][$i2-1][$i3-1] + $g + p($aa2,$aa3));
	}
    }
    

    ### Fill in interior of cube.
    foreach my $i1 (1..length($s1)) {
	my $aa1 = substr($s1,$i1-1,1);
	foreach my $i2 (1..length($s2)) {
	    my $aa2 = substr($s2,$i2-1,1);
	    my $p12 = p($aa1,$aa2);
	    foreach my $i3 (1..length($s3)) {
		my $aa3 = substr($s3,$i3-1,1);
		my $p13 = p($aa1,$aa3);
		my $p23 = p($aa2,$aa3);
		my @L = ($M[$i1-1][$i2-1][$i3-1] + $p12 + $p13 + $p23,
			    $M[$i1-1][$i2-1][$i3] + $p12 + $g + $g,
			    $M[$i1-1][$i2][$i3-1] + $g + $p13 + $g,
			    $M[$i1][$i2-1][$i3-1] + $g + $g + $p23,
			    $M[$i1][$i2][$i3-1] + 0 + $g + $g,
			    $M[$i1][$i2-1][$i3] + $g + 0 + $g,
			    $M[$i1-1][$i2][$i3] + $g + $g + 0);
		my @I = sort {$L[$b] <=> $L[$a]} (0..6);
		$M[$i1][$i2][$i3] = @L[$I[0]];
		$how[$i1][$i2][$i3] = $I[0];
	    }
	}
    }
    print "$M[length $s1][length $s2][length $s3]\n";
    my ($i1,$i2,$i3) = (length $s1, length $s2, length $s3);
    my ($t1,$t2,$t3);
    while ($i1 * $i2 * $i3 > 0) {
	my $h = $how[$i1][$i2][$i3];
#	print "$i1 $i2 $i3 $h -> ";
	if ($h == 0 || $h == 1 || $h == 2 || $h == 6) {
	    $i1--; $t1 = substr($s1,$i1,1) . $t1;
	} else {  $t1 = "-" . $t1; }
	if ($h == 0 || $h == 1 || $h == 3 || $h == 5) {
	    $i2--; $t2 = substr($s2,$i2,1) . $t2;
	} else {  $t2 = "-" . $t2; }
	if ($h == 0 || $h == 2 || $h == 3 || $h == 4) {
	    $i3--; $t3 = substr($s3,$i3,1) . $t3;
	} else {  $t3 = "-" . $t3; }
#	print "$i1 $i2 $i3 $t1 $t2 $t3\n";
    }

    return ( "$t1\n$t2\n$t3\n" );
}


##MAIN
{
    readblosum();
    my $s1 = <DATA>; chomp($s1);
    my $s2 = <DATA>; chomp($s2);
    my $s3 = <DATA>; chomp($s3);

    print similarity($s1,$s2,$s3), "\n";
    print "\n";

    print "Similarity score: ", sim2($s1,$s2), "\n";
    print "Alignment: \n";
    foreach my $x (getAlignment2($s1,$s2)) {print "\n$x"; }
    print "\n";
    print "\n";

    print "Similarity score: ", sim2($s1,$s3), "\n";
    print "Alignment: \n";
    foreach my $x (getAlignment2($s1,$s3)) {print "\n$x"; }
    print "\n";
    print "\n";

    print "Similarity score: ", sim2($s2,$s3), "\n";
    print "Alignment: \n";
    foreach my $x (getAlignment2($s2,$s3)) {print "\n$x"; }
    print "\n";
    print "\n";
    
}

__END__
TNLWGKVNINELGGEALGRLLVVYPWTQRFFEAFGDLS
FSKIGGHAEEYGAETLERMFIAYPQTKTYFPHFDLS
VQSSWKAVSHNEVEILAAVFAAYPDIQNKFSQFAGK


FSKIGGHAEEYGAETLERMFIAYPQTKTYFPHFDLS


GVLTDVQVALVKSSFEEFNANIPKNTHRFFTLVLEIAPGAKDLFSFLKGSSVHLSGGE
KSAVTNLWGKVNINELGGEALGRLLVVYPWTQRFFEAFGDLS
VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSVLSAADKTN
VKGVFSKIGGHAEEYGAETLERMFIAYPQTKTYFPHFDLS
VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSSPLTADE
ASLVQSSWKAVSHNEVEILAAVFAAYPDIQNKFSQFAGK



