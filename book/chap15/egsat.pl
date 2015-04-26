#!/usr/bin/perl -I . -I /home/dwyer/book/perl

use strict;
use Util;

my ($maxSatLength,$minSatLength, $minSatRepeats,$epsilon)=(15,5,5,0.10);
my $trace;
my $traceString =  "gtcaa";

main();

sub main {
   $|=1;
   my $s = <STDIN>;  ## first line gives options
   eval $s;          ## sets options maxSatLength, etc.
   my $s = join('',<STDIN>); ## remainder of file is sequence
   $s =~ s/\s//sg;   ## remove all whitespace
   my $satellites = findSatellites(\$s);
   my @sats = sort { ($$a[1]+$$a[2])<=>($$b[1]+$$b[2]) } @$satellites;
   foreach (@sats) { print "($$_[1],$$_[2])~$$_[3]  $$_[4]x $$_[0]\n"; }
}

sub findSatellites {
    my ($sref) = @_;

    ### Build match set for empty-string suffix & use it to prime pump.
    my @suffix0 = ("");
    my $slen = length $$sref;
    foreach my $i (0..$slen-1) { push @suffix0, $slen-$i,0; }
    my @suffixes = (\@suffix0);
    my @satellites;
    
    ### Main loop:  
    while (@suffixes) {
	my $suffix = pop @suffixes;
	next if length $$suffix[0] >= $maxSatLength;
#	printf "Expanding /%s/\n", $$suffix[0];
	foreach my $base ('a','c','g','t') {
	    my $eSuffix = suffixCheck($base, $suffix, $sref);
	    if ($eSuffix) {
		push @suffixes, $eSuffix;
		push @satellites, @{satelliteCheck($eSuffix,$sref)};
	    }
	}
    }
    return \@satellites;
}

sub print_suffix_matches {
    my ($suffix) = @_;
    return unless $trace;
    my @hex = qw(0 1 2 3 4 5 6 7 8 9 A B C D E F G H I J K);
    my @s = ((" ") x ($$suffix[1]+5));
    for (my $p=1; $p<@$suffix; $p+=2) {
	$s[$$suffix[$p]] = $hex[$$suffix[$p+1] % 10];
    }
    print @s, $$suffix[0], "\n";
}

sub print_reps {
    my ($reps) = @_;
    return unless $trace;
    my @hex = (" ", qw(1 2 3 4 5 6 7 8 9 A B C D E F G H I J K));
    print map($hex[$_], @$reps), "\n";
}


sub suffixCheck {
    my ($base, $suffix, $sref) = @_;
    my @eSuffix = ($base.$$suffix[0]);  ## string for extended suffix
    my $eLen = length $eSuffix[0];
    my $maxErr = int($epsilon * $maxSatLength);  ## permissible error
    $trace = ($traceString =~ m/$eSuffix[0]$/);
    print_suffix_matches($suffix);

    for (my $p=1; $p<@$suffix; $p+=2) {
	my $i = $$suffix[$p]-1;
	last if $i<0;
	my $eq = (substr($$sref,$i,1) eq $base)? 0 : 1;
	my $eErr = $$suffix[$p+1]+$eq;   ## pair base to base
	$eErr = $$suffix[$p+3]+1         ## put gap in sequence 
	    if ($$suffix[$p+2] == $i) && ($$suffix[$p+3]+1 < $eErr);
	$eErr = $eSuffix[$#eSuffix]+1    ## put gap in satellite 
	    if ($eSuffix[$#eSuffix-1] == $i+1) && ($eSuffix[$#eSuffix]+1<$eErr);
	push @eSuffix,$i,$eErr if ($eErr <= $maxErr);
    }
    print_suffix_matches(\@eSuffix);

    my @jump = (max($eLen,$minSatLength)-$maxErr .. $maxSatLength+$maxErr);
    my $filteredESuffix = filterForRepeats(\@eSuffix, $maxErr, \@jump);
#   print "/$eSuffix[0]/ ",( $filteredESuffix ? "accepted\n" : "rejected\n");
    return $filteredESuffix;
}


sub filterForRepeats {
    my ($pattern, $maxErr, $jumps) = @_;
    my (@repeatsBefore, @repeatsAfter);
    $#repeatsBefore = $#repeatsAfter = $$pattern[1];  ## fast allocate
    for (my $p=$#$pattern-1; $p>0; $p-=2) {
	next if $$pattern[$p+1] > $maxErr;
	my $i = $$pattern[$p];
	$repeatsBefore[$i] = 1 + max(map($repeatsBefore[$i-$_],@$jumps));
    }
    my @filteredPattern = ($$pattern[0]);
    for (my $p=1; $p<$#$pattern; $p+=2) {
	next if $$pattern[$p+1] > $maxErr;
	my $i = $$pattern[$p];
	$repeatsAfter[$i] = 1 + max(map($repeatsAfter[$i+$_],@$jumps));
	push @filteredPattern,$i,$$pattern[$p+1]
	    if $repeatsBefore[$i]+$repeatsAfter[$i] > $minSatRepeats;
    }
    print "before:\n";
    print_reps(\@repeatsBefore);
    print "after:\n";
    print_reps(\@repeatsAfter);

    return $#filteredPattern ? \@filteredPattern : "";
}

sub satelliteCheck {
    my ($sat,$sref) = @_;
    my $satString = $$sat[0];
    my $satLength = length $satString;
    return [ ] if ($satLength < $minSatLength);

    ## We have a fixed satellite now; error should not depend on $maxSatLength
    ## Clean out any positions that do not meet the more stringent criteria.
    my $maxErr = int($epsilon * $satLength);
    my @jump = ($satLength-$maxErr .. $satLength+$maxErr);
    print "/$satString/ final check\n";
    $sat = filterForRepeats($sat, $maxErr, \@jump);
    return [ ] unless $sat;

    ## Now we must look precisely at the individual substrings matching
    ## the satellite at each position and their chaining potential.
    my (@repeatsAfter,@stringEnd,@stringErr,%span);
    $#repeatsAfter = $#stringErr = $#stringEnd = $$sat[1];  ## fast allocate

    for (my $p=1; $p<@$sat; $p+=2) {
	my $i = $$sat[$p];
	my $goodJ = goodEndings($satString, $sref, $i, $maxErr);
	for (my $q=0; $q<@$goodJ; $q+=2) {
	    my $j = $$goodJ[$q];
	    if ($repeatsAfter[$j] >= $repeatsAfter[$i]) {
		$repeatsAfter[$i] = 1+$repeatsAfter[$j];
		$stringEnd[$i] = $stringEnd[$j] || $j;
		$stringErr[$i] = $stringErr[$j] + $$goodJ[$q+1];
	    }
	}
	next if ($repeatsAfter[$i] < $minSatRepeats);
	my $end = $stringEnd[$i];
	$span{$end} = $i
	    if !defined($span{$stringEnd[$i]})
		|| ($end-$i-$stringErr[$i] ## length minus error
		    > $end-$span{$end}-$stringErr[$span{$end}]);
    }
    print "/$satString/ final check\n";
    print_reps(\@repeatsAfter);
    print "/$satString/ final check\n";

    ## We return a list of the maximal chains for the satellite.
    return [map([$satString,$span{$_},$_,
		 $stringErr[$span{$_}],$repeatsAfter[$span{$_}]],
		(keys %span)
		)];
}

sub goodEndings {
    my ($satString, $sref, $i, $maxErr) = @_;
    my $satLength = length $satString;
    my @M;
    $M[0] = [0..$maxErr+1];
    foreach my $p (1..$satLength) {
	$M[$p][$p+$maxErr+1] = $maxErr+1; ## too big to matter; allocates row.
	my $lo = max(1,$p-$maxErr);
	$M[$p][$lo-1] = $p-($lo-1);
	foreach my $q ($lo..$p+$maxErr) {
	    if (substr($$sref,$i+$q-1,1) eq substr($satString,$p-1,1)) {
		$M[$p][$q] = $M[$p-1][$q-1];
	    } else { 
		$M[$p][$q] = 1 + min($M[$p-1][$q-1],$M[$p][$q-1],$M[$p-1][$q]);
	    }
	}
    }
    my @good;
    for (my $q=$satLength+$maxErr; $q>=$satLength-$maxErr; $q--) {
	push @good, $i+$q, $M[$satLength][$q] if $M[$satLength][$q]<=$maxErr;
    }
    return \@good;
}



