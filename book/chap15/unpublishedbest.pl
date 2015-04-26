#!/usr/bin/perl -I . -I /home/dwyer/book/perl

use strict;
use Util;

my ($maxSatLength,$minSatLength, $minSatRepeats,$epsilon)=(15,5,5,0.10);

main();

sub main {
   $|=1;
   my $s = <STDIN>;  ## first line gives options
   eval $s;

   my $s = join('',<STDIN>); ## remainder is sequence
   $s =~ s/\s//sg;
   my $satellites = findSatellites(\$s);
   my @sats = sort { ($$a[1]+$$a[2])<=>($$b[1]+$$b[2]) } @$satellites;
   foreach (@sats) { print "($$_[1],$$_[2])~$$_[3]  $$_[4]x $$_[0]\n"; }
}

sub findSatellites {
    my ($sref) = @_;

    ### Build hash table for empty-string suffix.
    my @suffix0 = ("");
    my $slen = length $$sref;
    foreach my $l (1..$slen) { 
	push @suffix0, $slen-$l+1,0;
    }
    
    ### prime pump with empty suffix.
    my @suffixes = (\@suffix0);
    my $satellites = [];  # ref to list of results.
    
    ### main loop: 
    while (@suffixes) {
	my $suffix = pop @suffixes;
	next if length $$suffix[0] >= $maxSatLength;
#	printf "Expanding /%s/\n", $$suffix[0];
	foreach my $base ('a','c','g','t') {
	    my $eSuffix = suffixCheck($base, $suffix, $sref);
	    if ($eSuffix) {
		push @suffixes, $eSuffix;
		mergeSatellites($satellites, 
				satelliteCheck($eSuffix,$sref));
	    }
	}
    }
    return $satellites;
}

sub suffixCheck {
    my ($base, $suffix, $sref) = @_;
    my @eSuffix = ($base.$$suffix[0]);
    my $eLen = length $eSuffix[0];
    my $maxErr = int($epsilon * $maxSatLength);  # permissible error

    for (my $i=1; $i<@$suffix; $i+=2) {
	my $l = $$suffix[$i]-1;
	last if $l<0;
	my $eq = (substr($$sref,$l,1) eq $base)? 0 : 1;
	my $eErr = $$suffix[$i+1]+$eq;                      ## pair base to base
	$eErr = min($eErr,$$suffix[$i+3]+1) if ($$suffix[$i+2] == $l);  ## gap in sequence
	$eErr = min($eErr,$eSuffix[$#eSuffix]+1) if $eSuffix[$#eSuffix-1] == $l+1;
                                                         ## gap in satellite
	push @eSuffix,$l,$eErr if ($eErr <= $maxErr);
    }

    my @jump = (max($eLen,$minSatLength)-$maxErr .. $maxSatLength+$maxErr);
    
    my (@repeatsBefore, @repeatsAfter);
    $#repeatsBefore = $#repeatsAfter = $eSuffix[1];  ## fast allocate

    foreach (my $i=$#eSuffix-1; $i>0; $i-=2) {
	my $l = $eSuffix[$i];
	$repeatsBefore[$l] = 1 + max(map($repeatsBefore[$l-$_],@jump));
    }

    my @finalESuffix = ($eSuffix[0]);
    foreach (my $i=1; $i<$#eSuffix; $i+=2) {
	my $l = $eSuffix[$i];
	$repeatsAfter[$l] = 1 + max(map($repeatsAfter[$l+$_],@jump));
	push @finalESuffix,$l,$eSuffix[$i+1]
	    if $repeatsBefore[$l]+$repeatsAfter[$l] > $minSatRepeats;
    }

    if (@finalESuffix>1) {
#	print "/$eSuffix[0]/ accepted\n";
	return \@finalESuffix;
    } else { 
#	print "/$eSuffix[0]/ rejected\n";
	return "";
    }
}

sub satelliteCheck {
    my ($sat,$sref) = @_;
    my $satString = $$sat[0];
    my $satLength = length $satString;
    return [] if ($satLength < $minSatLength);

    ## We have a fixed satellite now; error should not depend on $maxSatLength
    ## Clean out any positions that do not meet the more stringent criteria.
    my $maxErr = int($epsilon * $satLength);
    my @jump = ($satLength-$maxErr .. $satLength+$maxErr);
    my (@repeatsBefore,@repeatsAfter);
    $#repeatsAfter = $#repeatsBefore = $$sat[1];  ## fast allocate
    foreach (my $i=$#$sat-1; $i>0; $i-=2) {
	next if $$sat[$i+1] > $maxErr;
	my $l = $$sat[$i];
	$repeatsBefore[$l] = 1 + max(map($repeatsBefore[$l-$_],@jump));
    }
    my @cleanSat = ($$sat[0]);
    foreach (my $i=1; $i<$#$sat; $i+=2) {
	next if $$sat[$i+1] > $maxErr;
	my $l = $$sat[$i];
	$repeatsAfter[$l] = 1 + max(map($repeatsAfter[$l+$_],@jump));
	push @cleanSat,$l,$$sat[$i+1]
	    if $repeatsBefore[$l]+$repeatsAfter[$l] > $minSatRepeats;
    }
    $sat = \@cleanSat;
    return [] if (@$sat == 1);


    ## Now we must look precisely at the individual substrings matching
    ## the satellite at each position and their chaining potential.
    my (@repeatsAfter,@stringEnd,@stringErr,%span);
    $#repeatsAfter = $#stringErr = $#stringEnd = $$sat[1];  ## fast allocate

    for (my $index=1; $index<@$sat; $index+=2) {
	my $i = $$sat[$index];
	next if $$sat[$index+1] > $maxErr;
	## We need to align each substring beginning at position $i
        ## with all of the satellite using dynamic programming.
        ## We can ignore any entries guaranteed to be more than $maxErr.
	my @M_p = (0..$maxErr);
	foreach my $p (1..$satLength) {
	    my $M_p1_q = ($p>$maxErr+1) ? $M_p[$p-$maxErr-1] : $p-1 ;
	    foreach my $q (max(1,$p-$maxErr)..$p+$maxErr) {
		(my $M_p1_q1, $M_p1_q) = ($M_p1_q, $M_p[$q]);
		my $eq = (substr($$sref,$i+$q-1,1) eq substr($satString,$p-1,1))? 0 : 1;
		$M_p[$q] = min($M_p1_q1+$eq, $M_p[$q-1]+1, $M_p1_q+1);
	    }
	    $M_p[$p+$maxErr+1] = $maxErr+1;  ## too big to matter.
	}
	## With the edit distances from the last row of the dynamic programming
	## programming matrix, we can look at precise chaining potential.
	foreach my $j ($i+$satLength-$maxErr ..  $i+$satLength+$maxErr) {
	    next unless $M_p[$j-$i]<=$maxErr;
	    if ($repeatsAfter[$j] >= $repeatsAfter[$i]) {
		$repeatsAfter[$i] = 1+$repeatsAfter[$j];
		$stringEnd[$i] = $stringEnd[$j] || $j;
		$stringErr[$i] = $stringErr[$j] + $M_p[$j-$i];
	    }
	}
	next if ($repeatsAfter[$i] < $minSatRepeats);
	my $end = $stringEnd[$i];
	$span{$end} = $i
	    if !defined($span{$stringEnd[$i]})
		|| ($end-$i-$stringErr[$i] > $end-$span{$end}-$stringErr[$span{$end}]);
	           ## length minus error
    }
#   print "/$satString/ final check\n";

    ## We return a list of the maximal chains for the satellite.
    my $result = 
	[map
	 ([$satString,$span{$_},$_,$stringErr[$span{$_}],$repeatsAfter[$span{$_}]],
	  (keys %span)
	  )];
#   print "$satString passed final check\n" if @$result;
    return $result;
}

sub mergeSatellites {
    my ($satellites, $newOnes) = @_;
    push @$satellites, @$newOnes;
    return;

    foreach my $newOne (@$newOnes) {
	for (my $i=0; $i<@$satellites;  ) {
	    my $oldOne = $$satellites[$i];
	    my $relation = compareSatellites($newOne,$oldOne);
	    if ($relation == 0) { $i++; next; }
	    splice(@$satellites,$i,1); ## remove old
	    if ($relation == 1) { next; }  ## new beats old
	    $newOne = $oldOne; ## old beats new
	    last;
	    }
	unshift(@$satellites,$newOne);
    }
}

sub compareSatellites {
    my ($satA, $satB) = @_;
    my ($stringA, $iA, $jA, $errA, $repA) = @$satA;
    my ($stringB, $iB, $jB, $errB, $repB) = @$satB;
    my $signA = 0<=>(max(0,$iA-$iB)+max(0,$jB-$jA)+$errA-$errB);
    my $signB = 0<=>(max(0,$iB-$iA)+max(0,$jA-$jB)+$errB-$errA);
    return ($signA <=> $signB) if $signA || $signB;
    return ((length $stringB - length $stringA) || 1);  ## favor short if tied
}

sub redundant {
    my ($string) = @_;
    my $s = reverse $string;
    for (my $i=(length $s)-1; $i>0; $i--) {
	return $i if substr($s,$i) ge $s;
    }
    return 0;
}
