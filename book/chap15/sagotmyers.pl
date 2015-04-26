#!/usr/bin/perl -I . -I /home/dwyer/book/perl

use strict;
use Util;
my %L;
my %D;
my ($maxModelLength, $minModelLength, $minModelRepeats,$maxJump,$epsilon)
    = (15,5,5,2,0.10);
my $e = int($epsilon * $maxModelLength);
my $s =
    "aaaaccccggggttttaaaaccccggggttttaaaaccccggggttttaaaaccccggggtttt" 
    . "aggtcaacct"
    . "agtcaacct"
    . "aggtcaaacct"
    . "aggtcaacct"
    . "aggtcaacct"
    . "aggtcaacct"
    . "aggtcaacct"
    . "agcgaacct"
    . "aggtcaacct"
    . "agtcaacct"
    . "aggtcaagct"
    . "aggtcaactt"
    . "aggatcaacct"
    . "aaaaccccggggttttaaaaccccggggttttaaaaccccggggttttaaaaccccggggtttt";
## the sequence in which satellites are sought.


my (%jump, %gap);
foreach my $x (1..$maxJump) {
    foreach my $y ($minModelLength..$maxModelLength) {
	$jump{$x*$y}++;
	$gap{($x-1)*$y}++;
    }
}
my @jump = sort {$a<=>$b} (keys %jump);
my @gap = sort {$a<=>$b} (keys %gap);
undef (%jump,%gap);


my @L0 = reverse(0..length $s);
my %D0;
foreach my $i (0..length $s) { 
    foreach my $delta (0..$e) {
	$D{$i,$delta} = $delta;
    }
}

my $model0 = ["", \@L0, \%D0 ];
my @models = ( $model0 );

while (@models) {
    my ($model, $L, $D) = @{pop @models};
#   print "Processing Model: /$model/\n";
    next if length $model >= $maxModelLength;
    foreach my $base ('a','c','g','t') {
	my $eModel = $base.$model;
	my (@eL, %eD);
	foreach (@$L) {
	    my $i = $_ - 1;
	    my $iFlag = 0;
	    last if $i<0;
	    my $match = (substr($s,$i,1) eq $base)? 0 : 1;
	    foreach my $delta (-min($e,length $eModel) .. $e) {
		my $eD = $$D{$i+1,$delta} + $match;
		$eD = min($eD,$$D{$i,$delta+1}+1) 
		    if defined($$D{$i,$delta+1});
		$eD = min($eD,$eD{$i+1,$delta-1}+1)
		    if defined($eD{$i+1,$delta-1});
		$eD{$i,$delta} = $eD;
		$iFlag=1 if ($eD <= $e);
	    }
	    push @eL,$i if $iFlag;
	}

	my (%Lcnt, %Rcnt);
	foreach my $i (reverse @eL) { ## ascending order
	    $Lcnt{$i} = 1 + max(map { $Lcnt{$i-$_} } @jump);
#	    print "Lcnt{$i} = $Lcnt{$i}\n";
	}
	
	my @eLFinal;
	foreach my $i (@eL) { ## descending order
	    $Rcnt{$i} = 1 + max(map { $Rcnt{$i+$_} } @jump);
#	    print "Rcnt{$i} = $Rcnt{$i}\n";
	    push @eLFinal,$i if $Lcnt{$i}+$Rcnt{$i} > $minModelRepeats;
	}
	
	if (@eLFinal) {
	    push @models, [$eModel, \@eLFinal, \%eD];
	    record($eModel, \@eLFinal, \%eD)
		if (length $eModel >= $minModelLength);
	} else {
#	    print "Rejecting Model $eModel\n";
	}
    }
}


    my ($I,$J);
    my (%R,%L);
    my (%Right, %Left);
    my (%Rcnt,%Lcnt);
    my (%Rmrk,%Lmrk);
    my %mark;
    my $debug;

sub record {
    my ($model, $L, $D) = @_;
#    print "Recording $model\n";

    undef $I; undef $J;
    undef %R; undef %L;
    undef %Rcnt; undef %Lcnt;
    undef %Rmrk; undef %Lmrk;
    undef %Right; undef %Left;
    undef %mark;
    $debug =  ($model eq "aggtcaacct");

    sub RMark {   ## Declared inside scope of sub record !
	my ($j,$space) = @_;
	print "$space>Rmark($j) ($I,$J), (", join(',',@{$Left{$j}}), 
	      ") (", join(',',(map {$j+$_} @gap)), "\n"
	if $debug;

	$mark{$j}=1;
	$J = max($j,$J);
	foreach my $i (@{$Left{$j}}, (map {$j+$_} @gap)) {
	    LMark($i, " $space") 
		if $L{$i} && (!$mark{$i}) && 
		    ($Lcnt{$i}+$Rcnt{$i}>2*$minModelRepeats);
	}
	print "$space<Rmark($j) $I,$J\n" if $debug;
    }

    sub LMark {   ## Declared inside scope of sub record !
	my ($i,$space) = @_;
	print "$space>Lmark($i) ($I,$J), (", join(',',@{$Right{$i}}), 
	      ") (", join(',',(map {$i-$_} @gap)), "\n"
	if $debug;
	$mark{$i}=1;
	$I = min($i,$I);
	foreach my $j (@{$Right{$i}}, (map {$i-$_} @gap)) {
	    RMark($j, " $space") 
		if $R{$j} && (!$mark{$j}) && 
		    ($Lmrk{$j}+$Rmrk{$j}>2*$minModelRepeats);
	}
	print "$space<Lmark($i) $I,$J\n" if $debug;
    }


    foreach my $l (@$L) {
	$L{$l} = 1;
	foreach my $delta (-$e..$e) {
	    my $r = $l+$delta+(length $model);
	    if (defined($D{$l,$delta})) {
		$R{$r}=1;
		$Right{$l} ||= [];
		push @{$Right{$l}}, $r;
		$Left{$r} ||= [];
		push @{$Left{$r}}, $l;
	    }
	}
    }

    my @LR = sort {$b<=>$a} (@$L, (keys %R));

    foreach my $i (@LR) {
	$Rmrk{$i} = 1 + max(map {$Rcnt{$i+$_}} @gap) if $R{$i};
	$Rcnt{$i} = 1 + max(map {$Rmrk{$_}} @{$Right{$i}})if $L{$i};
    }

    foreach my $i (reverse @LR) {
	$Lcnt{$i} = 1 + max(map{$Lmrk{$i-$_}} @gap) if $L{$i};
	$Lmrk{$i} = 1 + max(map {$Lcnt{$_}} @{$Left{$i}}) if $R{$i};
    }

    if ($debug) {
	foreach my $i (reverse @LR) {
	    print "$i (", join(',',@{$Left{$i}||[]});
	    print ") (", join(',',@{$Right{$i}||[]});
	    print ") mrk=$Lmrk{$i},$Rmrk{$i}";
	    print "; cnt=$Lcnt{$i},$Rcnt{$i}\n";
	}
	foreach (0..20) { print "\n"; }
    }


    foreach my $j (keys %R) {
	print "j=$j; ", 0+$mark{$j},"; ", $Lmrk{$j}+$Rmrk{$j},"\n" if $debug;
	if (!defined($mark{$j}) && ($Lmrk{$j}+$Rmrk{$j}>2*$minModelRepeats)) {
	    $I = $J = $j;
	    RMark($j);
	    print "($I,$J) $model ", substr($s,$I, $J-$I+1), "\n" if $debug;
	}
    }
    die if $debug;
}






