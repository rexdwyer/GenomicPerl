#!/usr/bin/perl -I . -I /home/dwyer/book/perl
use strict;
use Util;

my %distAvail;
my @solutions;
foreach (map(split(/\D+/,$_), <STDIN>)) { $distAvail{$_}++;};
my @distances = sort {$a <=> $b} keys %distAvail;
my $numDists;
foreach (values %distAvail) { $numDists += $_; }
my @sites = (("-") x (1+int(sqrt(2*$numDists))));
$sites[0] = 0;
$sites[-1] = $distances[-1];
$distAvail{$distances[-1]}--;
placeSite(1,$#sites-1);
foreach (@solutions) { print "Solution: @$_\n"; }
exit();

sub placeSite {
    my ($l,$r) = @_;
    my $i=$#distances;
    while (($i>=0) && !$distAvail{$distances[$i]}) { $i--; }
    if ($i<0) {   ## All distances accounted for.
	push @solutions, [@sites];  ## Makes copy of @sites.
	return;
    }
    my $dist = $distances[$i];
    my $newSite = $sites[-1]-$dist;
    verifyDistances($l, $r, "L", $newSite);
    verifyDistances($l, $r, "R", $dist);
}

sub verifyDistances {
    my ($l,$r,$LR,$newSite) = @_;

    my $n = ($LR eq "L")? $l : $r;
    my $allAvail=1;
    foreach my $j ((0..$l-1), ($r+1..$#sites)) {
	$allAvail = $distAvail{abs($newSite-$sites[$j])}-- && $allAvail;
    }
    if ($allAvail) {  ### We have found all the distances we need!
	$sites[$n] = $newSite;
	if ($LR eq "L") { placeSite($l+1, $r); }
	else { placeSite($l, $r-1); }
	$sites[$n] = "-";
    }
    foreach my $j ((0..$l-1), ($r+1..$#sites)) {
	$distAvail{abs($newSite-$sites[$j])}++;
    }
}
