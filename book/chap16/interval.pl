#!/usr/bin/perl -I . -I /home/dwyer/book/perl
use strict;
use Util;
use IntegerInterval; 

my $error = 0.05;
$|=1;
my @distances;
my %distAvail;
my @solutions;
my @sites;
my @d;
### Turn distances into intervals and set up availability counts.
my %tmp;
foreach (map(split(/\D+/,$_), <STDIN>)) { $tmp{$_}++;};
my @distances = sort {$a<=>$b} keys %tmp;
for (my $i = 0; $i<@distances; $i++) {
    my $d = $distances[$i];
    $distances[$i] = IntegerInterval->new((1-$error)*$d,(1+$error)*$d);
    $distAvail{$distances[$i]->toString()} = $tmp{$d};
}
undef %tmp;

### Set up first and last site.
my $numDists;
foreach (values %distAvail) { $numDists += $_; }
@sites = (("-") x (1+int(sqrt(2*$numDists))));
$sites[0] = IntegerInterval->new(0,0);
$d[0][$#sites] = $d[$#sites][0] = $sites[-1] = $distances[-1];
$distAvail{$distances[-1]->toString()}--;
placeSite(1,$#sites-1);
foreach (@solutions) {
    print "Solution: ";
    foreach (@$_) { $_->print(" "); }
    print "\n";
}
exit();

sub placeSite {
    my ($l,$r) = @_;
    my $i=$#distances;
    while (($i>=0) && !$distAvail{$distances[$i]->toString}) { $i--; }
    if ($i<0) {   ## All distances accounted for.
	my $revisedSites = checkGlobalConstraints();
	push @solutions, $revisedSites if $revisedSites;
	return;
    }
    my $dist = $distances[$i];
    my $newSite = $sites[-1]->minus($dist);
    verifyDistances($l, $r, "L", $newSite);
    verifyDistances($l, $r, "R", $dist);
}

sub verifyDistances {
    verifyDistance(@_,0);
}
    
sub verifyDistance {
    my ($l,$r,$LR,$newSite,$j) = @_;
    ### Trying to place newSite at position l.
    ### Positions 0..l-1 and r+1,last already filled.
    ### If j<l, we have fit distances between newSite
    ### and sites[0],...,sites[j-1].
    ### If j>r, we have fit distances between newSite
    ### and sites[0],...,sites[l-1] and sites[r+1],...,sites[j]

    my $n = ($LR eq "L")? $l : $r;
    if ($j==@sites) {  ### We have found all the distances we need!
	$sites[$n] = $newSite;
	if ($LR eq "L") { placeSite($l+1, $r); }
	else { placeSite($l, $r-1); }
	$sites[$n] = "-";
	return;
    }

    ### Look for a distance that can fit the distance between
    ### sites[j] and newSite.
    my $distSought = 
	($j<$n) ? $newSite->minus($sites[$j]) : $sites[$j]->minus($newSite);
    my $i = @distances;
    while (--$i>=0) {
	next unless $distAvail{$distances[$i]->toString};
	next unless $distances[$i]->intersect($distSought);
	$distAvail{$distances[$i]->toString}--;
	$d[$j][$n] = $d[$n][$j] = $distances[$i];
	my $range = ($j<$n) ? 
	    $sites[$j]->plus($distances[$i]) : $sites[$j]->minus($distances[$i]);
	my $nextJ = ($j==$l-1) ? $r+1 : $j+1;
	verifyDistance($l, $r, $LR, $newSite->intersect($range), $nextJ);
	$distAvail{$distances[$i]->toString}++;
    }
}

sub checkGlobalConstraints {
    my @revisedSites = @sites;
    foreach my $j (1..$#sites) {
	foreach my $i (0..$j-1) {
	    $revisedSites[$j] = 
		$revisedSites[$j]->
		    intersect($revisedSites[$i]->plus($d[$i][$j]));
	    $revisedSites[$j] or return undef;
	}
	foreach my $i (reverse(0..$j-1)) {
	    foreach my $k ($i+1..$j) {
		$revisedSites[$i] = 
		    $revisedSites[$i]->
			intersect($revisedSites[$k]->minus($d[$i][$k]));
		$revisedSites[$i] or return undef;
	    }
	}
    }
    return \@revisedSites;
}
    
