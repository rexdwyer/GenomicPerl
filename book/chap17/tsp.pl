#!/usr/bin/perl -I . -I /home/dwyer/book/perl

use strict;

my (%graph, %probes, %clones);

sub addWeightedEdge {
    my ($v1, $v2, $weight) = @_;
    push @{$graph{$v1}}, [$v2, $weight];
}

## Read in strings, build suffix trees, create weighted graph.
my $clone;
while (my $in=<STDIN>) {
    $clone = $1 if $in =~ s/^(.*):\s*//;
    chomp($in);
    next unless $in =~ /\S/;   # skip whitespace lines;
    $probes{$clone} ||= [];
    $in =~ s/^\s+//;
    foreach my $probe (split /\s+/, $in) {
	$clones{$probe} ||= [];
	push @{$clones{$probe}}, $clone;
    }
}

print (keys %clones), "\n";
foreach my $p1 (keys %clones) {
    foreach my $p2 (keys %clones) {
	last if $p1 eq $p2;
	my %tmp = map(($_=>2), @{$clones{$p1}});
	my $dist = @{$clones{$p1}} + @{$clones{$p2}};
	foreach my $cl (@{$clones{$p2}}) { $dist -= $tmp{$cl}; }
	addWeightedEdge($p1, $p2, $dist);
	addWeightedEdge($p2, $p1, $dist);
	print "$p1 $p2 $dist\n" if ($p1.$p2) =~ /F/;
    }
}

## Sort adjacency lists to put small weights first.
foreach (keys %graph) {
    my @temp = sort {$$a[0] cmp $$b[0] } @{$graph{$_}};
    $graph{$_} = \@temp;
}

foreach (keys %graph) {
    print $_, ": ",  map($$_[0], @{$graph{$_}}), "\n";
}


my $bestWeight = 9E99;
my $bestString;
my %visited;
my @bestPath;
my $numProbes = scalar(keys %graph);
print "$numProbes\n";
foreach (keys %graph) { search($_, 0, ($_)); }
print "\n\nBest Probe Ordering ($bestWeight) is: @bestPath\n";
$|=1;

sub prospects {
    my ($newV,$viz) = @_;
    my @unviz = grep { !$$viz{$_} } (keys %graph);
    my %minIn;
    foreach my $v (@unviz) {
	foreach my $edge (@{$graph{$v}}) {
	    my ($u, $edgeWeight) = @$edge;
	    last if $edgeWeight == 0;
	    $minIn{$u} = $edgeWeight if $edgeWeight<$minIn{$u};
	}
    }
    my $result = 0;
    foreach my $v (@unviz) { $result += $minIn{$v}; }
    return ($result - $minIn{$newV});
}
	    
sub search {
    my ($v,$pathWeight,@path) = @_;
    print join(',', @_), "\n" if @path<=20 && @path%5==0;
    if (@path == $numProbes) {
	if ($pathWeight<$bestWeight) {
	    ($bestWeight,@bestPath)=($pathWeight,@path);
	    print "New Best ($pathWeight): @path\n";
	}
	return;
    }
    return if $visited{$v};
    return if ($pathWeight + prospects($v,\%visited)) >= $bestWeight;
    $visited{$v} = 1;
    foreach my $edge (@{$graph{$v}}) {
	my ($u, $edgeWeight) = @$edge;
	search($u,$pathWeight+$edgeWeight, (@path, $u));
    }
    $visited{$v} = 0;
}    
