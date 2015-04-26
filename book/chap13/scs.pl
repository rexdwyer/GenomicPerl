#!/usr/bin/perl -I . -I /home/dwyer/book/perl

use strict;
use SuffixTree;

my (%tree, %graph);

sub addWeightedEdge {
    my ($v1, $v2, $weight) = @_;
    push @{$graph{$v1}}, [$v2, $weight];
}

## Read in strings, build suffix trees, create weighted graph.
while (my $in=<STDIN>) {
    chomp($in);
    $graph{$in} = [];
    my $inTree = SuffixTree->new($in);
    foreach (keys %tree) {
	addWeightedEdge($in, $_, $inTree->maxSuffix($_));
	addWeightedEdge($_, $in, $tree{$_}->maxSuffix($in));
    }
    $tree{$in} = $inTree;
}

## Throw away suffix trees.
%tree = undef;
## Sort adjacency lists to put big weights first.
foreach (keys %graph) {
    my @temp = sort {$$b[1] <=> $$a[1]} @{$graph{$_}};
    $graph{$_} = \@temp;
}

my $bestWeight = 0;
my $bestString;
my %visited;
foreach (keys %graph) { search($_, 0, $_); }
print "\n\nShortest Common Superstring is:\n$bestString\n";

sub prospects {
    my ($newV,$viz) = @_;
    my @unviz = grep { !$$viz{$_} } (keys %graph);
    my %maxIn;
    foreach my $v (@unviz) {
	foreach my $edge (@{$graph{$v}}) {
	    my ($u, $edgeWeight) = @$edge;
	    last if $edgeWeight == 0;
	    $maxIn{$u} = $edgeWeight if $edgeWeight>$maxIn{$u};
	}
    }
    my $result = 0;
    foreach my $v (@unviz) { $result += $maxIn{$v}; }
    return ($result - $maxIn{$newV});
}
	    
sub search {
    my ($v,$pathWeight,$pathString) = @_;
    return if $visited{$v};
    return if ($pathWeight + prospects($v,\%visited)) <= $bestWeight;
    if ($pathWeight>$bestWeight) { 
	($bestWeight,$bestString) = ($pathWeight,$pathString);
	print "New Best ($pathWeight): $pathString\n";
    }
    $visited{$v} = 1;
    foreach my $edge (@{$graph{$v}}) {
	my ($u, $edgeWeight) = @$edge;
	search($u,$pathWeight+$edgeWeight,$pathString.substr($u,$edgeWeight)); 
    }
    $visited{$v} = 0;
}    
