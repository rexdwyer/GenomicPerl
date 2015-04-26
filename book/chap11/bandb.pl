#!/usr/bin/perl -I ../perl
use strict; 
use Util;
use Phylogeny;

my $count;

main();

sub main {
    $|=1;
    my @S = <STDIN>;    ## reads all lines at once.
    foreach (@S) { chomp($_); }
    my $greedy = build_tree_greedily(\@S);
    $greedy->printTable();
    my $best = build_trees_exhaustively(\@S);
    $best->printTable();
    $best = build_trees_efficiently(\@S);
    $best->printTable();
    $best = build_trees_efficiently(\@S, $greedy->countMutations());
    $best->printTable();
}


sub build_trees_exhaustively {
    my ($sequences) = @_;
    my $phyl = Phylogeny->new($sequences);
    $phyl->joinNewNeighbors(0,1);
    my $bound = 1E99;
    $count=0;
    my $best;
    build_all_trees($phyl,2, $phyl->countMutations(),\$bound,\$best);
    warn "$count calls to build_all_trees.\n";
    return $best;
}

sub build_all_trees {
    my ($phyl,$newLeaf,$mutations,$bound,$best) = @_;
    $count++;
    if (!($phyl->isLeaf($newLeaf))) {
	if ($mutations<=$$bound) {
	    $$bound = $mutations;
	    warn "Call no. $count: Found tree with $mutations mutations.\n";
	    $$best = $phyl->copy();
	}
	return;
    }
    foreach ($phyl->allEdges()) {
	(my $v1, my $v2) = @$_;
	my $changeList = [];
	my $deltaMut = $phyl->attachLeafBetween($newLeaf,$v1,$v2,$changeList);
	build_all_trees($phyl,$newLeaf+1,$mutations+$deltaMut,$bound,$best);
	$phyl->detachLeafBetween($newLeaf, $v1, $v2, $changeList);
    }
}
    
sub build_tree_greedily {
    my ($sequences) = @_;
    my $phyl = Phylogeny->new($sequences);
    $phyl->joinNewNeighbors(0,1);
    my $mutations = $phyl->countMutations();
    for (my $newLeaf=2; $phyl->isLeaf($newLeaf); $newLeaf++) {
	my ($bestDelta, $bestV1, $bestV2) = (1E99, -1, -1);
	foreach ($phyl->allEdges()) {
	    (my $v1, my $v2) = @$_;

	    my $changeList = [];
	    my $deltaMutations
		= $phyl->attachLeafBetween($newLeaf, $v1, $v2, $changeList);
	    ($bestDelta, $bestV1, $bestV2) = ($deltaMutations,$v1,$v2)
		if $deltaMutations < $bestDelta;
	    $phyl->detachLeafBetween($newLeaf, $v1, $v2, $changeList);
	}
	$phyl->attachLeafBetween($newLeaf, $bestV1, $bestV2, []);
	$mutations += $bestDelta;
    }
    return $phyl;
}

sub build_trees_efficiently {
    my ($sequences,$upperBound) = @_;
    my $phyl = Phylogeny->new($sequences);
    $phyl->joinNewNeighbors(0,1);
    $upperBound ||= 1E99;  ## assume the absurd if no bound given.
    $count=0;
    my $best;
    build_good_trees($phyl,2, $phyl->countMutations(), \$upperBound, \$best);
    warn "$count calls to build_good_trees.\n";
    return $best;
}

sub build_good_trees {
    my ($phyl,$newLeaf,$mutations,$bound,$best) = @_;
    $count++;
    if ($mutations > $$bound) {
#	warn "Pruning at Leaf $newLeaf...\n";
	return;
    }
    if (!($phyl->isLeaf($newLeaf))) {
	warn "Call no. $count: Found tree with $mutations mutations.\n";
	$$bound = $mutations;
	$$best = $phyl->copy();
	return;
    }
    foreach ($phyl->allEdges()) {
	(my $v1, my $v2) = @$_;
	my $changeList = [];
	my $deltaMut = $phyl->attachLeafBetween($newLeaf,$v1,$v2,$changeList);
	build_good_trees($phyl,$newLeaf+1,$mutations+$deltaMut,$bound,$best);
	$phyl->detachLeafBetween($newLeaf, $v1, $v2, $changeList);
    }
}



