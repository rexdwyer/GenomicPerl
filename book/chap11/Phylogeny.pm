package Phylogeny;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(new allEdges isLeaf neighbors printTable
	     attachLeafBetween removeLeafBetween joinNewNeighbors
	     countMutations);

use strict;

sub new {
    my ($this,$seqs,$names) = @_;
    my @characters;
    for (my $i=0; $i<@$seqs; $i++) {
	my @states = split(//, $$seqs[$i]);
	for (my $j=0; $j<@states; $j++) {
	    $characters[$j][$i] = $states[$j];
	}
    }
    $this =
	bless {characters=>\@characters,
	       names=>($names||[]),
	       numTraits=> scalar @characters,
	       numLeaves => scalar @{$seqs},
	       nextBranch => scalar @$seqs,
	       tree=>[]},
	(ref($this) || $this);
    return $this;
}

sub copy {  ## makes a "deep copy" of the phylogeny.
    my ($this) = @_;
    my @characters = map(\@$_, @{$this->{characters}});
    my @tree = map(\@$_, @{$this->{tree}});
    my $copy =
	bless {characters=> \@characters,
	       names     => \@{$this->{names}},
	       numTraits => $this->{numTraits},
	       numLeaves => $this->{numLeaves},
	       nextBranch=> $this->{nextBranch},
	       tree      => \@tree},
	ref($this);
    return $copy;
}

sub neighbors {   ## returns the list of neighbors of v
    my ($this,$v) = @_;
    return @{${$this->{tree}}[$v]};
}

sub joinNewNeighbors {    ## adds an edge between v1 and v2
    my ($this, $v1, $v2) = @_;
    push @{$this->{tree}[$v1]}, $v2;
    push @{$this->{tree}[$v2]}, $v1;
}

sub removeNeighbors {    ## removes v's neighbor list
    my ($this, $v) = @_;
    my $r = $this->{tree}[$v];
    $this->{tree}[$v] = undef;
    return @$r;
}


sub allEdges {  ## returns a list of all edges in the tree
    my ($this) = @_;
    my $tree = $this->{tree};
    my @result;
    for (my $i=0; $i<@$tree; $i++) {
	next unless defined($$tree[$i]);
	foreach (@{$$tree[$i]}) {
	    push @result, [$i,$_] if $i<$_;
	}
    }
    return @result;
}

sub replaceNeighbor {
    my ($this, $v, $oldNeighbor, $newNeighbor) = @_;
    my $r = ${$this->{tree}}[$v];
    for (my $j=0; $j<@$r; $j++) {
	$$r[$j] = $newNeighbor if $$r[$j] == $oldNeighbor; 
    }
}

sub isLeaf {  ## returns true if v is a leaf in the tree.
    my ($this,$v) = @_;
    return ($v < $this->{numLeaves});
}

sub printTable {  ## prints the phylogeny in a tabular form.
    my ($this) = @_;
    print "Phylogeny:";
#   foreach (keys %$this) { print "\n  $_ = $this->{$_}"; }
    my $tree = $this->{tree};
    my $chars = $this->{characters};
    for (my $i=0; $i<@$tree; $i++) {
	printf("\n%5d: (%2s,%2s,%2s) ", $i, @{$$tree[$i]||[]});
	for (my $j=0; $j<$this->{numTraits}; $j++) { print $$chars[$j][$i]; }
    }
    print "\n\n";
}

## creates a new branch node with neighbors leaf, old1, and old2.
## removes edge between old1 and old2.
sub attachLeafBetween {  
    my ($this, $leaf, $old1, $old2,$changeList) = @_;
    my $newBranch = $this->{nextBranch}++;
    $this->replaceNeighbor($old1,$old2,$newBranch);
    $this->replaceNeighbor($old2,$old1,$newBranch);
    ${$this->{tree}}[$newBranch] = [$leaf,$old1,$old2];
    ${$this->{tree}}[$leaf] = [$newBranch];
    return ($this->updateSeqs($leaf,$changeList));
}


## removes the branch node with neighbors leaf, old1, and old2.
## restores the edge between old1 and old2.
sub detachLeafBetween {
    my ($this, $leaf, $old1, $old2, $changeList) = @_;
    foreach (@$changeList) { ## undo changes to character states of branches
	my ($ref, $state) = @$_;
	$$ref = $state;
    }
    my $oldBranch = ($this->removeNeighbors($leaf))[0];
    $this->removeNeighbors($oldBranch);
    $this->{nextBranch}--;
    $this->replaceNeighbor($old1,$oldBranch,$old2);
    $this->replaceNeighbor($old2,$oldBranch,$old1);
}

## counts the total number of mutations for all characters in the phylogeny.
## normally, we update counts incrementally, but this method can be used for
## a final check.
sub countMutations {	      
    my ($this) = @_;
    my $chars = $this->{characters};
    my $n = $this->{numTraits};
    my $mutations = 0;
    foreach ($this->allEdges()) {
	my ($v1, $v2) = @$_;
	for (my $j=0; $j<$n; $j++) {
	    $mutations++ if $$chars[$j][$v1] ne $$chars[$j][$v2];
	}
    }
    return $mutations;
}

## updates the sequences of branch nodes to account for the addition of
## the leaf newLeaf.

sub updateSeqs {
    my ($this, $newLeaf,$changeList) = @_;
    my $totalDelta = 0;
    my $newBranch = ($this->neighbors($newLeaf))[0];
    my ($leaf, $old1, $old2) = $this->neighbors($newBranch);
    
    foreach my $column (@{$this->{characters}}) {
	$$column[$newBranch] = $$column[$old1];
	my $newState = $$column[$leaf];
	next if	$$column[$newBranch] eq $newState;
	my @bestChanges;
	my $traitDelta = $this->findBestChange($column,$leaf,$newBranch,$newState,
					       \@bestChanges);
	$this->makeBestChange($column,$leaf,$newBranch,$newState,
			      \@bestChanges,$changeList);
	$totalDelta += ($traitDelta+1);
    }
   return $totalDelta;
}

sub makeBestChange {
##
## Carries out the most advantageous change found by findBestChange.
## Constructs a list of the changes made for later reversal during backtrack.
##
    my ($this,       ## [in] ref to current phylogeny
	$column,     ## [in] ref to list of states of current character
	$prev,       ## [in] last tree node considered (parent of curr)
	$curr,       ## [in] current tree node
	$newState,   ## [in] state of the character in the new leaf
	$bestChanges,## [in] list of which nodes to change (if parent changes)
	$changeList  ## [out] list of all changes 
	) = @_;
    my ($next1,$next2) = grep($_ ne $prev, $this->neighbors($curr));
    return unless $$bestChanges[$curr];
    my $traitRef = \$$column[$curr];
    my $currState = $$traitRef;
#    print "push @$changeList, [$traitRef,$currState]\n";
    push @$changeList, [$traitRef,$currState];
    $$traitRef = $newState;
    $this->makeBestChange($column,$curr,$next1,$newState,$bestChanges,$changeList);
    $this->makeBestChange($column,$curr,$next2,$newState,$bestChanges,$changeList);
}

sub findBestChange {
##
## Investigates whether it is advantageous to change branch node $curr's 
## character to $newState, given that neighbor $prev's state is to be changed.
## To determine this, recursive calls to the other two neighbors are necessary.
## Return value is the total reduction in mutation count achievable in the
## subtree rooted at $curr by changing $curr's state to $newState.
##
## In calculating the change in number of mutations, we include the tree
## edge from $curr to $prev as well as all edges in the subtree.
##
    my ($this,       ## [in] ref to current phylogeny
	$column,     ## [in] ref to list of states of current character
	$prev,       ## [in] last tree node considered (parent of curr)
	$curr,       ## [in] current tree node
	$newState,   ## [in] state of the character in the new leaf
	$bestChanges ## [out] list of which nodes to change (if parent changes)
	) = @_;
    
    my $currState = $$column[$curr];
    my $mutBefore = ($currState eq $$column[$prev]) ? 0 : 1;
    $$bestChanges[$curr] = 0;

    ## Option 0:  don't change this vertex's state in subtree.
    my $mutAfter  = ($currState eq $newState) ? 0 : 1;
    my $delta0 = $mutAfter-$mutBefore;
    # changing is not always an option...
    return $delta0 if ($this->isLeaf($curr) || ($currState eq $newState));

    ## Option 1:  Change this vertex's state in subtree.
#   my $mutAfter  = ($newState eq $newState) ? 0 : 1;
    my ($next1,$next2) = grep($_ ne $prev, $this->neighbors($curr));
    my $b1 = $this->findBestChange($column,$curr, $next1, $newState,$bestChanges);
    my $b2 = $this->findBestChange($column,$curr, $next2, $newState,$bestChanges);
    my $delta1 = $b1 + $b2 - $mutBefore;

    ## Return, save best option.
    return $delta0 if $delta0 <= $delta1;

    $$bestChanges[$curr] = 1;
    return $delta1;
}










