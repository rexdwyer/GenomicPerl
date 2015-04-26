#!/usr/bin/perl -I /home/dwyer/book/perl

use strict;

sub findComponents {
    my ($mustvisit, $deleteProbe, $clones, $probes) = @_;
    my %visited;
    foreach my $pClone (@{$$clones{$deleteProbe}}) {
	foreach ($pClone,@{$$probes{$pClone}}) { $visited{$_}++; }
    }
    my @listOfComponents;
    foreach (@$mustvisit) {
	my @queue = ($_);
	my @component;
	while (@queue) {
	    my $current = shift @queue;
	    next if $visited{$current}++; 
	    push @component, $current if $$clones{$current};
	    push @queue, @{$$probes{$current} || $$clones{$current}};
	}
	push @listOfComponents, \@component if @component;
    }
    return @listOfComponents;
}

sub findSplitters {
    my ($component,$clones,$probes) = @_;
    my $bestdiff = 0.5;
    my $bestsplitter = -1;
    my @splitters = ();
#   print "\n\n\n";
    foreach my $p (@$component) {
	my @pcomponents = findComponents($component, $p, $clones, $probes);
#	print "probe $p components are (";
#	print join(';', map { join(',',@$_) } @pcomponents);
	if (@pcomponents == 2) {
	    push @splitters, $p;
	    my $size0 = @{$pcomponents[0]};
	    my $size1 = @{$pcomponents[1]};
	    my $diff = abs(0.5 - ($size0 / ($size0+$size1)));
#	    print ")  splitter\n";
#	    print "splitter $p ($size0,$size1) $diff\n";
	    ($bestsplitter, $bestdiff) = ($p,$diff) if ($diff <= $bestdiff);
	} else {
#	    print ")  non-splitter\n";
	}
    }
#   print "\n\nbestsplitter = $bestsplitter\n";
    return ($bestsplitter, \@splitters);
}

sub moreLikeWhich {
    my ($S, $L, $R) = @_;
    my %hash;
    my $bias = 0;
    foreach (@$S) { $hash{$_}=1; }
    foreach (@$L) { $bias -= $hash{$_}; }
    foreach (@$R) { $bias += $hash{$_}; }
    return ($bias <=> 0);
}

sub computeProbeBiases {
    my ($component,$clones,$probes,$bestSplitter,$splitters) = @_;
    my %bias;
    my ($L,$R) = findComponents($component, $bestSplitter, $clones, $probes);
    foreach my $p (@$splitters) {
	my @pcomponents = findComponents($component, $p, $clones, $probes);
	foreach my $S (@pcomponents) {
	    my $sign = moreLikeWhich($S,$L,$R);
	    foreach (@$S) { $bias{$_} += $sign; };
	}
    }
    return %bias;
}


sub groupByHashValue {
    my ($set, $hash) = @_;
    my @orderedSet = sort { $$hash{$a} <=> $$hash{$b} } @$set;
    my @groups;
    my $lastHashValue =  -1E99;
    foreach my $elem  (@orderedSet) {
	push @groups, [] if $$hash{$elem} > $lastHashValue;
	push @{$groups[$#groups]}, $elem;
	$lastHashValue = $$hash{$elem};
    }
    return \@groups;
}

sub orderComponent {
    my ($component, $clones, $probes) = @_;
    
    my ($bestSplitter,$splitters) = findSplitters($component,$clones,$probes);
    my %bias = 
       computeProbeBiases($component,$clones,$probes,$bestSplitter,$splitters);
    return (groupByHashValue($component, \%bias));
}

sub createNewClone {
    my ($oldclone, $suffix, $newProbes, $clones, $probes) = @_;
    my $newclone = "$oldclone~$suffix";
    print("Splitting off $newclone = (", join(',',@$newProbes),")\n");
    $$probes{$newclone} = $newProbes;
    foreach my $p (@$newProbes) {
	my @cList = @{$$clones{$p}};
	for (my $j=0; $j<@cList; $j++) {
	    ($cList[$j] = $newclone, last) if $cList[$j] eq $oldclone;
	}
	$$clones{$p} = \@cList;
    }
}

sub buildGroupIndex {
    my ($groupL) = @_;
    my %index;
    foreach my $i (0..$#$groupL) {
	foreach my $elem (@{$$groupL[$i]}) {
	    $index{$elem} = $i;
	}
    }
    return \%index;
}

sub splitChimerae {
    my ($groupL, $component, $clones, $probes) = @_;
    my $result = 0;
    my $groupH = buildGroupIndex($groupL);
    my %tested;
    foreach my $probe (@$component) {
	foreach my $clone (@{$$clones{$probe}}) {
	    next if $tested{$clone}++;
	    my @cProbes=sort {$$groupH{$a}<=>$$groupH{$b}} @{$$probes{$clone}};
	    my $gLo = $$groupH{$cProbes[0]}+1;
	    my $gHi = $$groupH{$cProbes[@cProbes-1]}-1;
	    next if $gLo>$gHi;  ## clones from 1 or 2 adjacent groupL
	    my %cloneHash;
	    foreach my $p (@{$$probes{$clone}}) { $cloneHash{$p}=1; }

	    my $suffix=1;
	    my $chimerical=0;
	    foreach my $g ($gLo..$gHi) {
		my $gBreak = 0;
		foreach my $p (@{$$groupL[$g]}) {
		    ($gBreak=1, last) if !$cloneHash{$p};
		}
		next unless $gBreak;

		my @newProbes;
		while ($$groupH{$cProbes[0]}<=$g) {
		    push @newProbes, shift @cProbes;
		}
		next unless @newProbes;
		createNewClone($clone,$suffix,\@newProbes,$clones,$probes);
		$suffix++;
		$result = $chimerical = 1;

	    }
	    next unless $chimerical;
	    createNewClone($clone, $suffix, \@cProbes,$clones,$probes);
	    delete $$probes{$clone};
	}
    }
    return $result;
}

sub orderSubsets {
    my ($groupL, $component, $clones, $probes) = @_;
    my $groupH = buildGroupIndex($groupL);
    my %tested;
    my (%lRank, %rRank);
    foreach (keys %$clones) { $lRank{$_} = $rRank{$_} = @{$$groupL[$$groupH{$_}]}; }
    foreach my $probe (@$component) {
	foreach my $clone (@{$$clones{$probe}}) {
	    next if $tested{$clone}++;
	    my @cProbes=sort {$$groupH{$a}<=>$$groupH{$b}} @{$$probes{$clone}};
	    print "Clone $clone is @cProbes\n";
	    my $gL = $$groupH{$cProbes[0]};
	    my $gR = $$groupH{$cProbes[@cProbes-1]};
	    next if $gL==$gR;  ## probes from a single groups
	    my @leftEnd = grep($$groupH{$_} == $gL, @cProbes);
	    my $newRank = +@leftEnd;
	    print "Left End of $clone is @leftEnd; setting ranks to $newRank\n";
	    foreach (@leftEnd) { 
		$lRank{$_} = $newRank
		    if !defined($lRank{$_}) || $lRank{$_} > $newRank;
	    }
	    my @rightEnd = grep($$groupH{$_} == $gR, @cProbes);
	    my $newRank = +@rightEnd;
	    print "Right End of $clone is @rightEnd; setting ranks to $newRank\n";
	    foreach (@rightEnd) {
		$rRank{$_} = $newRank
		    if !defined($rRank{$_}) || $rRank{$_} > $newRank;
	    }
	}
    }
    print "Ranks:\n";
    foreach (@$component) { print "$_ $lRank{$_} $rRank{$_}; "; };
    print "\n\n\n";
    die;
    
#    if (($lRank{$a} + $rRank{$b}) < groupSize) can't resolve this group.
#    if (($lRank{$a} < $lRank{$b}) && ($rRank{$a} >= $rRank{$b})) {$a;
#    } elsif (($rRank{$a} < $rRank{$b}) && ($lRank{$a} >= $lRank{$b})) {$b;
	
    my @newGroupL;
    foreach my $group (@$groupL) {
#	push @newGroupL, @{groupByHashValue($group, \%rank)};
    }
    return \@newGroupL;
}

sub readClones {
    my (%clones,%probes);
    my $clone=0;
    while (my $in = <STDIN>) {
	chomp($in);
	next unless $in;
	my @probes = split(//,$in);
	foreach my $p (@probes) {
	    $clones{$p} ||= [];
	    push @{$clones{$p}}, $clone;
	}
	$probes{$clone++}=\@probes;
    }
    return (\%clones, \%probes);
}

sub printHashOfLists {
    my ($h) = @_;
    foreach (keys %$h) {
	print("$_: ", join(',', @{$$h{$_}}), "\n");
    }
}

sub main {
    my ($clones, $probes) = readClones();
 
    print "Probe: clones\n";  printHashOfLists($clones);
    print "Clone: probes\n";  printHashOfLists($probes);

    my @components = findComponents([keys %$probes], undef, $clones, $probes);
    print "There are ", 0+@components, " components:\n";
    foreach (@components) {	print(join(",", @$_), "\n");    }
    
    foreach (@components) { 
	print("\n\nOrdering Component ",join(",", @$_), "\n");
	my $orderedComponent;
	while (1) {
	    $orderedComponent = orderComponent($_, $clones, $probes);
	    print join(' ; ', map { join(',', @$_) } @$orderedComponent), " (1)\n";
	    next if splitChimerae($orderedComponent, $_, $clones, $probes);
	    print join(' ; ', map { join(',', @$_) } @$orderedComponent), " (2)\n";
	    $orderedComponent = orderSubsets($orderedComponent, $_, $clones, $probes);	
	    print join(' ; ', map { join(',', @$_) } @$orderedComponent), " (3)\n";
	    last unless splitChimerae($orderedComponent, $_, $clones, $probes);
	}
	print join(' ; ', map { join(',', @$_) } @$orderedComponent), "\n";
    }
    print "Final Clones: \n";  printHashOfLists($probes);
}

main();
