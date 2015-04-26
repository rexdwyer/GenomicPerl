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
    print "Probe Biases:\n";
    foreach (keys %bias) { print "  $_ $bias{$_}\n"; };
    print "\n\n\n";
    return %bias;
}

sub orderComponent {
    my ($component, $clones, $probes) = @_;

    my ($bestSplitter,$splitters) = findSplitters($component,$clones,$probes);
    my %bias = computeProbeBiases($component,$clones,$probes,$bestSplitter,$splitters);

    my @order = sort { $bias{$a} <=> $bias{$b} } @$component;

    my @groups = ( [ $order[0] ] );
    my %group;
    $group{$order[0]} = @groups;
    my $lastbias = $bias{$order[0]};

    for (my $i=1; $i<@order; $i++) {
	push @groups, [] if $bias{$order[$i]} > $lastbias;
	push @{$groups[@groups-1]}, $order[$i];
	$lastbias = $bias{$order[$i]};
	$group{$order[$i]} = @groups-1;
    }

    print join(' ; ', map { join(',', @$_) } @groups), "\n";
    return (\%group, \@groups);
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

sub splitChimerae {
    my ($groupH, $groups, $component, $clones, $probes) = @_;
    my $result=0;
    my %tested;
    foreach my $probe (@$component) {
	foreach my $clone (@{$$clones{$probe}}) {
	    next if $tested{$clone}++;
	    my @cProbes=sort {$$groupH{$a}<=>$$groupH{$b}} @{$$probes{$clone}};
	    my $glo = $$groupH{$cProbes[0]}+1;
	    my $ghi = $$groupH{$cProbes[@cProbes-1]}-1;
	    next if $glo>$ghi;  ## clones from 1 or 2 adjacent groups
	    my %cloneHash;
	    foreach my $p (@{$$probes{$clone}}) { $cloneHash{$p}=1; }

	    my $suffix=1;
	    my $chimerical=0;
	    foreach my $g ($glo..$ghi) {
		my $gBreak = 0;
		foreach my $p (@{$$groups[$g]}) {
		    ($gBreak=1, last) unless $cloneHash{$p};
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
 
    print "Probes: clones\n";  printHashOfLists($clones);
    print "Clones: probes\n";  printHashOfLists($probes);

    my @components = findComponents([keys %$probes], undef, $clones, $probes);
    print "There are ", 0+@components, " components:\n";
    foreach (@components) {	print(join(",", @$_), "\n");    }
    
    foreach (@components) { 
	print("\n\nOrdering Component ",join(",", @$_), "\n");
	while (splitChimerae(orderComponent($_, $clones, $probes),
			     $_, $clones, $probes)) {
	    next;
	}
    }
    print "Final Clones: \n";  printHashOfLists($probes);
}

main();
