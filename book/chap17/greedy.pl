#!/usr/bin/perl -I . -I /home/dwyer/book/appC -I /home/dwyer/book/perl

use strict;
use UnionFind;

my $clones = readClonesIntoHash();
my $amities = computeAmities($clones);
my $neighbors = assignNeighbors($amities);
showNeighbors($neighbors);
exit(0);

sub assignNeighbors {
    my ($amities) = @_;
    $amities = [sort {$$b[0] <=> $$a[0]} @$amities];
    my %neighbors;
    my $uf = new UnionFind;
    foreach (@$amities) {
	my ($w, $p1, $p2) = @$_;
#	print "[$p1, $p2]; $w\n";
	next if $uf->inSameSet($p1, $p2);
	$neighbors{$p1} ||= [ ];
	$neighbors{$p2} ||= [ ];
	next if @$neighbors{$p1}==2;
	next if @$neighbors{$p2}==2;
	print "[$p1, $p2]; $w\n";
	$uf->union($p1, $p2);
	push @{$neighbors{$p1}}, $p2;
	push @{$neighbors{$p2}}, $p1;
    }
    return \%neighbors;
}

sub showNeighbors {
    my ($neighbors) = @_;
    my %visited;
    my $componentCount;
    foreach my $p (sort keys %$neighbors) {
	next if $visited{$p};
	next if @{$$neighbors{$p}}==2;  ## not an endpoint
	$componentCount++;
	print "$p";
	$visited{$p}++;
	while (my @pp = grep(!$visited{$_}, @{$$neighbors{$p}})) {
##	    print amity of $p and $pp[0] somehow... 
	    $p = $pp[0];
	    print " $p";
	    $visited{$p}++;
	}
	print "\n\n";
    }
    print "\# There are $componentCount sets of ordered probes.\n";
}

sub readClonesIntoHash { ## Read in clones, record probe/clone connections.
    my ($clone, %clones);
    while (my $in=<STDIN>) {
	if ($in =~ /^\#/) { # Echo comment lines: # in column 1.
	    warn $in; 
	    next;
	}
	$clone = $1 if $in =~ s/^(.*):\s*//;
	chomp($in);
	next unless $in =~ /\S/;   # skip whitespace lines;
	$in =~ s/^\s+//;
	foreach my $probe (split /\s+/, $in) {
	    $clones{$probe} ||= [];
	    push @{$clones{$probe}}, $clone;
	}
    }
    return \%clones;
}

sub computeAmities {
    my ($clones) = @_;
    ## For each pair of probes (i,j), compute 
    ##    | Ci intersect Cj | / | Ci union Cj |  
    ## where Ci = clones of probe i
    my @amities;
    foreach my $p1 (keys %$clones) {
	my %tmp = map(($_=>1), @{$$clones{$p1}});
	foreach my $p2 (keys %$clones) {
	    last if $p1 eq $p2;
	    my $commonClones;
	    foreach my $cl (@{$$clones{$p2}}) { $commonClones += $tmp{$cl}; }
	    my $amity = 
	      $commonClones/(@{$$clones{$p1}}+@{$$clones{$p2}}-$commonClones);
	    push @amities, [$amity, $p1, $p2] if $amity>0;
	    if ($amity==1) {
		print "$p1 and $p2 are indistinguishable; removing $p1\n";
		delete $$clones{$p1};
		while (@amities && $amities[-1][1] eq $p1) {
		    pop @amities;
		}
		last;
	    }
	}
    }
    return \@amities;
}

