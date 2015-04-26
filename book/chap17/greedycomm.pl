#!/usr/bin/perl -I . -I /home/dwyer/book/appC -I /home/dwyer/book/perl

use strict;
use UnionFind;

my (@affinities, %clones);
my $clones = readClonesIntoHash();
my $affinities = computeAffinities($clones);
my $path = buildPaths($affinities);
showPaths($path);
exit(0);


sub buildPaths {
    my ($affinities) = @_;
    $affinities = [sort {$$b[0] <=> $$a[0]} @$affinities];
    my (%degree, %path);
    my $uf = new UnionFind;
    foreach (@$affinities) {
	my ($w, $p1, $p2) = @$_;
#	print "[$p1, $p2]; $w\n";
	next if $uf->inSameSet($p1, $p2);
	next if $degree{$p1}==2;
	next if $degree{$p2}==2;
	$uf->union($p1, $p2);
	$path{$p1} ||= [ ];
	push @{$path{$p1}}, $p2;
	$degree{$p1}++;
	$path{$p2} ||= [ ];
	push @{$path{$p2}}, $p1;
	$degree{$p2}++;
    }
    return \%path;
}

sub showPaths {
    my ($path) = @_;
    my %visited;
    my $componentCount;
    foreach my $p (keys %$path) {
	next if $visited{$p};
	next if @$path{$p}==2;  ## not an endpoint
	$componentCount++;
	print "$p";
	$visited{$p}++;
	while (my @pp = grep(!$visited{$_}, @{$$path{$p}})) {
##	    print affinity of $p and $pp[0] somehow... 
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

sub computeAffinities {
    my ($clones) = @_;
    ## For each pair of probes (i,j), compute 
    ##    | Ci intersect Cj | / | Ci union Cj |  
    ## where Ci = clones of probe i
    my @affinities;
    foreach my $p1 (keys %$clones) {
	foreach my $p2 (keys %$clones) {
	    last if $p1 eq $p2;
	    my %tmp = map(($_=>1), @{$$clones{$p1}});
	    my $commonClones;
	    foreach my $cl (@{$$clones{$p2}}) { $commonClones += $tmp{$cl}; }
	    my $fractionShared =
	      $commonClones/(@{$$clones{$p1}}+@{$$clones{$p2}}-$commonClones);
	    push @affinities, [$commonClones, $p1, $p2] if $fractionShared>0;
	    if ($fractionShared==1) {
		print "$p1 and $p2 are indistinguishable; removing $p1\n";
		delete $$clones{$p1};
		while (@affinities && $affinities[-1][1] eq $p1) {
		    pop @affinities;
		}
		last;
	    }
	}
    }
    return \@affinities;
}
