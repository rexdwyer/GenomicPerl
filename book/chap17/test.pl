#!/usr/bin/perl -I . -I /home/dwyer/book/appC -I /home/dwyer/book/perl

use strict;
use UnionFind;

my $clones = readClonesIntoHash();
my $amities = computeAmities($clones);
showAmities($amities);
exit(0);
my $path = buildPaths($amities);
showPaths($path);

sub showAmities {
    my ($amity) = @_;
    my (%a,%b);
    foreach (@$amity) {
	my ($w, $p1, $p2) = @$_;
	$b{$p1} = $b{$p2} = 1;
	$a{$p1,$p2} = $a{$p2,$p1} = $w;
    }
    foreach my $p1 (sort keys %b) {
	print "$p1";
	foreach my $p2 (sort keys %b) {
	    print "&-- " and next if $p1 eq $p2;
	    print "&   " and next if !defined($a{$p1,$p2});
	    printf "&.%02d", int(0.5+100*$a{$p1,$p2});
	}
	print "\\\\\n";
    }
}



sub buildPaths {
    my ($amities) = @_;
    $amities = [sort {$$b[0] <=> $$a[0]} @$amities];
    my (%degree, %path);
    my $uf = new UnionFind;
    foreach (@$amities) {
	my ($w, $p1, $p2) = @$_;
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
	next if @{$$path{$p}}==2;  ## not an endpoint
	$componentCount++;
	print "$p";
	$visited{$p}++;
	while (my @pp = grep(!$visited{$_}, @{$$path{$p}})) {
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
	foreach my $p2 (keys %$clones) {
	    last if $p1 eq $p2;
	    my %tmp = map(($_=>1), @{$$clones{$p1}});
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

