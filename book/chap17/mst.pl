#!/usr/bin/perl -I . -I /home/dwyer/book/appC -I /home/dwyer/book/perl

use strict;
use UnionFind;

my (%weight, %clones, %probes, $numProbes, %score, %from);
my $fractionChimeric = 1.0/3.0;

## Read in strings, build suffix trees, create weighted graph.
$|=1;
my $clone;
while (my $in=<STDIN>) {
    $clone = $1 if $in =~ s/^(.*):\s*//;
    chomp($in);
    next unless $in =~ /\S/;   # skip whitespace lines;
    $in =~ s/^\s+//;
    foreach my $probe (split /\s+/, $in) {
	$clones{$probe} ||= [];
	push @{$clones{$probe}}, $clone;
	push @{$probes{$clone}}, $probe;
    }
}

# For each pair of probes (i,j), compute 
#    | Ci intersect Cj | / | Ci union Cj |  
# where Ci = clones of probe i
foreach my $p1 (keys %clones) {
    print "p";
    foreach my $p2 (keys %clones) {
	last if $p1 eq $p2;
	my %tmp = map(($_=>1), @{$clones{$p1}});
	my $commonClones;
	foreach my $cl (@{$clones{$p2}}) { $commonClones += $tmp{$cl}; }
	my $fractionShared = 
	    $commonClones / (@{$clones{$p1}}+@{$clones{$p2}}-$commonClones);
	if ($fractionShared==1) {
	    print "$p1 and $p2 are indistinguishable; removing $p1\n";
	    delete $clones{$p1};
	    last;
	}
	$weight{$p1,$p2} = $weight{$p2,$p1} = $fractionShared 
	    if $fractionShared>0;
    }
}

my @thresholds = map([$_, weakestLink($_)], (keys %probes));
   # find weakest link in each clone = threshold of chimerism.
@thresholds = sort { $$a[1]<=>$$b[1] } @thresholds;
#print map("$$_[1],",@thresholds),"\n";
$#thresholds = int($fractionChimeric*@thresholds)-1;
my $tooWeak = $thresholds[-1][1];
foreach (@thresholds) { splitChimera($$_[0], $tooWeak); }

foreach (sort keys %probes) {
    print "$_:\n";
    my $s = join(' ', @{$probes{$_}}) . ' ';
    while ($s =~ s/(.{1,70}) //) {  print "  $1\n"; }
    print "$s\n" if $s;
}


sub weakestLink {
    print "W";
    my ($clone) = @_;
    my @probes = grep($clones{$_}, @{$probes{$clone}});
    $probes{$clone} = \@probes;
    my @edges;
    foreach my $p1 (@probes) {
	foreach my $p2 (@probes) {
	    last if $p1 eq $p2;
	    push @edges, [$weight{$p1,$p2}, $p1, $p2];
	}
    }
    @edges = sort { $$b[0]<=>$$a[0] } @edges;
    my $uf = new UnionFind;
    my $weakest;
    foreach (@edges) {
	my ($w, $p1, $p2) = @$_;
	if ($uf->find($p1) ne $uf->find($p2)) {
	    $weakest = $w;
	    $uf->union($p1,$p2);
	}
    }
    return $weakest;
}
    
sub splitChimera {
    my ($clone, $limit) = @_;
    print "S";
    my @edges;
    foreach my $p1 (@{$probes{$clone}}) {
	foreach my $p2 (@{$probes{$clone}}) {
	    last if $p1 eq $p2;
	    push @edges, [$weight{$p1,$p2}, $p1, $p2];
	}
    }
    @edges = sort { $$b[0]<=>$$a[0] } @edges;
    my $uf = new UnionFind;
    foreach (@edges) {
	my ($w, $p1, $p2) = @$_;
	last if $w<=$limit;
	$uf->union($p1,$p2);
    }
    my %set;
    foreach my $p (@{$probes{$clone}}) {
	my $fp = $uf->find($p);
	$set{$fp} ||= [];
	push @{$set{$fp}}, $p;
    }
    return if scalar(keys %set)==1;
    my $i=1;
    foreach my $fp (keys %set) {
	my $nuClone = "$clone~$i";
	$i++;
	$probes{$nuClone} = $set{$fp};
	foreach my $p (@{$probes{$nuClone}}) {
	    my $ref = $clones{$p};
	    foreach my $j (0..$#$ref) {
		($$ref[$j] = $nuClone, last) if $$ref[$j] eq $clone;
	    }
	}
    }
    delete $probes{$clone};
}


