#!/usr/bin/perl -I . -I /home/dwyer/book/perl

use strict;
use Util;

main();

sub main {
    my $ranSites = generateRandomSites(20,100000);
#   my $ranSites = [0,2,4,7,10,12,16,18,20];
    print "sites: ",  join(',', @$ranSites), "\n";
    printPairwiseDistances($ranSites);
}

sub generateRandomSites{
    my ($count, $range) = @_;
    my @sites;
    push @sites,0;
    foreach (2..$count-1) { push @sites, int(rand($range)); }
    push @sites,$range;
    my @sites = sort { $a <=> $b } @sites;
    my $j=0;
    my $last = $range+1;
    foreach (my $i=1; $i<@sites; $i++) {  ## remove duplicates
	$sites[++$j] = $sites[$i] unless $sites[$i]==$sites[$j];
    }
    $#sites = $j;
    return \@sites;
}

sub printPairwiseDistances {
    my ($sites) = @_;
    foreach (my $sep=1; $sep<@$sites; $sep++) {
	foreach (my $i=$#$sites-$sep; $i>=0; $i--) {
	    print $$sites[$i+$sep]-$$sites[$i], " ";
	}
	print "\n";
    }
}

