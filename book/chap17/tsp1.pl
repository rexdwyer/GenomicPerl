#!/usr/bin/perl -I . -I /home/dwyer/book/perl

use strict;

my (@weight, %probeIndex, %clones, $numProbes, @probes, %score,%from);

## Read in strings, build suffix trees, create weighted graph.
$|=1;
my $clone;
while (my $in=<STDIN>) {
    $clone = $1 if $in =~ s/^(.*):\s*//;
    chomp($in);
    next unless $in =~ /\S/;   # skip whitespace lines;
    $in =~ s/^\s+//;
    foreach my $probe (split /\s+/, $in) {
	$probeIndex{$probe} ||= ++$numProbes;
	$clones{$probe} ||= [];
	push @{$clones{$probe}}, $clone;
    }
}

print ((keys %clones), "\n");
foreach my $p1 (keys %clones) {
    foreach my $p2 (keys %clones) {
	last if $p1 eq $p2;
	my %tmp = map(($_=>2), @{$clones{$p1}});
	my $dist = @{$clones{$p1}} + @{$clones{$p2}};
	foreach my $cl (@{$clones{$p2}}) { $dist -= $tmp{$cl}; }
	$weight[$probeIndex{$p1}][$probeIndex{$p2}] = $dist;
	$weight[$probeIndex{$p2}][$probeIndex{$p1}] = $dist;
    }
}

my @q = map("$_,0", (1..$numProbes));
push @q, "clean";

while (@q) {
    my $item = shift @q;
    if ($item eq "clean") {
	my %useful;
	foreach (@q) { $useful{$_} = 1; }
	print "cleaning all but ", scalar keys %useful, " items.",
	      "  Next item is $q[0].\n";
	foreach (keys %score) {
	    delete $score{$_} unless $useful{$_};
	}
	push @q, $item if @q;
	next;
    }
#    print "$item\n";
    my ($second,$setInt) = split(',', $item);
#    print "@set = $setInt\n";
#    printf "%o,", $setInt;  # octal
    my $newSetInt = $setInt + 1<<$second;
    my $saveIt=0;
    foreach my $first (1..$numProbes) {
	next if $newSetInt & (1<<$first);
	my $newScore = $score{$item} + $weight[$first][$second];
	my $newSet = "$first,$newSetInt";
	if (!defined($score{$newSet})) {
	    $score{$newSet} = $newScore;
	    push @q, $newSet;
	    $saveIt=1;
	} elsif ($newScore < $score{$newSet}) {
	    $score{$newSet} = $newScore;
	    $saveIt=1;
	}
    }
    if (!$saveIt) {
#	print "\n";
	delete $score{$item};
    }
}

my $all = (1<<($numProbes+1))-1;
foreach my $first (1..$numProbes) {
    my $allButFirst = $all - (1<<$first);
    my $item = "$first,$allButFirst";
    print "$item:   $score{$item}\n";
}


