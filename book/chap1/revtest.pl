#!/usr/bin/perl
use strict;   ## request conservative error checking

my $x = "foobar";
my $x = reverse $x, $x;
print "$x\n";

my $x = "foobar";
my $x = reverse($x,$x);
print "$x\n";

my $x = "foobar";
my $x = scalar reverse $x,$x;
print "$x\n";

my $x = "foobar";
my $x = scalar reverse($x,$x);
print "$x\n";
