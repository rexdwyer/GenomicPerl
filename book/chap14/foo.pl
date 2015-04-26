#!/usr/bin/perl -I . -I /home/dwyer/book/perl

use strict;
use Util;
use CSPredictor;

my $csp = new CSPredictor "foo.freqs";

my $s;
my @b = qw(a c g t);
foreach (0..999) { $s .= $b[rand(4)]; }
print "$s\n";

my $predictions = $csp->predictCodingSequences($s);

foreach (@$predictions) { print join(",", @$_), "\n"; }

exit(0);
