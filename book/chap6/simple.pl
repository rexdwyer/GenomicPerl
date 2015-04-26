#!/usr/bin/perl -I . -I /home/dwyer/book/perl

use strict;
use SimpleReader;


my $r = SimpleReader->new("file.txt");
while (my $seq = $r->readSeq()) {
    my $id = $r->seqId();
    print "\n\n", $r->fileName(), ", $id:\n";
    print "$seq\n";
}
$r->close();
