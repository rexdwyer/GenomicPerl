#!/usr/bin/perl -I /home/dwyer/book/perl

use strict;
use Util;

$|=1;
while (my $score = <STDIN>) {
    my $id1 = <STDIN>;
    my $id2 = <STDIN>;
    my $seq1 = <STDIN>;
    my $seq2 = <STDIN>;
    my $cnt=0;
    foreach my $i (0..((length $seq1)-1)) {
	$cnt++ if (substr($seq1,$i,1) eq substr($seq2,$i,1));
    }
    print($score,$id1,$id2,$seq1,$seq2) if (2*$cnt >= length $seq1);
}
