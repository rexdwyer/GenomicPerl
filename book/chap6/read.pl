#!/usr/bin/perl -I . -I /home/dwyer/book/perl

use strict;
use SeqReader;

my @readers;
#push @readers, new SeqReader("file.gb");
push @readers, new SeqReader("turkey.fasta");
#push @readers, new SeqReader("file.txt");

foreach my $r (@readers) {
    print "\n\n", $r->fileName(), ":  ", $r->fileType(), "\n";
    my $seq = $r->readSeq();
    my $id = $r->seqId();
    print "$id:\n$seq\n\n";
    my $seq = $r->readSeq();
    my $id = $r->seqId();
    print "$id:\n$seq\n\n";
    $r->close();
}

push @readers, new GenBankReader("file.fasta");



