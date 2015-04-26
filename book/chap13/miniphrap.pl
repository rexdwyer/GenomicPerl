#!/usr/bin/perl

use strict;   ## request conservative error checking
### These define classes that may be repeatedly instantiated.
use DnaRead;
use Alignment;
use Contig;

### These simply define groups of related subroutines (packages).
use JudiciousAligner;
use QualityReviser;
use ContigBuilder;
use ConsensusBuilder;

### See Chapters 4 and 6.
use Util;
use SeqReader;

main();

sub readDnaFile {
    my ($filename) = @_;
    my @readList;
    my $reader = SeqReader->new($filename);
    warn "No quality information in this file; supplying defaults"
	if ($reader->fileType() != "Fastq");
    while (my $seq = $reader->readSeq()) {
	my $id = $reader->seqId();
	my $qual = $reader->seqQuality();
	push @readList, DnaRead->new($id, $seq, $qual);
    }
    return \@readList;
}

sub main {
    my $readList = readDnaFile($ARGV[0]);
    my $alignmentsByPair = makeAlignmentMaps($readList);
    reviseQualities($readList, $alignmentsByPair, 0);
    reviseQualities($readList, $alignmentsByPair, 1);
    my $contigList = assignReadsToContigs($readList, $alignmentsByPair);
    reviseQualities($readList, $alignmentsByPair, 2);
    constructConsensusSequences($contigList);
    $contigList = [sort {$b->Score <=> $a->Score} @$contigList];
    foreach (@$contigList) { $_->dump(); }
    warn "\a\a\a\a\a\n";
}
