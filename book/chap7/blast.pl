#!/usr/bin/perl -I /home/dwyer/book/perl

use strict;
use Util;
use SeqReader;

my %lod;

sub fillLod {
    my ($matrixFile) = @_;
    open LOD, $matrixFile;
    my ($trash,@residues) = split /\s+/, <LOD>;
    while (<LOD>) {
	my ($r,@scores) = split;
	foreach (0..19) {
	    $lod{$r.$residues[$_]} = $scores[$_];
	}
    }
}

my @residues = split //,"acdefghiklmnpqrstvwy";   ## list of 20 amino acids

sub findSimilarWords {
    my ($word) = @_;
    my ($w1,$w2,$w3) = split //,$word;
#   print "Similar to $word:\n";
    return [] if $lod{$w1.$w1}+$lod{$w2.$w2}+$lod{$w3.$w3}<11;
    my @similar;
    foreach my $r1 (@residues) {
	foreach my $r2 (@residues) {
	    my $t = 11-$lod{$w1.$r1}-$lod{$w2.$r2};
	    foreach my $r3 (@residues) {
		push @similar, "$r1$r2$r3" if $lod{$w3.$r3}>=$t;
	    }
	}
    }
#   foreach (@similar) { print " $_"; }
#   print "\n";
    return \@similar;
}
	
my %similarWords;

sub preprocessQuery {
    my ($query) = @_;
    my %similarPositions;
    for (my $i=0; $i<(length $query)-2; $i++) {
	my $word = substr($query,$i,3);
	$similarWords{$word} ||= findSimilarWords($word);
	foreach (@{$similarWords{$word}}) {
	    $similarPositions{$_} ||= [ ];
	    push @{$similarPositions{$_}}, $i+1;
	}
    }
#   print scalar keys %similarPositions, " keys:\n";
#   foreach (sort keys %similarPositions) {
#       print "$_: ", join(",", @{$similarPositions{$_}}), "\n";
#   }
    return \%similarPositions;
}


sub extendHit {
    my ($query,$qPos,$qId,$target,$tPos,$tId) = @_;
    my @target = ("-", split //,$target);
    my @query = ("-", split //,$query);
    my ($lo,$hi) = (0,2);

    my $maxscore = my $score =
	$lod{$target[$tPos].$query[$qPos]}
        + $lod{$target[$tPos+1].$query[$qPos+1]}
        + $lod{$target[$tPos+2].$query[$qPos+2]};

    ### Try to grow gap-free alignment to right.
    my $ilim = min(@target - $tPos, @query - $qPos);
    for (my $i = 3; $i<$ilim; $i++) {
	$score += $lod{$target[$tPos+$i].$query[$qPos+$i]};
	last if $score < 0;
	($maxscore,$hi) = ($score,$i) if $score > $maxscore;
    }

    ### Try to grow gap-free alignment to left.
    my $score = $maxscore;
    my $ilim = min($tPos, $qPos);
    for (my $i = 1; $i<$ilim; $i++) {
	$score += $lod{$target[$tPos-$i].$query[$qPos-$i]};
	last if $score < 0;
	($maxscore,$lo) = ($score,$i) if $score > $maxscore;
    }

    ### Return best alignment as triple.
    my $len = $hi+$lo+1;
    return 
	[$maxscore, 
	 substr($query, $qPos-$lo-1, $len)." at ".($qPos-$lo-1)." in $qId",
	 substr($target, $tPos-$lo-1, $len)." at ".($tPos-$lo-1)." in $tId"];
}	
	 

sub scanTarget {
    my ($target,$tId,$query,$qId,$queryIndex) = @_;
    my @alignments;
    for (my $i=0; $i<(length $target)-2; $i++) {
	my $word = substr($target,$i,3);
	my $tPos = $i+1;
	foreach my $qPos (@{$$queryIndex{$word}}) {
#	    print "Hit $word at $tPos in target and $qPos in query.\n";
	    my $alignment = extendHit($query,$qPos,$qId,$target,$tPos,$tId);
#	    print(join("\n", @$alignment), "\n\n");
	    push @alignments, $alignment;
	}
    }
    @alignments = sort { $$b[0] <=> $$a[0] } @alignments;
    $#alignments = min( $#alignments, 24);

#   foreach (@alignments) { print "$$_[0] "; }
#   print "\n";
    return @alignments;
}

### MAIN PROGRAM
fillLod($ARGV[2] || "blosum62");
my $queryFile = new SeqReader $ARGV[0];
while (my $query=$queryFile->readSeq()) {
    my $qId = $queryFile->seqId();
    my $queryIndex = preprocessQuery($query);
    my $targetFile = new SeqReader $ARGV[1];
    my @bestAlignments;
    while (my $target = $targetFile->readSeq()) {
	my $tId = $targetFile->seqId();
	print "\nQuery: $qId\nTarget: $tId\n";
	
	my @newAlignments = 
	    scanTarget($target,$tId,$query,$qId,$queryIndex);
	### Add new alignments to sorted list.
	push @bestAlignments, @newAlignments;
	@bestAlignments = 
	    sort { $$b[0] <=> $$a[0] } @bestAlignments;
	### Cut off list at 25 alignments.
	$#bestAlignments = min($#bestAlignments,24);
    }
    $targetFile->close();
    print "********** 25 Best Alignments for $qId:\n\n";
    foreach (@bestAlignments) {
	print(join("\n", @$_), "\n\n");
    }
}
$queryFile->close();

