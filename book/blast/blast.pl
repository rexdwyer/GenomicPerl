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

my @residues = split //,"acdefghiklmnpqrstvwy";   # list of 20 amino acids

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
	
sub preprocessQuery {
    my ($query) = @_;
    my (%similarWords, %similarPositions);
    for (my $i=0; $i<(length $query)-2; $i++) {
	my $word = substr($query,$i,3);
	$similarWords{$word} ||= findSimilarWords($word);
	foreach (@{$similarWords{$word}}) {
	    push @{$similarPositions{$_}}, $i+1;
	}
    }
#   print scalar keys %similarPositions, " keys:\n";
#   foreach (sort keys %similarPositions) {
#       print "$_: ", join(",", @{$similarPositions{$_}}), "\n";
#   }
    return \%similarPositions;
}


sub expandHit {
    my ($query,$qpos,$qid,$target,$tpos,$tid) = @_;
    my @target = ("-", split //,$target);
    my @query = ("-", split //,$query);
    my ($lo,$hi) = (0,2);

    my $maxscore = my $score =
	$lod{$target[$tpos].$query[$qpos]}
        + $lod{$target[$tpos+1].$query[$qpos+1]}
        + $lod{$target[$tpos+2].$query[$qpos+2]};

    ### Try to grow gap-free alignment to right.
    my $ilim = min(@target - $tpos, @query - $qpos);
    for (my $i = 3; $i<$ilim; $i++) {
	$score += $lod{$target[$tpos+$i].$query[$qpos+$i]};
	last if $score < 0;
	($maxscore,$hi) = ($score,$i) if $score > $maxscore;
    }

    ### Try to grow gap-free alignment to left.
    my $score = $maxscore;
    my $ilim = min($tpos, $qpos);
    for (my $i = 1; $i<$ilim; $i++) {
	$score += $lod{$target[$tpos-$i].$query[$qpos-$i]};
	last if $score < 0;
	($maxscore,$lo) = ($score,$i) if $score > $maxscore;
    }

    ### Return best alignment as triple.
    my $len = $hi+$lo+1;
    return [$maxscore, 
	    substr($query, $qpos-$lo-1, $len)." at ".($qpos-$lo-1)." in $qid",
	    substr($target, $tpos-$lo-1, $len)." at ".($tpos-$lo-1)." in $tid"];
}	
	 

sub scanTarget {
    my ($target,$tid,$query,$qid,$queryMap) = @_;
    my @alignments;
    for (my $i=0; $i<(length $target)-2; $i++) {
	my $word = substr($target,$i,3);
	my $tpos = $i+1;
	foreach my $qpos (@{$$queryMap{$word}}) {
#	    print "Hit $word at $tpos in target and $qpos in query.\n";
	    my $alignment = expandHit($query,$qpos,$qid,$target,$tpos,$tid);
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

fillLod($ARGV[3] || "blosum62");
my $queryFile = new SeqReader $ARGV[0];
while (my @query=$queryFile->readSeq()) {
    my ($query,$qid,$qnotes) = @query;
    my $queryMap = preprocessQuery($query);
    my $targetFile = new SeqReader $ARGV[1];
    my @bestAlignments;
    while (my @target = $targetFile->readSeq()) {
	my ($target,$tid,$tnotes) = @target;
	print "\nQuery: $qid\nTarget: $tid\n";
	
	my @newAlignments = scanTarget($target,$tid,$query,$qid,$queryMap);
	### Add new alignments to sorted list.
	push @bestAlignments, @newAlignments;
	@bestAlignments = sort { $$b[0] <=> $$a[0] } @bestAlignments;
	### Cut off list at 25 alignments.
	$#bestAlignments = min($#bestAlignments,24);
    }
    $targetFile->close();
    print "****************** 25 Best Alignments for $qid: *****************\n\n";
    foreach (@bestAlignments) {
	print(join("\n", @$_), "\n\n");
    }
}
$queryFile->close();

