package JudiciousAligner;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(makeAlignmentMaps);

use strict;
use Alignment;
use Util;

my $SUBSTPENALTY = -2;
my $INDELPENALTY = -3; ## -8;

sub makeAlignmentMaps {
    my ($readList) = @_;
    my $fourteenmerHash = build14merHash($readList);
    my $candidateBands = findCandidateBands($fourteenmerHash);    
    my $alignmentsByPair = alignBands($readList, $candidateBands);
    ## returns a hash mapping pairs of reads to their best alignments.
    ## Index the alignments by read instead of just by pair.
    my @alignmentsByRead;
    foreach my $pair (keys %$alignmentsByPair) {
	my ($r1, $r2) = split($;, $pair);
	my $alignments =$$alignmentsByPair{$pair};
	$alignmentsByRead[$r1] ||= [ ];
	push @{$alignmentsByRead[$r1]}, @$alignments;
	$alignmentsByRead[$r2] ||= [ ];
	push @{$alignmentsByRead[$r2]}, @$alignments;
    }
    foreach my $read (@$readList) {
	$read->setAlignments
	    ($alignmentsByRead[$read->SerialNumber()] || [ ]);
    }
    return $alignmentsByPair;
}

sub build14merHash {
    my ($readList) = @_;
    my %hash14;
    foreach my $read (@$readList) {
	my $seq = $read->Bases();
	my $readIndex = $read->SerialNumber();
	my $lasti = (length $seq) - 15;
	for (my $i = 0; $i < $lasti; $i++) {
	    my $word = substr($seq, $i, 14);
	    $hash14{$word} ||= [ ];
	    push @{$hash14{$word}}, [$readIndex, $i];
	}
    }
    return \%hash14;
}

sub findCandidateBands {
    my ($hash14mers) = @_;
    ## Build hash mapping pairs of reads to offsets with matching 14-mers.
    my %pairHash;
    foreach (values %$hash14mers) { ## list of all locations with some 14-mer.
	next if $#$_ < 1;  ## not strictly necessary.
	my @wordLocs = sort {$$a[0]<=>$$b[0] || $$a[1]<=>$$b[1]} @$_;

	for (my $i = 0; $i < @wordLocs; $i++) {
	    my ($read1, $loc1) = @{$wordLocs[$i]};
	    for (my $j = $i+1; $j < @wordLocs; $j++) {
		my ($read2, $loc2) = @{$wordLocs[$j]};
		$pairHash{$read1,$read2} ||= [ ];
		push @{$pairHash{$read1,$read2}}, ($loc2-$loc1);
	    }
	}
    }
    
    ## From above, extract list of bands to search for each pair.
    foreach my $pair (keys %pairHash) {
	my @bands;
	my @endPoints;
	foreach my $offset (@{$pairHash{$pair}}) {
	    ## +1 is left endpoint; -1 is right endpoint.
	    push @endPoints, [$offset-14,1], [$offset+28,-1];
	}
	@endPoints =
	    sort { $$a[0] <=> $$b[0] || $$b[1] <=> $$a[1] } @endPoints;
	
	my $nesting = 0;
	my $currLeft;
	foreach (@endPoints) {
	    my ($endPoint, $deltaNesting) = @$_;
	    $currLeft = $endPoint
		if ($nesting == 0 && $deltaNesting == 1);
	    push @bands, [$currLeft, $endPoint]
		if ($nesting == 1 && $deltaNesting == -1);
	    $nesting += $deltaNesting;
	}
	$pairHash{$pair} = \@bands;
    }
    return \%pairHash;
}

sub alignBands {
    my ($readList, $candidateBands) = @_;
    my $MINSCORE = 30; 
    my %alignmentsByPair;
    ## Try to align each band for each pair.
    foreach my $pair (keys %$candidateBands) {
	my ($r1, $r2) = @$readList[split($;, $pair)];
	my ($sn1, $sn2) = ($r1->SerialNumber(), $r2->SerialNumber());
	my @alignments; 
	foreach my $band (@{$$candidateBands{$pair}}) {
	    my $alignment = alignBand($r1,$r2, $band);
#	    $alignment->dump();
	    push @alignments, $alignment 
		if $alignment->RawScore() > $MINSCORE
	}
	if (@alignments) {  ## Did any alignments meet MINSCORE criterion?
	    @alignments =   ## Mark best & save all.
		sort {$b->RawScore() <=> $a->RawScore()} @alignments;
	    $alignmentsByPair{$sn1,$sn2} = \@alignments;
	}
    }
    return \%alignmentsByPair;
}

sub alignBand {
    my ($read1, $read2, $band) = @_;
#   print "Enter alignBand @_ @{$_[2]}\n";
    my ($lo, $hi) = @$band;
    my ($s, $t) = ($read1->Bases(), $read2->Bases());

    my %M; 

    ## Add offset to position in $s to get position in $t,
    ## OR substract offset from position in $t to get position in $s.
    ## Since position in $t can't be negative, 
    my $firstI = max(0, -$hi);
    ## Similarly for other end:
    my $lastI = min(length($s), length($t)-$lo)-1;

    ## Fill in part of scoring matrix associated with this band.
#   print "$firstI <= i <= $lastI\n";
    foreach my $i ($firstI ..$lastI) {
	## Reasoning similar to above:
	my $firstJ = max(0, $i + $lo);
	my $lastJ = min(length($t)-1, $i + $hi);
#	print "i=$i: $firstJ <= j <= $lastJ\n";
	foreach my $j ($firstJ ..$lastJ) {
	    $M{$i,$j} = 
                max(0,
		    $M{$i-1,$j} + $INDELPENALTY,
                    $M{$i,$j-1} + $INDELPENALTY,
                    $M{$i-1,$j-1} 
		    +((substr($s,$i,1) eq substr($t,$j,1))||$SUBSTPENALTY));
	}
    }

    ## Find best local alignment score anywhere in this band.
    my ($bestKey, $bestScore);
    foreach my $key (keys %M) {
	($bestKey, $bestScore) = ($key, $M{$key}) if $M{$key} > $bestScore;
    }

    ## Now trace back through matrix to produce a string representing the
    ## best local alignment.
    my ($i,$j) = split $;, $bestKey;
    my @path;
    my ($finalI, $finalJ);
    while ($M{$i,$j} > 0) {
	($finalI, $finalJ) = ($i, $j);
	my ($sBase, $tBase) = (substr($s,$i,1), substr($t,$j,1));
	if ($sBase eq $tBase) {
	    push @path, "M";  ## dial "M" for "MATCH"
	    $i--; $j--;
	} elsif ($M{$i,$j} == $M{$i-1,$j-1} + $SUBSTPENALTY) {
	    push @path, "S";  ## substitution in $t w.r.t. $s
	    $i--; $j--;
	} elsif ($M{$i,$j} == $M{$i-1,$j} + $INDELPENALTY) {
	    push @path, "D";  ## deletion in $t w.r.t. $s
	    $i--;
	} else { ## ($M{$i,$j}==$M{$i,$j-1}+$INDELPENALTY), by elimination.
	    push @path, "I";  ## insertion in $t w.r.t. $s
	    $j--;
	} 
    }
    ## Return an alignment object.
    return Alignment->new($read1, $read2, $bestScore, $finalI, $finalJ, 
			  join("", reverse @path));
}

