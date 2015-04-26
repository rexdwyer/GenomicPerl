package PSSM;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(new maxScore setBackground setProbabilities dumpDefinition
	     findWatsonSignals findCrickSignals);

use strict;

sub max {
    my ($m,@l) = @_;
    foreach (@l) { $m=$_ if $_>$m };
    $m;
}

sub sum {
    my ($s,@l) = @_;
    foreach (@l) { $s+=$_ };
    $s;
}

sub lg { log(shift)/log(2.0); }

sub percent { int(0.5+100*shift); }

sub new {
    my ($this,$description,$signal,$bgndref,$probref) = @_;
    printf STDERR "Creating PSSM: $description\n";
    $this =
	bless {bgndarr=>[],
	       bgndhsh=>{},
	       probarr=>[],
	       probhsh=>{},
	       lodarr=>[],
	       size=>0,
	       signal=>$signal,
	       description=>$description},
	(ref($this) || $this);
    $this->setBackground($bgndref);
    $this->setProbabilities($probref);
    return $this;
}


my @basename = qw(a c g t);
my %basenumber = (a=>0, c=>1, g=>2, t=>3);

sub setBackground {
    my ($this,$hash) = @_;
    foreach (keys(%$hash)) {
	$%{$this->{bgndhsh}}{$_} =
	    $@{$this->{bgndarr}}[$basenumber{$_}] =
		0.01 * $$hash{$_};
    }
}


sub setProbabilities {
    my ($this,$hash) = @_;
    foreach my $k (keys(%$hash)) {
	my $pos;
	my $bn = $basenumber{$k};
	my $bgndk = $@{$this->{bgndarr}}[$bn];
	$%{$this->{probhsh}}{$k} = map {0.01 * $_} @$hash{$k};
	foreach (@{$$hash{$k}}) {
	    $@{$this->{probarr}}[$pos][$bn] = 0.01 * $_;
	    $@{$this->{lodarr}}[$pos][$bn] =
		&lg( max((0.01*$_/$bgndk), 0.000001) );
	    ++$pos;
	}
	$this->{size} = $pos;
    }

}

sub dumpDefinition {
    my ($this) = @_;

    my ($desc,$sgnl) = ($this->{description},$this->{signal});
    print
	"PSSM->new(\"$desc\",$sgnl,\n          {",
	join(", ",
	     map(sprintf($basename[$_]."=>".percent($@{$this->{bgndarr}}[$_])),
		 (0..3))),
	"},\n",
	"          {";
    my @rows;
    foreach my $bn (0..3) {
	my @entries;
	foreach my $col (0..$this->{size}-1) {
	    push
		@entries,
		sprintf "%2d", percent($@{$this->{probarr}}[$col][$bn]);
	}
	push
	    @rows,
	    "$basename[$bn]=>[". join(",", @entries) . "]";
    }
    print join(",\n           ", @rows);
    print "}\n 	);\n";

    print "# LOD scores (bits):\n";
    foreach my $bn (0..3) {
	my @entries;
	foreach my $col (0..$this->{size}-1) {
	    push
		@entries,
		sprintf "%5.1f", $@{$this->{lodarr}}[$col][$bn];
	}
	print "# $basename[$bn]=>[". join(",", @entries) . "]\n";
    }
}

sub findCrickSignals {
    my ($this,$seq,$threshold) = @_;  ## threshold is a fraction of maxScore
    $seq=scalar reverse $seq;
    $seq =~ tr/acgt/3210/;
    my @signals = $this->findWatsonSignals($seq,$threshold);
    foreach (@signals) { $$_[0] = - $$_[0]; }
    return @signals;
}



sub findWatsonSignals {
    my ($this,$seq,$threshold) = @_;  ## threshold is a fraction of maxScore
    $seq =~ tr/acgt/0123/;
    my @seq= (-1, map((int($_)), (split //, $seq)));
    my $maxscore=$this->maxScore();
    my @signals=();
    for (my $i=1; $i<=@seq-$this->{size}; $i++) {
        my $score = $this->scorePosition(\@seq,$i);
	push @signals, [$i+$this->{signal},$score]
	    unless $threshold && $score<$maxscore*$threshold;
    }
    return @signals;
}

sub scorePosition {
    my ($this,$seqref,$position) = @_;
    unless (ref($seqref)) {   ## if it's a string, not an array ref
	$seqref =~ tr/acgt/0123/;
	$seqref = [-1, map((int($_)), (split //, $seqref))];
    }
    my $score;
    return 0 if (@$seqref-$position<$this->{size});
    for (my $i=0, my $p=$position; $i<$this->{size}; $i++,$p++) {
	$score += $@{$this->{lodarr}}[$i][$$seqref[$p]];
    }
    return $score;
}

sub maxScore {
    return &sum(map {  &max(@$_) } @{$@{$_[0]->{lodarr}}});
}
	    
