#########################################
package CSPredictor;
#########################################
## This package identifies probable coding sequence in a DNA sequence based
## on dicodon (hexagram) frequencies.
#########################################
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(new predictCodingSequences);

use strict;
use Util;

#########################################
sub new
## Creates a new coding sequence predictor object.
#########################################
{
    my ($this,         ## class name CSPredictor or an instance
	$modelFile,      ## ref to file of dicodon frequencies
	) = @_;
    printf STDERR "Creating CSPredictor: $modelFile\n";

    my @lods = ((0) x 4096);
    open FREQS, "<$modelFile" or die "can't open $modelFile\n";
    my $threshold = <FREQS>;
    $threshold =~ s/[^0-9]//;	
    my $fpPat = '\d*(?:\.\d*)?(?:[Ee]-?\d+)?';
    while (<FREQS>) {
	my ($k, $kprob) = m/([acgt]{6})\s*($fpPat)/;
	$k =~ tr/acgt/0123/;
	my $hexinx = 0;
	foreach (split //, $k) {
	    $hexinx = $hexinx * 4 + int($_);
	}
	$lods[$hexinx] = ($kprob==0.0) ? -1E99 : (12+lg($kprob));
	print "lods[$hexinx] $lods[$hexinx]\n";
    }
    close FREQS;

    return bless {modelFile=>$modelFile, 
                  lodArr=>\@lods, 
                  threshold=>int($threshold)},
	         (ref($this) || $this);
}

#########################################
sub modelFile
## Returns name of file from which model was read.
#########################################
{
    my ($this) = @_;
    return $this->{modelFile};
}

#########################################
sub prepStrand
## This routine does the initial creation and filling of the dynamic programming
## arrays.  This routine is called twice, once for watson and once for crick.
## It returns a list including all three dynamic programming arrays for the
## sequence.
#########################################
{
    my ($this,    ## reference to this CSPredictor
	$seq      ## DNA sequence to be analyzed
	) = @_;
    $$seq =~ tr/acgt/0123/;
    my @lod= (-1, (split //, $$seq),0,0,0,0,0,0);    ## pad
    my @best;  $#best = $#lod-6;       ## allocate in one shot for speed.
    my @start;  $#start = $#lod-6;     ## allocate in one shot for speed.
    my $lodArr = $this->{lodArr};      ## take this out of the loop for speed.

    for (my $i=1; $i<=3; $i++) {  ## initialize three frames.
	$lod[$i] = 16*$lod[$i] + 4*$lod[$i+1] + $lod[$i+2];
	$lod[$i+3] = 16*$lod[$i+3] + 4*$lod[$i+4] + $lod[$i+5];
	$lod[$i] = $$lodArr[64*$lod[$i] + $lod[$i+3]];
	my $t = $lod[$i];
	($best[$i],$start[$i]) = $t>0?($t,$i):(0,0);
    }
    for (my $i=4; $i<@start; $i++) {  ## the "real" loop.
	$lod[$i+3] = 16*$lod[$i+3] + 4*$lod[$i+4] + $lod[$i+5];
	$lod[$i] = $$lodArr[64*$lod[$i] + $lod[$i+3]];
	my $t = $best[$i-3]+$lod[$i];
	($best[$i],$start[$i]) = $t>0?($t,($start[$i-3] || $i)):(0,0);
    }
    $#lod = $#start;      ## remove padding from @lod.
    return (\@lod, \@best, \@start, $seq, length $$seq);
}

#########################################
sub predictCodingSequences
## Returns a reference to a list of predicted coding sequences.
## Format of each list item is:
##   [ start, end, score ]
## where
##   * start and end are the indices of the first and last bases of the c.s.
##      (start>end if coding sequence is in crick.)
##   * protein is the translation of the c.s. into 1-letter codes.
##   * score is the sum of scores for the hexagrams of the c.s.
#########################################
{
    my ($this,$watson) = @_;
    (my $crick = reverse $watson) =~ tr/acgt/tgca/;
    my @watson = $this->prepStrand(\$watson);
    my @crick = $this->prepStrand(\$crick);
    my $seqlen = length $watson;
    my @pcs = $this->greedyPCS(1, $seqlen, $seqlen, \@watson, \@crick);
    return \@pcs;
}

#########################################
sub greedyPCS
## This recursive routine uses dynamic programming to identify and report the
## subsequences that maximize sum of hexagram scores.
#########################################
{
    my ($this, $lo, $hi, $size, $watson, $crick) = @_;
    return () if $lo>$hi-15;
    my ($winx,$wval) = reviewStrand($watson,$lo,$hi);
    my ($clo,$chi) = ($size+1-$hi,$size+1-$lo);
    my ($cinx,$cval) = reviewStrand($crick,$clo,$chi);
    my $bestv = ($wval>$cval) ? $wval : $cval;
    return() if $bestv < $this->{threshold};
    my ($beststart,$bestend,$bestlo,$besthi);
    if ($wval>=$cval) {
	($bestlo,$besthi) =
	    ($beststart,$bestend) = retrieveBest($watson,$winx,$lo,$hi);
    } else {
	($beststart,$bestend) = retrieveBest($crick,$cinx,$clo,$chi);
	($besthi,$bestlo) =
	    ($beststart,$bestend) = ($size+1-$beststart,$size+1-$bestend);
    }
    return ($this->greedyPCS($lo,$bestlo-1,$size,$watson,$crick),
	    [$beststart,$bestend,$bestv],
	    $this->greedyPCS($besthi+1,$hi,$size,$watson,$crick));
}
#########################################
sub reviewStrand
## This routine revises the contents of the dynamic programming arrays to
## account for the removal of a subsequence.
#########################################
{
    my ($arr,     ## reference to the list of dynamic programming arrays
	          ## constructed by prepStrand
	$lo,$hi   ## range of current interest in the sequence.
	) = @_;
    my ($lod, $best, $start,$seq,$size) = @$arr;
    $hi = min($hi,$size-2);
    foreach my $i ($lo..$lo+2) {   ## for each reading frame
	my $t = $$lod[$i];
	($$best[$i],$$start[$i]) = ($t>0)?($t,$i):(0,0);
	for (my $j=$i+3; $j<=$hi; $j+=3) {
	    ## stop as soon as entire sequence lies within range of interest
	    last unless $$start[$j] && $$start[$j]<$lo;
	    my $t = $$best[$j-3]+$$lod[$j];
	    ($$best[$j],$$start[$j]) = $t>0?($t,($$start[$j-3] || $j)):(0,0);
	}
    }
    ## find the best prediction's value and its location
    my ($besti, $bestv) = ($hi,$$best[$hi]);
    for (my $i=$hi-1; $i>=$lo; $i--) {
	($besti, $bestv) = ($i,$$best[$i]) if $$best[$i]>$bestv;
    }
    return ($besti, $bestv);
}

#########################################
sub retrieveBest
## Called once the best prediction has been determined.
## This routine tries to expand the prediction to nearest start and stop codons.
## Start is tricky because it could lie in either direction, but expansion is
## favored.  How far do you look?  The cutoff was set at 500 based on some
## experiments.
#########################################
{
    my ($arr,     ## the list of dynamic programming arrays
	$besti,   ## the index of the last codon of the best prediction
	$lo,$hi   ## current range of interest (don't exceed it!)
	) = @_;
    my ($lod, $best, $start,$seq,$size) = @$arr;

    ### Feel forward for stop codon.
    my ($begin,$end) = ($$start[$besti]+3, $besti+2);
    $hi = min($hi,$besti+500,$size-2);
    for (my $i=$besti; $i<=$hi; $i+=3) {
	($end = $i-1), last if $$lod[$i]==-1E99;
    }

    ### Feel backward for start (M) codon.
    $lo = max($lo,$begin-500);
    for (my $i=$begin-1; $i>$lo; $i-=3) {
	my $codon = substr($$seq,$i,3);
	last if $codon =~ /(302|320|300)/; ## don't back up over STOP:tag,tga,taa
	return ($i+1,$end) if $codon eq "032";  ## atg - Met
    };
    ### Feel forward for start (M) codon.
    $hi = min($begin+500,$end);
    for (my $i=$begin-1; $i<$hi; $i+=3) {
	my $codon = substr($$seq,$i,3);
	return ($i+1,$end) if $codon eq "032";    ## atg - Met
    };
    return ($begin,$end);  ## no start nearby; return initial guess.
}
