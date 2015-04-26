#########################################
package CSPredictor;
#########################################
## This package identifies probable coding sequence in a DNA sequence based
## on hexamer (dicodon) frequencies.
#########################################
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(new predictCodingSequences dumpDefinition translateSeq
	     translate translateSubseq);


use strict;
use Util;

my $infinity = 1E99;
my $threshold = 45;  ## We don't report predictions with smaller scores.

my @basename = qw(a c g t);
my %basenumber = (a=>0, c=>1, g=>2, t=>3);

#########################################
sub new
## Creates a new coding sequence predictor object.
#########################################
{
    my ($this,         ## class name CSPredictor or an instance
	$description,  ## brief textual description (species, etc.)
	$bgndref,      ## reference to hash of background probs for a,c,g,t
	$probref,      ## ref to hash of hexamer frequencies OR filename
	$translationHash ## reference to hash for translation, if non-standard
	) = @_;
    printf STDERR "Creating CSPredictor: $description\n";
    $this =
	bless {bgndarr=>[],
	       bgndhsh=>{},
	       probarr=>[(0) x 4096],
	       probhsh=>{},
	       lodarr=>[],
	       description=>$description,
	       codonhsh=>{}},
	(ref($this) || $this);
    $this->setBackground($bgndref);
    $probref = readHash($probref) unless ref($probref);
    $this->setProbabilities($probref);
    $this->setTranslation($translationHash);
    return $this;
}


#########################################
sub setTranslation
## Sets up the codon translation hash for this CSPredictor.
## The standard table is used unless another is explicitly provided.
## Also sets up table to translate codons given in numeric form.
#########################################
{
    my ($this,        ## reference to CSPredictor object
	$hashref      ## reference to a codon translation hash; optional
	) = @_;
    unless (scalar (keys %$hashref)) {
	$hashref =
	{ttt=>"F", tct=>"S", tat=>"Y", tgt=>"C",
	 ttc=>"F", tcc=>"S", tac=>"Y", tgc=>"C",
	 tta=>"L", tca=>"S", taa=>"#", tga=>"#",    ##  "#" is "Stop"
	 ttg=>"L", tcg=>"S", tag=>"#", tgg=>"W",

	 ctt=>"L", cct=>"P", cat=>"H", cgt=>"R",
	 ctc=>"L", ccc=>"P", cac=>"H", cgc=>"R",
	 cta=>"L", cca=>"P", caa=>"Q", cga=>"R",
	 ctg=>"L", ccg=>"P", cag=>"Q", cgg=>"R",

	 att=>"I", act=>"T", aat=>"N", agt=>"S",
	 atc=>"I", acc=>"T", aac=>"N", agc=>"S",
	 ata=>"I", aca=>"T", aaa=>"K", aga=>"R",
	 atg=>"M", acg=>"T", aag=>"K", agg=>"R",

	 gtt=>"V", gct=>"A", gat=>"D", ggt=>"G",
	 gtc=>"V", gcc=>"A", gac=>"D", ggc=>"G",
	 gta=>"V", gca=>"A", gaa=>"E", gga=>"G",
	 gtg=>"V", gcg=>"A", gag=>"E", ggg=>"G"}
    }
    foreach my $codon (keys %$hashref) {
	(my $c = $codon) =~ tr/acgt/0123/;
	$$hashref{int($c)} = $$hashref{$c} = $$hashref{$codon};
    }
    $this->{codonhsh} = $hashref;
}


#########################################
sub readHash
## Reads hexamer probabilities from a file and stores in a hash.
## Format of file is like:
##   aaaaaa 0.0003095
##   aaaaac 0.0004724
## etc.
## These files are produced in training by buildCSP.pl.
## Returns a reference to the hash.
#########################################
{
    my ($filename) = @_;    ## name of the file
    $filename = "<" . $filename;
    $filename =~ s/^\<\</\</;
    open HASH, $filename or die "can't open $filename\n";
    my @hash;
    while (<HASH>) {
	push @hash, m/([acgt]{6})\s*(\d*(?:\.\d*)?)/;
    }
    my %hash = @hash;
    return \%hash;
}


#########################################
sub setBackground
## Sets up the background probability hash and array.
## Translates percents into fractions.
#########################################
{
    my ($this,    ## reference to CSPredictor object
	$hash)    ## ref to hash table with probability for each of 4 bases.
	= @_;
    $hash = {a=>25,c=>25,g=>25,t=>25} unless ref($hash);
    my $ba = $this->{bgndarr};
    my $bh = $this->{bgndhsh};
    foreach (keys(%$hash)) {
	$%{$bh}{$_} = $@{$ba}[$basenumber{$_}] = 0.01 * $$hash{$_};
    }
}



#########################################
sub setProbabilities
## Uses the hexamer probabilities in input hash to set up arrays of lod scores
## and probabilities.
#########################################
{
    my ($this,$hash) = @_;
    $this->{probhsh} = $hash;
    my $pa = $this->{probarr};          ## These three
    my @ba = @{$@{$this->{bgndarr}}};    ## lines make things
    my $la = $this->{lodarr};             ## considerably faster...
    foreach my $k (keys(%$hash)) {
	my $kprob = $$hash{$k};
	$k =~ tr/acgt/0123/;
	my ($hexinx, $hexbgnd) = (0,1);
	foreach (split //, $k) {
	    $hexinx = $hexinx * 4 + int($_);
	    $hexbgnd *= $ba[int($_)];
	}
	$@{$pa}[$hexinx] = $kprob;
	$@{$la}[$hexinx] = ($kprob==0.0)?-$infinity:&lg($kprob/$hexbgnd);
    }
}
#########################################
sub dumpDefinition
## Writes a description of this CSPredictor to a file for later execution.
#########################################
{
    my ($this) = @_;

    my $desc = $this->{description};
    print
	"CSPredictor->new(\"$desc\",\n          {",
	join(", ",
	     map(sprintf($basename[$_]."=>".percent($@{$this->{bgndarr}}[$_])),
		 (0..3))),
	"},\n";

    my $pa = $this->{probarr};
    my $ba = $this->{bgndarr};
    my $la = $this->{lodarr};
    printf("          {aaaaaa=>%9.7f, #lod score=%+5.1f\n",
	   $@{$pa}[0], $@{$la}[0]);
    foreach (1..4094) {
	my $hexamer;
	my $bn=4*$_;
	foreach (1..6) { $hexamer = $basename[($bn/=4) % 4] . $hexamer; }
	printf("           $hexamer=>%9.7f, #lod score=%+5.1f\n",
	       $@{$pa}[$_], $@{$la}[$_]);
    }
    printf("           tttttt=>%9.7f} #lod score=%+5.1f\n     );\n\n",
	   $@{$pa}[4095], $@{$la}[4095]);
}

#########################################
sub predictCodingSequences
## Returns a reference to a list of predicted coding sequences.
## Format of each list item is:
##   [ start, end, score, protein ]
## where
##   * start and end are the indices of the first and last bases of the c.s.
##      (start>end if coding sequence is in crick.)
##   * protein is the translation of the c.s. into 1-letter codes.
##   * score is the sum of scores for the hexamers of the c.s.
#########################################
{
    my ($this,$watson) = @_;
    (my $crick = reverse $watson) =~ tr/acgt/tgca/;
    my @watson = $this->prepStrand(\$watson);
    my @crick = $this->prepStrand(\$crick);

    my $seqlen = length $watson;
    my @pcs = greedyPcs(1, $seqlen, $seqlen, \@watson, \@crick);
    foreach (@pcs) {
	push @$_, $this->translateSubseq($$_[0],$$_[1],\$watson,\$crick);
    }
    return \@pcs;
}

#########################################
sub greedyPcs
## This recursive routine uses dynamic programming to identify and report the
## subsequences that maximize sum of hexamer scores.
#########################################
{
    my ($lo, $hi, $size, $watson, $crick) = @_;
    return (()) if $lo>$hi-15;
    my ($winx,$wval) = reviewStrand($watson,$lo,$hi);
    my ($clo,$chi) = ($size+1-$hi,$size+1-$lo);
    my ($cinx,$cval) = reviewStrand($crick,$clo,$chi);
    my $bestv = $wval>$cval?$wval:$cval;
    return (()) unless $bestv >= $threshold;
    my ($beststart,$bestend,$bestlo,$besthi);
    if ($wval>=$cval) {
	($bestlo,$besthi) =
	    ($beststart,$bestend) = retrieveBest($watson,$winx,$lo,$hi);
    } else {
	($beststart,$bestend) = retrieveBest($crick,$cinx,$clo,$chi);
	($besthi,$bestlo) =
	    ($beststart,$bestend) = ($size+1-$beststart,$size+1-$bestend);
    }
    return (greedyPcs($lo,$bestlo-1,$size,$watson,$crick),
	    [$beststart,$bestend,$bestv],
	    greedyPcs($besthi+1,$hi,$size,$watson,$crick));
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
    my $lodarr = $this->{lodarr};      ## take this out of the loop for speed.

    for (my $i=1; $i<=3; $i++) {  ## initialize three frames.
	$lod[$i] = 16*$lod[$i] + 4*$lod[$i+1] + $lod[$i+2];
	$lod[$i+3] = 16*$lod[$i+3] + 4*$lod[$i+4] + $lod[$i+5];
	$lod[$i] = $@{$lodarr}[64*$lod[$i] + $lod[$i+3]];
	my $t = $lod[$i];
	($best[$i],$start[$i]) = $t>0?($t,$i):(0,0);
    }
    for (my $i=4; $i<@start; $i++) {  ## the "real" loop.
	$lod[$i+3] = 16*$lod[$i+3] + 4*$lod[$i+4] + $lod[$i+5];
	$lod[$i] = $@{$lodarr}[64*$lod[$i] + $lod[$i+3]];
	my $t = $best[$i-3]+$lod[$i];
	($best[$i],$start[$i]) = $t>0?($t,($start[$i-3] || $i)):(0,0);
    }
    $#lod = $#start;      ## remove padding from @lod.
    return (\@lod, \@best, \@start, $seq, length $$seq);
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
## experiments with known genes in Ashbya Chromosome I.
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
	($end = $i-1), last if $$lod[$i]==-$infinity;
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



#########################################
sub translate
## Uses the codon translation table to translate a codon.
## All arguments are concatenated to get the codon, so it may be given as
## a number, a string of 3 characters, 3 single characters, etc.
#########################################
{
    my ($this,@codon) = @_;
    my $codon = join("",@codon);
    return (${$this->{codonhsh}}{$codon} || "Z");
}

#########################################
sub translateSeq
## Uses the codon translation table to translate a nucleotide string.
#########################################
{
    my ($this,$seq) = @_;
    my $seq1 = $seq;
    my $check = length $seq;
    my $hash = $this->{codonhsh};
    my $polypeptide;
    while ($seq =~ s/^(...)//) {
	$polypeptide .= (${$hash}{$1} || "Z");
    }
##    ($check == 3*length $polypeptide) or
##       die "translation error\n$seq1\n$polypeptide";
    return $polypeptide;
}
#########################################
sub translateSubseq
## Uses the codon translation table to translate a subsequence
## All arguments are concatenated to get the codon, so it may be given as
## a number, a string of 3 characters, 3 single characters, etc.
#########################################
{
    my ($this,$start,$end,$refWatson,$refCrick) = @_;
    my $subseq;
    if ($start<$end) {   ## watson
	$subseq = substr($$refWatson, $start-1, $end-$start+1);
    } else {
	unless (ref($refCrick)) {
	    (my $crick = reverse $$refWatson) =~ tr/acgt0123/tgca3210/;
	    $refCrick = \$crick;
	}
	$subseq = substr($$refCrick, -$start, $start-$end+1);
    }
    return $this->translateSeq($subseq);
}
