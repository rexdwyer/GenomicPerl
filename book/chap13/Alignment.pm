#########################################
package Alignment;
#########################################
## This package ...
#########################################
use strict;
use Util;

my $SerialNumber;


my @ImmutableFields = 
    ("Read1", "Read2", ## Refs to objects describing the two reads aligned.
     "Start1","Start2",## Starting positions of alignment in each of above.
     "End1","End2",    ## Ending positions of alignment in each of above.
     "RawScore",       ## Raw Smith-Waterman score of the local alignment.
     "Path",           ## A string (MSDI)* for path through alignmt matrix.
     "MeanOffset",     ## Mean offset of aligned positions.
     "SerialNumber");

my @MutableFields = 
    ("WeightedScore", ## Smith-Waterman score weighted by column base quality
     "SupportedContig"); ## Contig whose construction this alignment supports.

foreach (@MutableFields) {
    eval "sub set$_ { \$_[0]->{$_} = \$_[1] } ;" ;
}

foreach (@MutableFields, @ImmutableFields) {
    eval "sub $_ { \$_[0]->{$_} } ;" ;
}


#########################################
sub new
#########################################
{
    my ($this,     ## literal"Alignment" or ref to existing Alignment
	$read1,    ## references to a DnaRead object
	$read2,    ## references to another DnaRead object
	$rawScore, ## raw Smith-Waterman score
	$start1,   ## starting position in $read1's sequence
	$start2,   ## starting position in $read2's sequence
	$path      ## optimal path through alignment matrix
                   ##    encoded by string of M,S,I,D.
	) = @_;

    if ($read1->SerialNumber() > $read2->SerialNumber()) {
	($read1, $read2) = ($read2, $read1);
	($start1, $start2) = ($start2, $start1);
	$path =~ tr/ID/DI/;
    }

    ## Find average offset of paired positions in the two reads.
    ## We count and compute offsets only in columns with two bases,
    ## not columns with a gap character.  Simultaneously, we are computing
    ## where the alignment ends in each read.

    my ($pos1, $pos2) = ($start1-1, $start2-1);
    my ($totalOffset, $numColumns) = (0,0);
    foreach my $column (split //, $path) {
	$pos1++ unless ($column eq "I");
	$pos2++ unless ($column eq "D");
	if ($column =~ /[MS]/) {
	    $totalOffset += ($pos2 - $pos1);
	    $numColumns++;
	}
    }
    my $meanOffset = int($totalOffset / $numColumns);

    my $alignment
	= {Read1 => $read1,
	   Read2 => $read2,
	   Start1 => $start1,
	   End1 => $pos1,
	   Start2 => $start2,
	   End2 => $pos2,
	   RawScore => $rawScore,
	   Path => $path,
	   MeanOffset => $meanOffset,
	   SerialNumber => ++$SerialNumber
	   };

    return bless $alignment, ($this || ref $this);
}


sub dump
{
    my ($this) = @_;
    $this->summarize();
    foreach ((sort @ImmutableFields),(sort @MutableFields)) {
	print "   $_ = $$this{$_}\n";
    }
}

sub summarize {
    my ($this) = @_;
    printf("Alignment number %3d between reads %3d and %3d\n",
	   $this->SerialNumber(),
	   $this->Read1()->SerialNumber(),
	   $this->Read2()->SerialNumber());
}

1;
