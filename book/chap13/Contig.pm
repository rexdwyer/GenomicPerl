#########################################
package Contig;
#########################################
## This package ...
#########################################
use strict;

my $SerialNumber;

my @ImmutableFields = 
    ("SerialNumber");  ## Records order in which Contig objects are created.

my @MutableFields = 
    ("Bases",       ## The consensus sequence of the contig; a string.
     "Quals",       ## Ref to list of qualities caller assigned to each base.
     "Reads",       ## Ref to list of included reads.
     "AlignedPairs",## Ref to list of alignments pertaining to this contig.
     "Length",      ## Estimated length of contig.
     "FirstStart",  ## Where leftmost read starts, wrt contig's 0. (Often <0).
     "LastEnd",     ## Where rightmost read end, relative to contig's 0.
     "Score");      ## Sum of nucleotide qualities for this contig.

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
    my ($this,     ## literal "Contig" or ref to an existing Contig.
	$reads,
	$numReads,
	$firstStart,
	$lastEnd,
	) = @_;

    return bless {
	Reads => $reads,
	NumReads => $numReads,
	FirstStart => $firstStart,
	LastEnd => $lastEnd,
	AlignedPairs => [ ],
	SerialNumber => $SerialNumber++}, ($this || ref $this);
}


#########################################
sub dump
#########################################
{
    my ($this) = @_;  ## ref to an existing Contig.
    print "****************************************************************\n";
    print "Contig $$this{SerialNumber} at $this\n";
    print "****************************************************************\n";
    foreach ((sort @ImmutableFields),(sort @MutableFields)) {
	print "   $_ = $$this{$_}\n";
    }
    print "This Contig contains the following reads:\n";
    foreach my $read (@{$this->Reads()}) { $read->summarize(); }
}

1;
