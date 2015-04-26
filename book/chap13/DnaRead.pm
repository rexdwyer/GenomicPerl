package DnaRead;

use strict;
use Util;

my @AllReads;

my @ImmutableFields = 
    ("Bases",       ## The nucleotide sequence of the read; a string.
     "Length",      ## Number of nucleotides in the read.
     "Name",        ## The name of the sequence; a string.
     "OrigQual",    ## Ref to list of qualities caller assigned to each base.
     "SerialNumber");  ## Records index of object in @AllReads.

my @MutableFields = 
    ("Alignments",  ## Ref to list of alignments to this read.
     "AdjQual",     ## Ref to list of current quality estimates.
     "Contig",      ## Ref to contig to which read is currently assigned.
     "ContigStart", ## Approx. interval of contig covered by this read.
     "ContigEnd");   ## 

foreach (@MutableFields) {
    eval "sub set$_ { \$_[0]->{$_} = \$_[1] } ;" ;
}

foreach (@MutableFields, @ImmutableFields) {
    eval "sub $_ { \$_[0]->{$_} } ;" ;
}


#########################################
sub new {
    my ($this,     ## literal"DnaRead" or ref to an existing DnaRead
	$name,     ## identifier of this sequences
	$bases,    ## the bases
	$qualities ## qualities corresponding to bases
	) = @_;

    my $this = bless {Bases => $bases, 
		  Length => length($bases),
		  Name => $name,
		  OrigQual => $qualities,
		  AdjQual => \@$qualities,  ## Makes a copy of the array.
		  SerialNumber => scalar(@AllReads)}, ($this || ref $this);
    push @AllReads, $this;
    return $this;
}

#########################################
sub GetBySerial {
    my ($this, $num) = @_;
    return $AllReads[$num];
}

#########################################
sub dump {
    my ($this) = @_;
    
    print "================ Dump of $this =======================\n";
    foreach ((sort @ImmutableFields),(sort @MutableFields)) {
	print "   $_ = $$this{$_}\n";
    }
    print "   OrigQual = ", join('',@{$this->{OrigQual}}), "\n";
    print "    AdjQual = ", join('',@{$this->{AdjQual}}), "\n";

    print "=======================================================\n";
}

#########################################
sub summarize {
    my ($this) = @_;
    printf("**> read %3d at (%6d,%6d) in contig %3d ( %s )\n",
	   @$this{"SerialNumber","ContigStart","ContigEnd"},
	   $this->Contig()->SerialNumber(), 
	   $this->Name());
}

1;

