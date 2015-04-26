#########################################
package SeqReader;
#########################################
## This package allows DNA or protein sequences to be transparently read from
## files in fasta, GenBank, or one-sequence-per-line format.
#########################################
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(new readSeq prefix fileType close);

use strict;
use FileHandle;

#########################################
sub new
## Creates a new "attachment" to the named sequence file.
## Determine the file type - fasta, GenBank, or one-per-line - for future reads.
#########################################
{
    my ($this,     ## literal"SeqReader" or ref to an existing SeqReader
	$filename) ## string giving operating system's name for file
	= @_;

    my $prefix = ($filename || "STDIN");
    $prefix =~ s/^[\<\&]//;        ## remove "<" or "&" from start of name
    $prefix = $prefix || "STDIN";
    $prefix =~ s/\.[^\.]*$//;      ## remove one ".xxx" from end.
    my $fh = *STDIN;
    ($fh = new FileHandle "<$filename" or die "can't open $filename\n")
	if $filename ne "STDIN";
    my $buff = <$fh>;              ## read first line to determine type
    chomp $buff;
    ### Save the correct reading routine according to type:
    my $readSeq = \&readOnePerLine;
    $readSeq = \&readFasta if $buff =~ m/^\>/;
    $readSeq = \&readGenBank if $buff =~ m/^LOCUS/;
    ### Save the type as a string, too.
    my $fileType = "OneLiners";
    $fileType = "fasta" if $buff =~ m/^\>/;
    $fileType = "GenBank" if $buff =~ m/^LOCUS/;
    warn "Reading from $fileType file $filename...\n";

    return
	bless {prefix=>$prefix,        ## save prefix of name (to name outputs)
	       fh=>$fh,                ## save filehandle
	       buff=>$buff,            ## save line read for next read.
	       readSeq=>$readSeq,    ## save read routine
	       fileType=>$fileType}, ## save file type
	   (ref($this) || $this);
}

#########################################
sub readSeq {
## Invokes the correct read routine for this file
## Each returns a list of three items: (sequence, id, annotations).
## However, id and annotations may be undefined for some files.
#########################################
    return &{$_[0]->{readSeq}}(@_);
}

#########################################
sub compress {
## Squeezes out extra whitespace.
#########################################
    $_[0] =~ s/\s+/ /g;
    return $_[0];
}

#########################################
sub prefix {
## Returns prefix of file name
#########################################
    return $_[0]->{prefix};
}

#########################################
sub fileType {
## Returns type of file (fasta, GenBank, one-per-line)
#########################################
    return $_[0]->{fileType};
}


#########################################
sub readOnePerLine {
## Reads, returns sequence from one-sequence-per-line file.
#########################################
    my ($this) = @_;
    return () unless $this->{buff};
    my $seq = lc $this->{buff};
    $seq =~ s/[^a-z]//g;
    my $fh = $this->{fh};
    $this->{buff} = <$fh>;
    chomp $this->{buff};
    return ($seq,undef,undef);
}

#########################################
sub readFasta {
## Reads, returns sequence and id from fasta file.
#########################################
    my ($this) = @_;
    return () unless $this->{buff};
    my $id = compress($this->{buff});
    $id =~ s/^\>\s*//;
    my ($seq, $tbuff);
    my $fh = $this->{fh};
    while (($tbuff = <$fh>) && ($tbuff !~ /^\>/)) {
	chomp $tbuff;
	$tbuff = lc $tbuff;
	$tbuff =~ s/[^a-z]//g;
	$seq .= $tbuff;
    }
    chomp $tbuff;
    $this->{buff} = $tbuff;
    return ($seq,$id,undef);
}

#########################################
sub readGenBank {
## Reads, returns sequence, id and annotations from GenBank file.
#########################################
    my ($this) = @_;
    return () unless $this->{buff};
    my $id = compress($this->{buff});
    $id =~ s/^LOCUS\s\s*//;
    my ($seq, $notes, $tbuff);
    my $fh = $this->{fh};
    while (($tbuff = <$fh>) && ($tbuff !~ /^ORIGIN/)) {
	$notes .= $tbuff;
    }
    $notes .= $tbuff;
    while (($tbuff = <$fh>) && ($tbuff !~ /^LOCUS/)) {
	chomp $tbuff;
	$tbuff = lc $tbuff;
	$tbuff =~ s/[^a-z]//g;
	$seq .= $tbuff;
    }
    $this->{buff} = $tbuff;
    return ($seq,$id,$notes);
}

#########################################
sub close {
## Closes the file.
#########################################
    my ($this) = @_;
    $this->{fh}->close();
}
