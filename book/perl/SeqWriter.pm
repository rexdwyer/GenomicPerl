#########################################
package SeqWriter;
#########################################
## This package allows DNA or protein sequences to be transparently written to
## files in fasta, GenBank, or one-sequence-per-line format.
#########################################
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(new writeSeq prefix fileType close fileName writeSneak);

use strict;
use FileHandle;

#########################################
sub new
## Create a new output file and an object to manipulate it.
#########################################
{
    my ($this,     ## literal "SeqWriter" or ref to an existing SeqWriter
	$filetype, ## type of file to write: "fasta", "GenBank", or "OneLiners"
	$filename)  ## operating system's name for new file (type is appended)
	= @_;
    my $prefix = ($filename || "STDOUT");
    $prefix =~ s/^[\>\&]//;
    $prefix = $prefix || "STDOUT";
    my $fh = *STDOUT;
    ($fh = new FileHandle ">$filename.$filetype" or
     die "can't open $filename\n") if $filename ne "STDOUT";
    ### Save the correct writing routine for this type:
    my $writeSeq = \&writeSneak;
    $writeSeq = \&writeFasta if lc($filetype) eq "fasta";
    $writeSeq = \&writeGenBank if lc($filetype) eq "genbank";
    warn "Writing to $filetype file $filename.$filetype...\n";

    return
	bless {prefix=>$prefix,        ## prefix of filename
	       fh=>$fh,                ## filehandle
	       writeSeq=>$writeSeq,  ## routine for writing
	       fileType=>$filetype},  ## file type: fasta, genbank, etc.
	   (ref($this) || $this);
}

#########################################
sub close
## Closes the file to further output.
#########################################
{
    my ($this) = @_;
    $this->{fh}->close();
}

#########################################
sub writeSeq {
## Invokes the correct write routine to write a sequence entry.
#########################################
    return &{$_[0]->{writeSeq}}(@_);
}

#########################################
sub compress {
## Squeezes out extra whitespace
#########################################
    $_[0] =~ s/\s+/ /g;
    return $_[0];
}

#########################################
sub prefix {
## Returns prefix of file's operating system name.
#########################################
    return $_[0]->{prefix};
}

#########################################
sub fileType {
## Returns file's type: fasta, GenBank, OneLiners.
#########################################
    return $_[0]->{fileType};
}

#########################################
sub fileName {
## Returns file's entire operating system name.
#########################################
    my ($this) = @_;
    return ($this->{prefix} . "." . $this->{fileType});
}


#########################################
sub writeSneak {
## Writes sequence to a OneLiners file.
## Can also be used to sneak HTML, etc., into other files.
#########################################
    my ($this,$seq,$id,$notes) = @_;
    $this->{fh}->print("$seq\n");
}


#########################################
sub writeFasta {
## Writes sequence and id to a fasta file.
## Sequence is divided into lines and they are numbered.
#########################################
    my ($this,$seq,$id,$notes) = @_;
    $this->{fh}->print(">$id\n");
    $seq =~ s/(.{80})/$1\n/g;
    $this->{fh}->print("$seq\n");
}

#########################################
sub writeGenBank {
## Writes sequence, id, and annotations to a GenBank file.
## Sequence is divided into lines and they are numbered.
## Annotations must already be correctly formatted; they are just printed.
#########################################
    my ($this,$seq,$id,$notes) = @_;
    $this->{fh}->print("LOCUS $id\n");
    $this->{fh}->print($notes);
    $this->{fh}->print("ORIGIN\n");
    my $index = 1;
    my $line;
    while ($seq =~ s/(.{60})//) {
	$line = $1;
	$line =~ s/(.{10})/ $1/g;
	$this->{fh}->print(substr("         ".$index,-10), "$line\n");
	$index += 60;
    }
    if ($seq) {
	$seq =~ s/(.{10})/$1 /g;
	$this->{fh}->print(substr("         ".$index,-10), " $seq\n");
    }
    $this->{fh}->print("//\n");
}
