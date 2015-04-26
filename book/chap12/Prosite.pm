#########################################
package Prosite;
use strict;
use FileHandle;

#########################################
sub new
## Creates a new "attachment" to a PROSITE database.
#########################################
{
    my ($this,     ## literal "Prosite" or ref to an existing Prosite
	$fileName) ## string giving operating system's name for file
	= @_;

    my $fh = *STDIN;
    if ($fileName ne "STDIN") {
	$fh = new FileHandle "<$fileName" 
	    or die "can't open $fileName ";
    }
    bless {fh=>$fh}, (ref($this) || $this);
}


sub readMotif {
    my ($this) = @_;
    my $fh = $this->{fh};
    my $line;
    while (($line = <$fh>) && ($line !~ /^ID/)) { };
    return undef unless $line;
    my %entry;
    $entry{ID} = substr($line,5);
    while (($line = <$fh>) && ($line !~ m|^//|)) {
	$line =~ s/(..)...//;
	$entry{$1} .= $line;
    }
    return \%entry;
}

sub close {
    my ($this) = @_;
    $this->{fh}->close();
    delete @{$this}{keys %$this};
}

sub perlizePattern {
    my ($this, $pattern) = @_;
    $pattern =~ s/\.$//;
    my ($left,$right);
    my $left = '^' if ($pattern =~ s/^\<//);
    my $right = '$' if ($pattern =~ s/\>$//);
    my $newpattern = join('', map(perlizeElement($_), split('-',$pattern)));
    return $left . $newpattern . $right;
}

sub perlizeElement {
    my ($elem) = @_;
    my ($residue,$repetition) = ($elem =~ /^([^\(]*)(\(.*\))?$/);
    $repetition =~ tr|()|{}|;
    $residue =~ tr/xX/../;
    return "[^$1]$repetition" if $residue =~ /\{(.*)\}/;
    return "$residue$repetition";
}

1;
