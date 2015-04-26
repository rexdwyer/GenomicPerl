package FastqReader;
require SeqReader;
@ISA = qw(SeqReader);
use strict;

sub verifyFormat {
    my ($this, $hash) = @_;
    return (bless $hash) if $$hash{buff} =~ m/^[qQ]\>/;
    return "";
}

#########################################
sub readSeq {
## Reads, returns sequence and id from fasta file.
#########################################
    my ($this) = @_;
    return () unless $this->{buff};
    my $id = $this->{buff};
    $id =~ s/^\>\s*//;
    $this->{seqId} = $id;
    my ($seq, $tbuff);
    my $fh = $this->{fh};
    while (($tbuff = <$fh>) && ($tbuff !~ /^[Qq]\>/)) {
	chomp $tbuff;
	$tbuff = lc $tbuff;
	$tbuff =~ s/[^a-z0-9]//g;
	$seq .= $tbuff;
    }
    chomp $tbuff;
    my $quals = $seq;
    $seq =~ s/[0-9]//g;
    $quals =~ s/[a-z]//g;
    warn "Incorrect number of qualities for sequence $id in $$this{fileName}"
	if (length $quals > 0 && (length $quals != length $seq));
##    warn "No qualities for sequence $id in $$this{fileName}; using defaults"
##	if length $quals == 0;
    $quals = "1" x (length $seq) if (length $quals != length $seq);
    $this->{seqQuality} = [map {int $_} split(//, $quals)];
    $this->{buff} = $tbuff;
    return $seq;
}


sub seqQuality {
    my ($this) = @_;
    $this->{seqQuality};
}

1;
