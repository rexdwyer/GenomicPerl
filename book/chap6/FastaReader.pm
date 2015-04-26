package FastaReader;
require SeqReader;
@ISA = qw(SeqReader);
use strict;

sub verifyFormat {
    my ($this, $hash) = @_;
    return (bless $hash) if $$hash{buff} =~ m/^\>/;
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
    while (($tbuff = <$fh>) && ($tbuff !~ /^\>/)) {
	chomp $tbuff;
	$tbuff = lc $tbuff;
	$tbuff =~ s/[^a-z]//g;
	$seq .= $tbuff;
    }
    chomp $tbuff;
    $this->{buff} = $tbuff;
    return $seq;
}



