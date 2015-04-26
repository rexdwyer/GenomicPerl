package SimpleReader;
require SeqReader;
@ISA = qw(SeqReader);
use strict;

sub verifyFormat {
    my ($this,$hash) = @_;
    return (bless $hash);
}

sub readSeq {
    my ($this) = @_;
    return () unless $this->{buff};
    my ($id, $seq) = ($this->{buff} =~ m/(.*):(.*)$/);
    $this->{seqId} = $id;
    $seq ||= $this->{buff};
    $seq =~ s/[^a-z]//g;
    my $fh = $this->{fh};
    $this->{buff} = <$fh>;
    chomp $this->{buff};
    return $seq;
}

1;
