package GenBankReader;
require SeqReader;
@ISA = qw(SeqReader);
use strict;

sub verifyFormat {
    my ($this, $hash) = @_;
    return (bless $hash) if $$hash{buff} =~ m/^LOCUS/;
    return "";
}

sub readSeq {
    my ($this) = @_;
    return () unless $this->{buff};
    my $id = $this->{buff};
    $id =~ s/^LOCUS\s\s*//;
    $this->{seqId} = $id;
    my ($seq, $notes, $tbuff);
    my $fh = $this->{fh};
    while (($tbuff = <$fh>) && ($tbuff !~ /^ORIGIN/)) {
	$notes .= $tbuff;
    }
    $notes .= $tbuff;
    $this->{seqNotes} = $id;
    while (($tbuff = <$fh>) && ($tbuff !~ /^LOCUS/)) {
	chomp $tbuff;
	$tbuff = lc $tbuff;
	$tbuff =~ s/[^a-z]//g;
	$seq .= $tbuff;
    }
    $this->{buff} = $tbuff;
    return $seq;
}

sub complement {
    my ($s) = @_;
    $s =~ tr/acgtACGT/tgcaTGCA/;
    return reverse $s;
}

sub ident {
    my ($s) = @_;
    return $s;
}

sub extract_feature {
    my ($this, $location, $s) = @_;
    print @_, "\n";
    ## Allow $s to be either a string or a reference to a string.
    unless (ref $s) {my $t = $s; $s = \$t; }
    $location = lc $location;
    $location =~ s/,/\./g;
    $location =~ s/join/ident/g;
    $location =~ s/(\d+)/N$1/g;
    $location =~ s/N(\d+)\.\.N(\d+)/substr(\$\$s,$1-1,$2+1-$1)/g;
    $location =~ s/N(\d+)/substr(\$\$s,$1-1,1)/g;
    return eval($location);
}





