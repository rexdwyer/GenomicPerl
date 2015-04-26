package IntegerInterval;
use strict;
use Util;

sub new {
    my ($this, $lo, $hi) = @_;
    my $iLo = int($lo);
    $iLo++ if $iLo < $lo;
    my $iHi = int($hi);
    $iHi-- if $iHi > $hi;
    bless [$iLo, $iHi];
}

sub intersect {
    my ($this, $other) = @_;
    my ($a,$b) = @$this;
    my ($c,$d) = @$other;
    $a>$c or $a=$c;
    $b<$d or $b=$d;
    $b>=$a or return undef;
    return $this->new($a,$b);
}

sub plus {
    my ($this, $other) = @_;
    my ($a,$b) = @$this;
    my ($c,$d) = @$other;
    return $this->new($a+$c,$b+$d);
}

sub minus {
    my ($this, $other) = @_;
    my ($a,$b) = @$this;
    my ($c,$d) = @$other;
    return $this->new($a-$d,$b-$c);
}
    
sub print {
    my ($this,@s) = @_;
    $$this[2] ||= "$$this[0]:$$this[1]";
    print "[$$this[2]]@s";
}

sub toString {
    my ($this) = @_;
    $$this[2] ||= "$$this[0]:$$this[1]";
}

1;
