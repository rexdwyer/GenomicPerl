package Recurrence;
use strict;

sub new
{
    my ($this,@terms) = @_;
    my %coefficients;
    while (@terms) {
	(my $d, my $c, @terms) = @terms;
	$coefficients{0-$d} += $c;
    }
    my @coefficients 
	= map([$_,$coefficients{$_}], sort {$a<=>$b} keys %coefficients);
    my @polynomial;
    my $offset = $coefficients[0][0];
    $polynomial[-$offset] = -1.0;
    foreach (@coefficients) {
	my ($d,$c) = @$_;
	$polynomial[$d-$offset] += $c;
    }
    my $hash = {coefficients => \@coefficients,
		polynomial => \@polynomial,
		memory =>  { }};
    bless $hash, (ref($this) || $this);
}

sub printPolynomial {
    my ($this) = @_;
    my $poly = $this->{polynomial};
    print "$$poly[-1]*x^$#$poly +";
    for (my $i=$#$poly-1; $i>=0; $i--) {
	next if $$poly[$i]==0;
	print " $$poly[$i]";
	if ($i>1) { print "*x^$i +"; }
	elsif ($i==1)  { print "*x +"; }
	else { print " "; }
    }
    print "= 0\n";
}

sub printRecurrence {
    my ($this) = @_;
    print "\nRecurrence Relation is\n";
    my $leader =  "P[s,m] = ";
    foreach (@{$this->{coefficients}}) {
	my ($d,$c) = @$_;
	print "$leader $c * P[s";
	printf "%+0d", $d if $d!=0;
	print ",m-1] \n";
	$leader =  "        +";
    }
}

sub polyRoot {  ## not a general routine for solving polynomials!
    my ($this) = @_;
    my $poly = $this->{polynomial};
    my ($neg,$pos) = (1,0);
    while ($neg-$pos > 1E-15) {
	my $x = ($neg+$pos)/2.0;
	my $value = 0;
	my $power = 1;
	foreach my $coeff (@$poly) {
	    $value += $coeff * $power;
	    $power *= $x;
	}
	if ($value>0) { $pos = $x; } else { $neg = $x; }
    }
    return ($neg+$pos)/2.0;
}

sub value {
    my ($this,$s,$m) = @_;
    my $mem = $this->{memory};
    return $$mem{$s,$m} if defined($$mem{$s,$m});
    return 1 if $s<=0 && $m>=0;
    return 0 if $m<=0;
    my $val;
    foreach (@{$this->{coefficients}}) {
	my ($d,$c) = @$_;
	$val += $c * $this->value($s+$d,$m-1)
    }
    return $$mem{$s,$m}=$val;
}

1;
