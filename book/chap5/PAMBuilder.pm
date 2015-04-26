package PAMBuilder;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(printMatrix perlFormatMatrix collectCounts 
	     addPseudocounts freqsFromCounts ratesFromFreqs
	     rescaleRates freqsFromRates lodsFromFreqs
	     averageMutationRate multiplyRates);

use strict;
use Util;
my @residues = split //,"acdefghiklmnpqrstvwy";

sub checkType {
    my ($M,$type) = @_;
    $$M{type} =~ /($type)/;
    $1 or die "Matrix is of type $$M{type}; $type is required.";
}

sub printMatrix {
    my ($M,$header) = @_; 
    my %format = (frequencies=>" %4.2f",'mutation rates'=>" %4.2f",
		  counts=>"%5d", 'lod scores'=>"%5d");
    my $format = $format{$$M{type}};
    print "\#\#\# $header\n\#  bkgnd";
    foreach (@residues) { printf "%5s", $_; }
    foreach my $r1 (@residues) {
	printf "\n# $r1$format", $$M{$r1};
	foreach my $r2 (@residues) {
	    printf $format, $$M{$r1.$r2};
	}
    }
    print "\n\n";
}

sub perlFormatMatrix {
    my ($M,$name) = @_;
    print "$name=\n";
    my $s = join(",", map("$_=>'$$M{$_}'", (keys %$M)));
    print "(";
    while ($s =~ s/^(.{0,70}\,)//) { print " $1\n";}
    print " $s);\n";
}

sub collectCounts {
    my %nu = (type=>'counts');
    while (chomp(my $seq1 = <STDIN>)) {
	chomp(my $seq2 = <STDIN>);    
	foreach my $i (0..((length $seq1)-1)) {
	    my ($r1,$r2) = (substr($seq1,$i,1), substr($seq2,$i,1));
	    next if ($r1.$r2) =~ "-";
	    $nu{$r1}++; $nu{$r2}++;
	    $nu{$r1.$r2}++; $nu{$r2.$r1}++;
	}
    }
    return \%nu;
}

sub addPseudocounts {
    my ($M, $pseudo) = @_;   ## pseudo optional; default 1.
    $pseudo ||= 1;
    checkType($M, 'counts');
    my %nu = %$M;
    foreach my $r1 (@residues) {
	foreach my $r2 (@residues) {
	    $nu{$r1} += $pseudo ; $nu{$r2} += $pseudo;
	    $nu{$r1.$r2} += $pseudo; $nu{$r2.$r1} += $pseudo;
	    last if ($r1 eq $r2);
	}
    }
    return \%nu;
}

sub freqsFromCounts {
    my ($M) = @_;
    checkType($M, 'counts');
    my %nu = (type=>'frequencies');
    my $total=0;
    foreach my $r (@residues) { $total += $$M{$r}; }
    foreach my $r (@residues) { $nu{$r} = $$M{$r} / $total; }
    $total=0;
    foreach my $r1 (@residues) {
	foreach my $r2 (@residues) { 
	    $total += $$M{$r1.$r2};
	}
    }
    foreach my $r1 (@residues) {
	foreach my $r2 (@residues) { 
	    $nu{$r1.$r2} = $$M{$r1.$r2}/$total;
	}
    }
    return \%nu;
}

sub ratesFromFreqs {
    my ($M) = @_;
    checkType($M, 'frequencies');
    my %nu = (type=>'mutation rates');
    foreach my $r1 (@residues) {
	my $q = $nu{$r1} = $$M{$r1};
	foreach my $r2 (@residues) { $nu{$r1.$r2} =  $$M{$r1.$r2} / $q; }
    }
    return \%nu;
}

sub freqsFromRates {
    my ($M) = @_;
    checkType($M, 'mutation rates');
    my %nu = (type=>'frequencies');
    foreach my $r1 (@residues) {
	my $q = $nu{$r1} = $$M{$r1};
	foreach my $r2 (@residues) { $nu{$r1.$r2} =  $$M{$r1.$r2} * $q; }
    }
    return \%nu;
}

sub averageMutationRate {
    my ($M) = @_;
    my $type = checkType($M, 'mutation rates|frequencies');
    my $rate = 1;
    if ($type eq "frequencies") {
	foreach my $r (@residues) { $rate -= $$M{$r.$r}; }
    } else {
	foreach my $r (@residues) { $rate -= $$M{$r} * $$M{$r.$r}};
    }
    return $rate;
}

sub rescaleRates {
    my ($M,$newrate) = @_;
    checkType($M, 'mutation rates');
    my %nu = (type=>'mutation rates');
    my $factor = $newrate/averageMutationRate($M);
    foreach my $r1 (@residues) {
	foreach my $r2 (@residues) { 
	    $nu{$r1.$r2} = $$M{$r1.$r2} * $factor;
	}
	$nu{$r1.$r1} += (1-$factor);
	$nu{$r1} = $$M{$r1};
    }
    return \%nu;
}

sub multiplyRates {
    my ($M1,$M2) = @_;
    checkType($M1, 'mutation rates');
    checkType($M2, 'mutation rates');
    my %nu = (type=>'mutation rates');
    foreach my $r1 (@residues) {
	foreach my $r2 (@residues) { 
	    my $prod;
	    foreach my $r (@residues) {
		$prod += $$M1{$r1.$r} * $$M2{$r.$r2};
	    }
	    $nu{$r1.$r2} = $prod;
	}
	$nu{$r1} = $$M1{$r1};
    }
    return \%nu;
}

sub lodsFromFreqs {
    my ($M) = @_;
    checkType($M, 'frequencies');
    my %nu = (type=>'lod scores');
    foreach my $r1 (@residues) {
	$nu{$r1} = round(2*lg($$M{$r1}));
	foreach my $r2 (@residues) { 
	    $nu{$r1.$r2} =
		round(2*lg($$M{$r1.$r2} / ($$M{$r1}*$$M{$r2})));
	}
    }
    return \%nu;
}

