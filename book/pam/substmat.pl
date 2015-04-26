#!/usr/bin/perl -I /home/dwyer/book/perl

use strict;
use Util;

my %pair;
my %bkgnd;

sub printMatrix {
    my ($format) = @_;     # optional printf-style entry for matrix entries.
    $format = $format || "%5d";   # default format is 5-digit integer.
    my @residues = split //,"cstpagndeqhrkmilvfyw";   # list of 20 amino acids
    foreach my $r1 (@residues) {
	print "# $r1";
	foreach my $r2 (@residues) {
	    printf $format, ($pair{$r1.$r2} || $pair{$r2.$r1});
	    last if $r1 eq $r2;      # lower triangular matrix
	}
	print "\n";
    }
    print "\#  ";
    foreach (@residues) { printf "%5s", $_; }
    print "\n\#  ";
    foreach (@residues) { printf $format, $bkgnd{$_}; }
    print "\n\n";
}



while (my $score = <STDIN>) {
    my $id1 = <STDIN>;
    my $id2 = <STDIN>;
    chomp(my $seq1 = <STDIN>);
    chomp(my $seq2 = <STDIN>);    
    foreach my $i (0..((length $seq1)-1)) {
	my ($r1,$r2) = (substr($seq1,$i,1), substr($seq2,$i,1));
	next if ($r1.$r2) =~ "-";
	$bkgnd{$r1}++; $bkgnd{$r2}++;
	($r1 ge $r2) ?  $pair{$r1.$r2}++ : $pair{$r2.$r1}++ ;
    }
}

print "\#\#\# Raw Frequencies:\n";
printMatrix();

my @residues = split //,"acdefghiklmnpqrstvwy";   # list of 20 amino acids
foreach my $r1 (@residues) {
    foreach my $r2 (@residues) {
	$pair{$r1.$r2}++;
	last if $r1 eq $r2;      # lower triangular matrix
    }
}
print "\#\#\# Frequencies plus Pseudocounts:\n";
printMatrix();

### Convert %bkgnd from counts to probabilities
my $bcount;
foreach my $r (keys %bkgnd) { $bcount += $bkgnd{$r}; }
foreach my $r (keys %bkgnd) { $bkgnd{$r} /= $bcount; }
### Convert %bkgnd from counts to probabilities
my $pcount;
foreach my $rr (keys %pair) { $pcount += $pair{$rr}; }
foreach my $rr (keys %pair) { $pair{$rr} /= $pcount; }

### Convert %bkgnd and %pair to percentages for display; display; convert back
foreach my $r (keys %bkgnd) { $bkgnd{$r} *= 100; }
foreach my $rr (keys %pair) { $pair{$rr} *= 100; }
print "\#\#\# Probabilities (Percent)\n";
printMatrix(" %4.1f");   # display to three digits, e.g., 17.4%
foreach my $r (keys %bkgnd) { $bkgnd{$r} /= 100; }
foreach my $rr (keys %pair) { $pair{$rr} /= 100; }

### Convert %pair to lod-scores in half-bits; display.
foreach my $rr (keys %pair) {
    my ($r1,$r2) = split //, $rr;
    $pair{$rr} = round(2*lg($pair{$rr}/(2*$bkgnd{$r1}*$bkgnd{$r2})));
    $pair{$rr} += 2 if $r1 eq $r2;   # cancel 2 in denominator above;
}
print "\#\#\# LOD scores in half-bits\n";
printMatrix();

### Now print as source for a Perl hash table.
print "my %lod = (\n";
my $t;
print join(",\n", 
	   map((($_ eq ($t=reverse $_)) ?
		"    $_ => $pair{$_}" :
		"    $_ => $pair{$_}, $t => $pair{$_}"),
	       (sort keys %pair)));
print ");\n";
