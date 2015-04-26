#!/usr/bin/perl -I /home/dwyer/book/perl
use strict;
use Util;
use Recurrence;

my ($p,$s) = readScoreScheme();
my $avg = averagePairScore($p,$s);
print "\nAverage Score per Pair is $avg\n";

my @terms;
foreach my $i (0..$#$p) {
    foreach my $j (0..$#$p) {
	push @terms, $$s[$i][$j], $$p[$i]*$$p[$j];
    }
}
my $Psm = Recurrence->new(@terms);
$Psm->printRecurrence();
print "\nCharacteristic Polynomial is\n";
$Psm->printPolynomial();
my $root = $Psm->polyRoot();
print "\nRelevant root is $root\n";

print "\nTable of Expected-s[m] and P[s,m]\n";
printTable($Psm, [1..20], [1..22]);
print "\nTable of Expected-s[m] and P[s,m]\n";
printTable($Psm, [1..20], [map(10*$_,(1..22))]);

print "\nTable of P[s,220] / ($root)^s\n";
foreach my $s (1..20) { printf " %6d", $s; }
print "\n";
my $power = 1;
foreach my $s (1..20) {
    $power *= $root;
    printf " %6.4f", $Psm->value($s,220)/$power;
}
print "\n";


sub averagePairScore {     ## Compute average score per pair.
    my ($p,$s) = @_;
    my $avg;
    foreach my $i (0..$#$p) {
	foreach my $j (0..$#$p) {
	    $avg += $$p[$i]*$$p[$j]*$$s[$i][$j];
	}
    }
    return $avg;
}    

sub readScoreScheme {
    ## Read background frequencies into @p.
    print "Background Frequencies (input)\n";
    $_ = <STDIN>;
    my @p = split;
    print;
    
    ## Read scoring matrix into @s; remember largest entry.
    print "\nScoring Matrix (input)\n";
    my @s;
    while (<STDIN>) {
	print;
	push @s, [split];
    }
    return (\@p, \@s);
}

sub printTable {
    my ($rec, $sbounds, $mbounds) = @_;
    print "\ns\\m";
    foreach my $m (@$mbounds) { 
	printf " %6d", $m;
    }
    print " ratio\n E:";
    foreach my $m (@$mbounds) { 
	my $sum=0;
	for (my $s=1; $rec->value($s,$m)>0; $s++) {
	    $sum += $rec->value($s,$m);
	}
	printf " %6.4f", $sum;
    }
    foreach my $s (@$sbounds) {
	printf "\n%2d:", $s;
	foreach my $m (@$mbounds) { 
	    printf " %6.4f", $rec->value($s,$m);
	}
	my $m = $$mbounds[-1];
	printf " %6.4f", 
            $rec->value($s,$m) / $rec->value($s-1,$m)
	    if $s>1;
    }
    print "\n";
}
