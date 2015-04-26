#!/usr/bin/perl

use strict; 

my @structure;
my @bases;

sub max {
    my ($x,$y) = @_;
    if ($x>$y) { return $x; } else { return $y};
}

sub evalRna {
    my ($l, $r) = @_;
    my %bonds = (GU=>1,UG=>1,AU=>2,UA=>2,CG=>3,GC=>3);
    my $energy = $bonds{$bases[$l].$bases[$r]};
    my $level=0;
    
    my $ii = $l;
    
    for (my $i=$l+1; $i<=$r; $i++) {
	$level-- if ($structure[$i] eq ")");
	if ($level==0) {
	    $energy += evalRna($ii,$i) if ($structure[$i] eq ")");
	    $ii = $i;
	}
	$level++ if ($structure[$i] eq "(");
    }
    return $energy;
}

sub evalRnaStructure {
    my ($basestring,$structurestring) = @_;
    @bases = split(//, 5 . $basestring . 3);
    @structure = split(//, "($structurestring)");
    return evalRna(0, $#structure);
}

my $basestring = <STDIN>;
chomp($basestring);
my $parenstring = <STDIN>;
chomp($parenstring);
print evalRnaStructure($basestring,$parenstring), 
     " hydrogen bonds in this structure.\n";
