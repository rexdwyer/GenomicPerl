#!/usr/bin/perl -I /home/dwyer/book/perl

use strict;

my @residues = split //,"acdefghiklmnpqrstvwy";   ## list of 20 amino acids

my %lod;

sub fillLod {
    my ($matrixFile) = @_;
    open LOD, $matrixFile or die "can't open";
    my ($trash,@residues) = split /\s+/, <LOD>;
    while (<LOD>) {
	my ($r,@scores) = split;
	foreach (0..19) {
	    $lod{$r.$residues[$_]} = $scores[$_];
	}
    }
}
    


sub countSimilarWords {
    my ($word) = @_;
    my ($w1,$w2,$w3) = split //,$word;
    return 0 if $lod{$w1.$w1}+$lod{$w2.$w2}+$lod{$w3.$w3}<11;
    my $similar;
    foreach my $r1 (@residues) {
	foreach my $r2 (@residues) {
	    my $t = 11-$lod{$w1.$r1}-$lod{$w2.$r2};
	    foreach my $r3 (@residues) {
		$similar++ if $lod{$w3.$r3}>=$t;
	    }
	}
    }
    return $similar;
}
       

fillLod("blosum62");
my $total=0;
my $min=8000;
my $max=0;
my @counts;
my ($maxs,$mins);
foreach my $r1 (@residues) {
    foreach my $r2 (@residues) {
	foreach my $r3 (@residues) {
	    my $cnt = countSimilarWords($r1.$r2.$r3);
	    push @counts, $cnt;
	    $total += $cnt;
	    if ($cnt > $max) { $max=$cnt; $maxs = $r1.$r2.$r3; }
	    if ($cnt < $min) { $min=$cnt; $mins = $r1.$r2.$r3; }
}}}
print "$mins $min  $maxs $max\n";
print $total/8000, " $total\n";
@counts = sort {$a <=>$b} @counts;
print "$counts[4000] $counts\[3999]\n";


