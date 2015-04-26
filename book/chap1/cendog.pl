#!/usr/bin/perl -w 
use strict;   ## request conservative error checking

my %codonMap;  ## declare a hash table

## transcribe translates DNA strings to RNA
sub transcribe {
    my ($dna) = @_;
    my $rna = scalar reverse $dna;
    $rna =~ tr/ACGT/UGCA/;
    return $rna;
}


## translate translates mRNA strings to proteins
sub translate {
    my ($mrna) = @_;
    my $pro = "";
    while ( $mrna =~ s/(...)// ) {
        $pro = $pro . $codonMap{$1};
    }
    return $pro;
}


## Construct hash that maps codons to amino acids by reading table
## from DATA at the end of the program, which have the form:
## Residue Codon1 Codon2 ...

while (my $in = <DATA> ) {  ## assigns next line of DATA to $in; fails if none
    chomp($in);             ## remove line feed at end of $in
    my @codons = split " ",$in;
    my $residue = shift @codons;  ##remove first item from @codons and assign

    foreach my $nnn (@codons) {
        $codonMap{$nnn} = $residue;
    }
}

## Now read DNA strands from input <STDIN> and print translations in all six
## possible reading frames

while ( my $dna = <STDIN> ) {
    chomp($dna);
    print "DNA: ", $dna, "\n";
    my $rna = transcribe($dna);
    print "RNA: ", $rna, "\n";
    my $protein = translate($rna);
    print "RF1: ", $protein, "\n";
    $rna =~ s/.//;
    $protein = translate($rna);
    print "RF2:  ", $protein, "\n";
    $rna =~ s/.//;
    $protein = translate($rna);
    print "RF3:   ", $protein, "\n\n";
}
    
## The lines below are not perl statements and are not executed as part of the 
## program.  Instead, they are available to be read as input by the program
## using the file handle DATA.
__END__
Ala GCU GCC GCA GCG
Arg CGU CGC CGA CGG AGA AGG
Asn AAU AAC
Asp GAU GAC 
Cys UGU UGC
Gln CAA CAG
Glu GAA GAG
Gly GGU GGC GGA GGG
His CAU CAC
Ile AUU AUC AUA
Leu UUA UUG CUU CUC CUA CUG
Lys AAA AAG
Met AUG
Phe UUU UUC
Pro CCU CCC CCA CCG
Ser UCU UCC UCA UCG AGU AGC
Thr ACU ACC ACA ACG
Trp UGG
Tyr UAU UAC
Val GUU GUC GUA GUG
... UAA UAG UGA
