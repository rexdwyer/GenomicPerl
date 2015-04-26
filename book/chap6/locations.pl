#!/usr/bin/perl -I . -I /home/dwyer/book/perl

use strict;

sub complement {
    my ($s) = @_;
    $s =~ tr/acgtACGT/tgcaTGCA/;
    return reverse $s;
}

sub ident {
    my ($s) = @_;
    return $s;
}

sub extract_location {
    my ($loc, $s) = @_;
    print @_, "\n";
    ## We don't want to create lots of copies of a long string.
    unless (ref $s) {my $t = $s; $t =~ s/\s//g; $s = \$t; }
    $loc = lc $loc;
    $loc =~ s/,/\./g;
    print "$loc\n";
    $loc =~ s/(\d+)/N$1/g;
    print "$loc\n";
    $loc =~ s/join/ident/g;
    print "$loc\n";
    $loc =~ s/N(\d+)\.\.N(\d+)/substr(\$\$s,$1-1,$2+1-$1)/g;
    print "$loc\n";
    $loc =~ s/N(\d+)/substr(\$\$s,$1-1,1)/g;
    print "$loc\n";
    print eval($loc),"\n";
    print $@, "\n";
    return eval($loc);
}



my $seq = "ourfatherwhoartinheavenhallowedbethyname";
my $location = "join(18..24,complement(join(5..10,6..11)),27,40..45)";

extract_location($location,$seq);

