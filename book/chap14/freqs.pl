#!/usr/bin/perl -I . -I /home/dwyer/book/perl

use strict;
use Util;
use CSPredictor;


print "10\n";
foreach my $b1 ('a','c','g','t') {
foreach my $b2 ('a','c','g','t') {
foreach my $b3 ('a','c','g','t') {
foreach my $b4 ('a','c','g','t') {
foreach my $b5 ('a','c','g','t') {
foreach my $b6 ('a','c','g','t') {
    print "$b1$b2$b3$b4$b5$b6 ", (rand()/1500), "\n";
}}}}}}    

exit(0);



