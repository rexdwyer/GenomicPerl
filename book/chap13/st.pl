#!/usr/bin/perl -I . -I /home/dwyer/book/perl

use strict;
use SuffixTree;

$|=1;

my $t = SuffixTree->new("banana");
$t->dump();

my $t = 
bless
{target=>'banana$',
 root=>
   {off=>0, len=>0,
    'a'=>{off=>1, len=>1,    
          'n'=>{off=>2, len=>2,
                'n'=>{off=>4, len=>3},
                '$'=>{off=>6, len=>1}
               },
          '$'=>{off=>6, len=>1}
         },
    'b'=>{off=>0, len=>7},
    'n'=>{off=>2, len=>2,
          'n'=>{off=>4, len=>3},
          '$'=>{off=>6, len=>1}
         },
    '$'=>{off=>6, len=>1}
   }
},
"SuffixTree";
$t->dump();



my $t = SuffixTree->new("abracadabra");
$t->dump();
my $t = SuffixTree->new("bookkeeper");
$t->dump();
print $t->lookUp("bookkeeper");
print $t->lookUp("ookkeeper");
print $t->lookUp("okkeeper");
print $t->lookUp("kkeeper");
print $t->lookUp("keeper");
print $t->lookUp("eeper");
print $t->lookUp("eper");
print $t->lookUp("per");
print $t->lookUp("er");
print $t->lookUp("r");
print $t->lookUp("");
print "\n";
print $t->lookUp("beekeeper");
print $t->lookUp("eekeeper");
print $t->lookUp("ekeeper");
print $t->lookUp("beekiller");
print $t->lookUp("eekiller");
print $t->lookUp("ekiller");
print $t->lookUp("killer");
print $t->lookUp("iller");
print $t->lookUp("ller");
print $t->lookUp("ler");
print "\n";

my $t = SuffixTree->new("baaaaaaaaaaaaaaa");
my $t = SuffixTree->new("misssisssisssisssippippi");


while (<STDIN>) {
    chomp;
    my $t = SuffixTree->new($_);
#    $t->dump();
     print "$_\n";
}






