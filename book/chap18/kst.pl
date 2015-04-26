#!/usr/bin/perl -I /home/dwyer/book/perl

## An implementation of the Kaplan-Shamir-Tarjan algorithm
## for sorting by reversals

use strict;
use Util;
use UnionFind;

my  @pi;

### Left and right endpoints of a wish.
sub L { $_[0][0]; }
sub R { $_[0][1]; }

### True/False if wish is twisted/simple
sub twisted { odd(L($_[0])+R($_[0])) }     ## = KST's unoriented
sub simple { !odd(L($_[0])+R($_[0])) }     ## = KST's oriented
sub odd { (1 & shift) }

### Prints a permutation with brackets indicating a reversal.
### Arguments $left, $right indicate limits of reversal.
### Third argument is a reference to a list containing the
### unsigned permutation.
sub printBracketedPerm {
    my ($left, $right, $perm) = @_;
    my $width = 3 + int(log(max(@$perm)/2)/log(10.0));
    print "\n";
    for	(my $i=1; $i<$#$perm; $i+=2) {
	my $number =  odd($$perm[$i]) ?  ($$perm[$i]+1)/2 :  -($$perm[$i]/2);
	if ($i==$left) {
	    print substr( "     [$number", -$width );
	} elsif ($i==$right+1) {
	    print "]",substr( "      $number", 1-$width );
	} else {
	    print substr( "      $number", -$width );
	}
    }
    print (($right==$#$perm-1) ? "]" : " ");
}

### Prints a permutation with no brackets indicating a reversal. 
### Single argument is a reference to a list containing the
### unsigned permutation.
sub printPerm {
    printBracketedPerm(-1,-1,$_[0]);
}

### Prints a list of pairs, such as a wish list.
sub printPairList {
    my ($label, @edges) = @_;
    print "$label: ";
    foreach my $e (@edges) {
	print "[$$e[0],$$e[1]] ";
    }
    print "\n";
}

### Draws a wish list as a series of bars representing
### intervals on the number line.
sub drawWishList {
    my (@edges) = @_;
    foreach my $e (@edges) {
	my ($lo,$hi) = @$e;
	my $c = twisted($e) ? "T" : "S";
	print(("  "x$lo), $c, ("--"x($hi-$lo)), "$c\n");
    }
}

### Returns a random signed permutation of 1..n
sub randomSigned {
    my ($n) = @_;
    my @spi = (1..$n);
    foreach my $i (1..$n) {
##	$spi[$i-1] *= -1 if (rand() < 0.5);
	my $j = int(rand($i));
	my $t = $spi[$i-1];
	$spi[$i-1] = $spi[$j];
	$spi[$j] = $t;
    }
    return @spi;
}

### Converts a signed permutation of 1..n into
### the corresponding unsigned permutation of 1..2n
### for sorting
sub signedToUnsigned {
    push my @pi,0;
    foreach my $i (0..$#_) {
	my $j = $_[$i];
	if ($j > 0) { 
	    push @pi, 2*$j-1, 2*$j;
	} else {
	    push @pi, 2*(-$j), 2*(-$j)-1;
	}			
    }
    push @pi,2*$#_+3;
    return @pi;
}

### Converts an unsigned permutation of 1..2n into
### the corresponding signed permutation of 1..n for display
sub unsignedToSigned {
    my @spi;
    for (my $i=1; $i<$#_; $i+=2) {
	if (odd($_[$i])) {
	    push @spi, ($_[$i]+1)/2;
	} else {
	    push @spi, -($_[$i]/2);
	}			
    }
    return @spi;
}

### Creates the "wish list" for a given unsigned permutation.
### Argument is a reference to a list containing the permutation.
sub makeWishList {
    my ($pi) = @_;
    my @piInv = (1) x @$pi;
    foreach my $i (0..$#pi) { $piInv[$$pi[$i]] = $i; }
    my @wishes = ();
    for (my $i=0; $i<$#piInv; $i+=2) {
	my ($j,$j1) = ($piInv[$i], $piInv[$i+1]);
	($j,$j1) = ($j1,$j) if ($j1<$j);
	next if ($j == $j1-1);
	push @wishes, [$j, $j1];
    }
    return(sort { L($a) <=> L($b) } @wishes);
}

### Finds a "happy clique" in a given wish list.
### Argument is the wish list (not a reference).
### Return value is the list of wishes in the happy clique.
sub findHappyClique {
    my @wishes= @_;
    my (@C, $t);
    foreach my $w (@wishes) {
	if (twisted($w)) { next; }         ## wish not simple
	elsif (@C==0) { @C = ($w); next;}   ## first simple wish
	elsif (R(@C[$#C]) < L($w)) { last; }  ## doesn't overlap clique; done
	elsif ($t && R($w) > R($t)) { next; } ## overlaps t; case 1
	elsif (R($w) < R($C[$#C])) { ## case 2c
	    $t = pop @C;  ## this one contains new one
	    @C = ($w);     ## start over; clique contains only new one.
	} elsif (L($w) <= R($C[0])) { ## case 2a; add
	    push @C, $w;
	} elsif (L($w) > R($C[0])) { ## case 2b; start over
	    @C = ($w);
	}
    }
    return @C;
}

### Given a happy clique and a wish list (both by reference),
### findMaxTwistedDegree returns the wish that overlaps the
### largest number of twisted wishes.
sub findMaxTwistedDegree {
    my ($clique, $wishes) = @_;
    my @bounds = ((map { L($_); } @$clique), (map { R($_); } @$clique));
#   print join(",", @bounds), " bounds\n";
    my $j = @$clique;
    my @o = (0) x ($j + 2);
    foreach my $w (@$wishes) {
	if (simple($w)) { next; }
	my ($l,$r) = (0,0);
	foreach (@bounds) {
	    $l++ if (L($w)>$_);
	    $r++ if (R($w)>$_);
	}
	if ($l == @bounds || $r == 0) { next;}
	if ($r == @bounds && $l == 0) { next;}
	if ($r == $l) { next;}
	if ($r<=$j) {
	    $o[$l+1]++;
	    $o[$r+1]--;
	} elsif ($l>=$j) {
	    $o[$l-$j+1]++;
	    $o[$r-$j+1]--;
	} else {
	    my $m = min($l, $r-$j);
	    if ($m>0) { $o[1]++; $o[$m+1]--; }
	    my $M = max($l, $r-$j);
	    if ($M<$j) { $o[$l+1]++; }
	}
#       print join(',',@$w), ", l=$l, r=$r\n";
#	print join(",", @o), "\n";
    }
#   print join(",", @o), "\n";
    my $sum=0;
    my $max=0;
    foreach (0..$#o) {
	$sum += $o[$_]; 
	$o[$_] = $sum;
	$max = $_ if ($sum > $o[$max]);
    }
    return $$clique[$max-1]; ## adjust down 1 for Perl array indexing
}

### Carries out an even permutation (twist) of a permutation,
### and prints the result.
### Arguments: $lo, $hi: indices of limits of twist.
###     $piref: reference to list containing permutation.
###     $note: annotation to add to printed version.
sub twist {
    my ($lo, $hi, $piref, $note) = @_;
    $lo++ if !odd($lo);
    $hi-- if odd($hi);
    printBracketedPerm($lo,$hi,$piref);
    for (0  ;$lo < $hi; $lo++, $hi--) {
	($$piref[$lo],$$piref[$hi]) = ($$piref[$hi],$$piref[$lo]);
    }
    print $note;
}

### Carries out an even permutation (twist) on the elements
### of a wish list.  Not used by the program to do its job,
### but sometimes helpful for debugging or tracing.
### Arguments: $lo, $hi: indices of limits of twist.
###     $wishes: reference to wish list
### Side Effect: the wish list is changed.
sub wishListTwist {
    my ($lo, $hi, $wishes) = @_;
    $lo++ if !odd($lo);
    $hi-- if odd($hi);
    foreach my $w (@$wishes) {
	$$w[0] = $hi - (L($w) - $lo) if ($lo <= L($w) && L($w) <= $hi);
	$$w[1] = $hi - (R($w) - $lo) if ($lo <= R($w) && R($w) <= $hi);
	($$w[0],$$w[1]) = ($$w[1],$$w[0]) if L($w)>R($w);
    }
}

### Finds the block (equivalence classes) of the transitive
### closure of the "overlaps" relation.
### Argument: reference to a wish list.
### Return value: reference to a union/find object.
sub findOverlapBlocks {
    my ($wishes) = @_;
    my $ov = new UnionFind();
    my @handle;
    my @events;
    foreach (@$wishes) {
	$ov->union(L($_),R($_));
	$handle[L($_)] = $handle[R($_)] = $_;
	push @events, ([L($_), "add"], [R($_), "remove"]);
    }
    @events = sort { $$a[0] <=> $$b[0] } @events;

#   printPairList("findOverlapBlocks: about to scan", @events);
    
    my @stack= (-1);
    foreach my $e (@events) {
#       print "stack: ", join(",",@stack), "\n";
#	print "e = [", join(",",@$e),  "]\n";
	my ($endpt, $action) = @$e;
	if ($action eq "add") { # entering interval of wish
	    push @stack, $endpt; # make interval active
	} elsif ($action eq "remove") { # leaving interval
	    my $block = $ov->find($endpt);
	    my $Lh = L($handle[$block]);
#	    print "block = $block; Lh = $Lh\n";
	    while ($Lh <= $stack[$#stack]) {
		my $block1 = $ov->find(pop @stack);
		$block = $ov->union($block, $block1);
		$handle[$block] = $handle[$block1] 
                    if R($handle[$block1]) > R($handle[$block]);
	    }
	    push @stack, L($handle[$block]) if R($handle[$block]) > $endpt;
	}
    }
    return $ov;
}

### Classifies the blocks of the overlap relation as "simple",
### "twisted", "simple hurdle", or "superhurdle".
### Arguments: references to a wish list and to a union-find object.
### Returns: counts of simple and super- hurdles, list of contiguous
### runs, and a list of classifications for the blocks.

sub classifyBlocks {
    my ($wishes, $ov) = @_;

    my @blockType;
    foreach my $e (@$wishes) {
	my $block = $ov->find(L($e));
	$blockType[$block] = (simple($e) && "simple") 
	    || $blockType[$block] || "twisted";
    }
    my @runs = (-1);   ## condense twisted runs to a single entry
    my @numRuns;       ## count number of runs for each block
    
    foreach (map { $ov->find($_) } (0..$#pi)) {
	push @runs, $_ 
	    if ($blockType[$_] eq "twisted") && ($runs[$#runs] != $_);
    }
    shift @runs;   ## remove the -1

    if (@runs==1) { ## special case for one giant twisted block...
	$blockType[$runs[0]] = "simple hurdle";
	return (0,1,\@runs,\@blockType);
    }
    
    foreach (@runs) { $numRuns[$_]++; }
    
    if ($runs[0] == $runs[$#runs]) { $numRuns[0]--; }
    else { push @runs, $runs[0] }
    push @runs, $runs[1];
    
    my ($superhurdles, $simplehurdles);
    for (my $i=1; $i<$#runs; $i++) {
	if ($numRuns[$runs[$i]] == 1) {    ## consecutive...a hurdle
	    if (($runs[$i-1] == $runs[$i+1]) && ($numRuns[$runs[$i-1]] == 2)) {
		$blockType[$runs[$i]] = "superhurdle";
		$superhurdles++;
	    } else {
		$blockType[$runs[$i]] = "simple hurdle";
		$simplehurdles++;
	    }
	}
    }
    $#runs -= 2;
    return ($superhurdles, $simplehurdles, \@runs, \@blockType);
}


### Searches the list of runs for blocks of given types to be merged.
### Accepts up to three types, to be found in succession.
### Arguments:  $runs, ref. to list of runs of blocks.
###     $blockType: ref. to list of types of the blocks
###     $type1,$type2,$type3: block types to search for.

sub blockSearch {
    my ($runs, $blockType, $type1, $type2, $type3) = @_;
    my ($h1,$h2,$h3);
    for ($h1=0; $h1<@$runs; $h1++) {
	last if ($$blockType[$$runs[$h1]] eq $type1);
    }
    for ($h2=$h1+1; $h2<@$runs; $h2++) {
	last if ($$blockType[$$runs[$h2]] eq $type2);
    }
    for ($h3=$h2+1; $h3<@$runs; $h3++) {
	last if ($$blockType[$$runs[$h3]] eq $type3);
    }
    return ($type1,$$runs[$h1]) if $h2 >= @$runs;   ## no h2
    return ("$type1+$type2", $$runs[$h1],$$runs[$h2]) if $h3 >= @$runs; ## no h3
    return ("$type1+$type3", $$runs[$h1],$$runs[$h3]);
}

### Clears at least one hurdle from the current permutation.
### Arguments:  $wishes, ref to wish list.
###     $ov, ref to union/find object describing blocks.
### Returns: 1 if successful; 0 if no more hurdles to clear.
sub clearHurdles {
    my ($wishes, $ov) = @_;
    my ($superhurdles, $simplehurdles, $runs, $blockType)
	= classifyBlocks($wishes,$ov);

#   printPairList("\nwishes", @$wishes);
#   printPairList("\nsets", map { [$_, $ov->find($_)] } (0..41));
#   printPairList("\nblocks", map { [$_, $$blockType[$_]] } @$runs);
    
    my ($note,$h1,$h2);
    if ($superhurdles+$simplehurdles == 0) { return 0; }
    if ($superhurdles==0 && $simplehurdles == 1) {
	($note,$h1) = blockSearch($runs,$blockType,"simple hurdle");
	### untwist wish with leftmost endpoint.
	my ($i,$j);
	for ($i=0; $i<@pi; $i++) { last if $ov->find($i) == $h1};
	foreach (@$wishes) { $j = R($_) if (L($_)==$i) }
	twist($i, $j,\@pi, " >twist simple hurdle>");
	return 1;
    }

    ### Choose two blocks to twist together.
    if ($superhurdles>=2) {
	($note,$h1,$h2) 
	    = blockSearch($runs,$blockType,"superhurdle","superhurdle");
    } elsif (($superhurdles==1) && ($simplehurdles>=1)) {
	($note,$h1,$h2) 
	    = blockSearch($runs,$blockType,"superhurdle","simple hurdle");
	($note,$h1,$h2) 
	    = blockSearch($runs,$blockType,"simple hurdle","superhurdle")
		unless $h2;
    } elsif ($superhurdles==1) {
	($note,$h1,$h2) 
	    = blockSearch($runs,$blockType,"superhurdle","twisted");
	($note,$h1,$h2) 
	    = blockSearch($runs,$blockType,"twisted","superhurdle")
		unless $h2;
    } elsif ($simplehurdles>=2) {
	($note,$h1,$h2) 
	    = blockSearch($runs,$blockType,
			  "simple hurdle","simple hurdle","simple hurdle");
    }

    print ((map { "[$_, $$blockType[$_]] " } ($h1,$h2)), "\n");
    
    ### Find boundaries of blocks.
    my $end1, my $start2;
    for (my $i=0; $i<@pi; $i++) { $end1=$i if $ov->find($i) == $h1};
    for (my $i=@pi-1; $i>=0; $i--) { $start2=$i if $ov->find($i) == $h2};
    if ($h1==$$runs[0]) {   ## First hurdle could wrap around
	for ($end1=$start2; $end1>=0; $end1--) {
	    last if $ov->find($end1) == $h1;
	}
    }
    if ($h2==$$runs[0]) {   ## Second hurdle could wrap around
	for ($start2=$end1; $start2<@pi; $start2++) {
	    last if $ov->find($start2) == $h2;
	}
    }
    ### Twist together.
    twist($end1, $start2,\@pi, " >twist $note>");
    return 1;
}

##MAIN
{
    srand(13);
    foreach (1..100) {
	my @spi = randomSigned(20);
#	@spi = (1,13,11,2,14,8,3,5,4,6,10,9,7,15,12,16,19,17,18,20);
	@pi = signedToUnsigned(@spi);
	printPerm(\@pi);
	
	my $count;
	while (1) {
	    my @wishes = makeWishList(\@pi);
	    my $ov = findOverlapBlocks(\@wishes);
	    last unless clearHurdles(\@wishes, $ov);
	    $count++;
	}
	
	while (1) {
	    my @wishes = makeWishList(\@pi);
	    my @happy = findHappyClique(@wishes);
#           printPairList("\nwishes",@wishes);
#           drawWishList(@wishes);
#	    printPairList("\nhappy",@happy);
	    last unless @happy;
	    $count++;
	    my @interval = @{findMaxTwistedDegree(\@happy, \@wishes)};
	    twist(@interval,\@pi," >grant wish>");
	}
	printPerm(\@pi);
	print ":    $count reversals required\n";
	my $i = 0;
	foreach(@pi) { print join(",",@spi), "   ERROR!!\n" if ($_ != $i++); }
#	print "\n\n\n";
    }
}
