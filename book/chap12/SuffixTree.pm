package SuffixTree;
use strict;

my $debugflag;
sub debug { print(@_) if $debugflag;}
sub debugf { printf(@_) if $debugflag;}

#########################################
sub new
## 
#########################################
{
    my ($this, $target) = @_;
    my $root = {off=>0,len=>0};
    $$root{sfx} = $root;
    $target = lc $target;
#   my $cpy = $target;
#   my $cnt = 0;
#   while ($cpy =~ s/(.{1,50})//) { print $1, " ", ($cnt+=50), "\n"; }
#   print length $target, " ", substr($target,-10), "\n";
    my ($i,$j,$p,$pDepth,$needSfx) = (0,0,$root,0,0);
    while ($j < length $target) {
#	print  "($i,$j) --> ";
	($i,$j,$p,$pDepth,$needSfx) = add($target,$i,$j,$p,$pDepth,$needSfx);
	my $ci = substr($target,$i,1);
	my $cj = substr($target,$j,1);
#	print  "($i,$j) $ci, $cj\n";
	debug "\n\nRoot Dump:\n";
	dumpNode($root, $target, "") if $debugflag;
	$pDepth=0 if $pDepth<0;
	$needSfx = undef if $needSfx == $root;
	$j++ if ($i>$j);
    }
    $|=1;
    $this = bless {target=>$target, root=>$root}, (ref($this)||$this);
#   $this->dump();
    return $this;
}


sub add {
    my ($target,   ## the target string
	$i,        ## the position in $target at which the substring starts
	$j,        ## the position in $target at which the substring ends
	$p,        ## the suffix-tree node from which searching is to begin 
	$pDepth,   ## number of letters in path from root through $p's parent
	$needSfx)  ## a new internal node waiting for a suffix pointer, if any
	= @_;
    
    shift @_;
    debugf "\n\nadd(%s) %s\n", join(',',@_), substr($target,$i,$j-$i+1);
    debug  "add($i,$j)\n";
    my $sLeft = ($j+1-$i)-$pDepth;  ## length of unmatched portion
    my $sLetter = substr($target,$j,1); ## new letter to append
    my ($q,$r,$pLen);
    while ($p && ($sLeft>($pLen=$$p{len}))) {  ## Move down as far as possible.
	($pDepth,$sLeft) = ($pDepth+$pLen,$sLeft-$pLen);
	($r,$q,$p) = ($q,$p,$p->{substr($target,$i+$pDepth,1)});
    }
    my $qDepth = $pDepth - $$q{len};
    if (!$p) { ## Substring is absent; we fell off the tree.
	debug "Case i: Add Tail to $q\n";
	$$q{$sLetter} = {off=>$j,len=>(length($target) - $j),startsAt=>$i};
	debug "$q now points to $$q{$sLetter}\n";
	$$needSfx{sfx} = $q if ref($needSfx);
	if ($$q{sfx}) {  ## parent already internal node.
	    debug "Case i.a: $q was already internal (or root)\n";
	    return ($i+1, $j, $$q{sfx}, $qDepth-1, undef);
	} else {   ## parent becomes internal node.
	    debug "Case i.b: $q has become internal\n";
	    return ($i+1, $j, $$r{sfx}, $qDepth-$$r{len}-1, $q);
	}
    } else { ## We're still on the tree, but the last letter may mismatch.
	my ($pOff, $pLen) = ($$p{off},$$p{len});
	my $pLetter = substr($target,$pOff+$sLeft-1,1);
	if ($sLetter eq $pLetter) { ## Suffix is already there
	    debug "Case ii: Found at node $p\n";
	    debug "\n\nNode Dump:\n";
	    dumpNode($p, $target, "") if $debugflag;
	    $$needSfx{sfx} = $q if ref($needSfx);
	    return ($i, $j+1, $q, $qDepth, undef);
	} else { ## Must split node to add suffix
	    debug "Case iii: Splitting node $p\n";
            ## One new node contains the rest of this suffix.
	    my $new1 = {off=>$j,len=>(length($target)-$j),startsAt=>$i};  

            ## The old node retains the latter part of its original key.
	    $$p{off} = $pOff+$sLeft-1;
	    $$p{len} = $pLen-$sLeft+1;

            ## A second new node contains common start of key and substring.
            ## It points to the two differing endings
	    my $new2 = {off=>$pOff, len=>$sLeft-1};
	    $$new2{$pLetter} = $p;
	    $$new2{$sLetter} = $new1;

            ## Parent of original node $p must now point to new one.
	    $q->{substr($target,$i+$pDepth,1)} = $new2;
	    $$needSfx{sfx} = $new2 if ref($needSfx);
#	    debug "$new2 now points to $p and $new1\n";
            ## We've put whole suffix in a new node; we're done.
	    return ($i+1, $j, $$q{sfx}, $qDepth-1, $new2); 
	}
    }
}

#########################################
sub isSubstring
## 
#########################################
{
    my ($this,$query) = @_;
    my $s = $this->{target};
    $query = lc $query;
    my $p = $this->{root};
    my $qLen = length $query;
    while ($query) {
	$p = $$p{substr($query,0,1)};
	return undef if !$p;
	my $patLen = min($$p{len},$qLen);
	my $pStr = substr($s, $$p{off}, $patLen);
	return undef if $query !~ s/^$pStr//;
	$qLen -= $patLen;
    }
    return $p;
}

#########################################
sub gatherStartPoints
## 
#########################################
{
    my ($this, $p) = @_;
    my @results;
    my @pending = ($p);
    while (@pending) {
	$p = shift @pending;
	push @results, $$p{startsAt} if defined($$p{startsAt});
	push @pending, map($$p{$_}, grep(/^.$/, keys %$p));
    }
    return @results;
}

#########################################
sub dump
## 
#########################################
{
    my ($this) = @_;
    print "Suffix Tree for $$this{target}:\n";
    dumpNode($this->{root}, $this->{target}, " ");
    print "\n\n";
}

sub dumpNode {
    my ($p,$s,$indent) = @_;
    my $ss = substr($s, $$p{off}, $$p{len});
    if (length $ss > 19) { $ss = substr($ss,0,13). "..." . substr($ss,-3); }
    print "$indent$ss", (" "x(20 - (length $indent) - (length $ss)));
    print "$p ($$p{off},$$p{len})";
    print "  -->$$p{sfx}" if $$p{sfx};
    print "\n";
    foreach (sort keys %$p) { 
	dumpNode($$p{$_}, $s, $indent.('"'x length $ss)) if 1==length $_;
    }
}

1;
