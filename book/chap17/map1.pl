#!/usr/bin/perl -I /home/dwyer/book/perl

## An implementation of the C1P algorithm described in Section 5.3

use strict;

## define clone object, access functions
sub newClone {
    return [$_[0], $_[1], undef, undef, undef];
}

sub id { $_[0][0]; }
sub probes { $_[0][1]; }
sub cComp { $_[0][2]; }
sub L { $_[0][3]; }
sub R { $_[0][4]; }

sub setcComp { $_[0][2] = $_[1]; }
sub setL { $_[0][3] = $_[1]; }
sub setR { $_[0][4] = $_[1]; }

## define (connected) component (of Gc) object, access functions
sub newComponent {
    return [shift,undef,[],[],1];
}

#sub id { $_[0][0]; }        # also works for components
#sub probes { $_[0][1]; }    # also works for components
sub clones { $_[0][2]; }
sub order  { $_[0][3]; }
sub offset { $_[0][4]; }

sub setId { $_[0][0] = $_[1]; }
sub setProbes { $_[0][1] = $_[1]; }
sub setClones { $_[0][2] = $_[1]; }
sub setOrder  { $_[0][3] = $_[1]; }
sub setOffset { $_[0][4] = $_[1]; }

sub printComponent {
    my $c = shift;
    printf "   C%d: %s [%s] {%s} off=%d\n", id($c), probes($c),
    join(",",@{order($c)}), join(",",@{clones($c)}),offset($c);
}

sub printComponentList {
    my ($label, $comps) = @_;
    print "$label:\n";
    map {printComponent($_)} @$comps;
    print "\n";
}

################################################################
########  ADDED CODE BEGINS HERE
sub forAndBackwards {
    return map { ($_, scalar reverse $_) } @_ ;
}

sub cartesianProd {
    my @lists = @_;
    return [] if (@lists == 0);
    return @{$lists[0]} if (@lists == 1);
    my $list1 = $lists[0];
    my @list2 = cartesianProd(@lists[1..$#lists]);
    my @ans =
	map { 
	    my $e2 = $_;
	    map { $_ . $e2 } @$list1
	    } @list2;
    return @ans;
}


sub permute {
    my @lists = @_;
    return [] if (@lists == 0);
    return $lists[0] if (@lists == 1);
    my @ans = ();
    for (my $i=0; $i<@lists; $i++) {
	my $choice = $lists[$i];
	my @shortlist = (@lists[0..$i-1],@lists[$i+1..$#lists]);
	@ans = (@ans, cartesianProd($choice, permute(@shortlist)));
    }
    return \@ans;
}


sub allOrders {
    my ($order) = @_;
    return [$order] if (!ref($order));
    if ($$order[0] eq "*SET*") {
	my @orders = map { allOrders($_) } @$order[1..$#$order];
	my $ans = permute(@orders);
	return $ans;
    } else {
	my @orders = map { allOrders($_) } @$order;
	my @ans = forAndBackwards(cartesianProd(@orders));
	return \@ans;
    }
}
########  ADDED CODE ENDS HERE
################################################################

sub orderString {
    my ($order) = @_;
    return $order if (!ref($order));
    if ($$order[0] eq "*SET*") {
	return "{" 
	    . join(",", map { orderString($_) } @$order[1..$#$order])
		. "}";
    } else {
	return "["
	    . join("", map { orderString($_) } @$order)
		. "]";
    }
}

sub printComponentOrder {
    my ($label, $cc) = @_;
    print "$label:\n", orderString(order($cc)), "\n";
}

sub addToComponent {
    my ($clone, $cc) = @_;
    setcComp($clone, id($cc));
    push @{clones($cc)}, id($clone);
    setProbes($cc, setUnion(probes($cc),probes($clone)));
}

sub setUnion {
    my ($s1,$s2) = @_;
    return $s1 unless $s2;
    $s1 =~ s/[$s2]//g;
    return $s1 . $s2;
}


sub setIntersect {
    my ($s1,$s2) = @_;
    $s1 =~ s/[^$s2]//g;
    return $s1;
}

sub setSubset {
    my ($s1,$s2) = @_;
    $s1 =~ s/[$s2]//g;
    return !$s1;
}

sub setOverlap {
    my ($s1,$s2) = @_;
    ($s1,$s2) = ($s2,$s1) if length($s1) > length($s2);
    my $sc = $s1;
    $s1 =~ s/[^$s2]//g;
    return "" if $s1 eq $sc;    ## c1 subset of c2
    return $s1;
}	

sub cloneOverlap {
    setOverlap(map { probes($_) } @_);
}

sub intersectSize {
    length setIntersect(map { probes($_) } @_);
}

sub cloneSize {
    length (probes(shift));
}

sub printClone {
    my ($c) = @_;
    my ($id,$cc,$l,$r,$p) = (id($c), cComp($c), L($c), R($c), probes($c));
    print ("  [\#$id: $cc($l,$r), $p]\n");
}

sub printCloneList {
    my ($label, $clones) = @_;
    print "$label:\n";
    map {printClone($_)} @$clones;
    print "\n";
}


sub placeLeft { ## place u to left of v
    my ($u,$v) = @_;
    print "placing u=$$u[0] to left of v=$$v[0]\n";
    setR($u, L($v) + intersectSize($u,$v) - 1);
    setL($u, R($u) - cloneSize($u) + 1);
}

sub placeRight { ## place u to right of v
    my ($u,$v) = @_;
    print "placing u=$$u[0] to right of v=$$v[0]\n";
    setL($u, R($v) - intersectSize($u,$v) + 1);
    setR($u, L($u) + cloneSize($u) - 1);
}

sub place {
    my ($u,$v,$w,$component) = @_;
    my $leftend = offset($component);
    if ($v && $w) {
	if (intersectSize($u,$w) < intersectSize($u,$v)
	    && intersectSize($u,$w) < intersectSize($v,$w)) {
	    if (L($v) < L($w)) { placeLeft($u,$v); } 
	    else {  placeRight($u,$v); } 
	} else {
	    if (L($v) > L($w)) { placeLeft($u,$v); } 
	    else {  placeRight($u,$v); } 
	}
    } elsif ($v) {
	if (L($v)==$leftend) {  placeLeft($u,$v); } 
	else { placeRight($u,$v); } 
    } else { ## first clone of component
	setL($u,1);
	setR($u,cloneSize($u));
    }
    setOffset($component,L($u)) if L($u)<$leftend;
    my ($l,$r) = (L($u),R($u));
    print "placing u=$$u[0] at ($l,$r)  [v=$$v[0], w=$$w[0]]\n";
}

sub depthFirst {
    my ($u,$v,$w,$component,$clones) = @_;
    return if defined(L($u));
    addToComponent($u,$component);
    place($u,$v,$w,$component);
    my $vv;
    foreach (@$clones) {
	if (cloneOverlap($u,$_) && !defined(L($_))) {
	    depthFirst($_, $u, ($v || $vv), $component, $clones);
	    $vv = $_;
	}
    }
}

sub assignProbes {
    my ($cc,$allclones,$allccs) = @_;  ##cc = connected component of Gc
    my $id = id($cc);
    my $probes = probes($cc);
    my $leftend = offset($cc);
    my @clones =  map { $$allclones[$_] } @{clones($cc)};
    my $i=0;
    my @events = sort { $$a[0] <=> $$b[0] }
    (map { ([setL($_,L($_)-$leftend), $i],
	    [setR($_,R($_)-$leftend)+0.5, $i++]) }
     @clones);
    setOffset($cc,0);
    
    my @out = (1) x (1+$#clones);
    for (my $i=0; ; ) {
	while (@events && $events[0]->[0] <= $i) {
	    my ($ii,$c) = @{shift @events};
	    $out[$c] = 1 - $out[$c];
	};
	last unless @events;
	my $s = $probes;
	for (my $j=0; $j<=$#out; $j++) {
	    my $p = probes($clones[$j]);
	    if ($out[$j]) { $s =~ s/[$p]//g; }
	    else { $s =~ s/[^$p]//g; }
	}
	
	my @set = ("*SET*");
	foreach my $occ (@$allccs) {
	    next if ($occ == $cc);
	    my $p = probes($occ);
	    if (setSubset($p,$probes)) { setId($occ, -1) };
	    last unless (length $p <= length $s);
	    if (setSubset($p,$s)) {
		push @set, order($occ);
		$s =~ s/[$p]//g;
	    }
	}
	map {push @set,$_} split(//,$s);
	if ($#set == 1) { push @{order($cc)}, $set[1]; }
	elsif ($#set == 2) { push @{order($cc)}, [$set[1],$set[2]]; }
	else { push @{order($cc)}, \@set; }
	
	while ($i < $events[0]->[0]) { $i++; }
    }
    printComponentOrder("\nAssigning Probes to Component $id", $cc);
}


sub main {
    my @clones = ();
    while (my $in = <DATA>) {
	chomp $in;
	push @clones, newClone($#clones+1,$in) if ($in);
    }
    printCloneList("\nInput", \@clones);
    
    my @components = ();
    foreach (@clones) {
	if (!defined(cComp($_))) {
	    my $newccomp = newComponent($#components+1);
	    push @components, $newccomp;
	    depthFirst($_,undef,undef,$newccomp,\@clones);
	}
    }
    
    printCloneList("\nPlaced Clones", \@clones);
    printComponentList("Components of Gc", \@components);
    
    @components = 
	sort {(length(probes($a)) <=> length(probes($b)))
		  || ($#{clones($b)} <=> $#{clones($a)})}
    @components;
    
    map {assignProbes($_,\@clones,\@components)} @components;
    
    my $answer = ["*SET*"];
    map { push @$answer, order($_) if (id($_)>=0); } @components;
    shift @$answer if (@$answer == 2 || @$answer == 3);
    $answer = $$answer[0] if (@$answer == 1);
    print "\nFinal Answer is: \n", orderString($answer), "\n";
    print "\nOr: \n"; 
    map { print "$_\n"; } (sort @{allOrders($answer)});
    
}

main();

__END__
ABDEGI
BCDEFGHI
BDEGI
CH
CF
DG
BG
DEI
acdef
cgh
bef
af
acdfg
XYZA

