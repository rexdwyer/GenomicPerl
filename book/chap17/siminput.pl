#!/usr/bin/perl -I . -I /home/dwyer/book/appC -I /home/dwyer/book/perl

use strict;

main();
exit();

sub main {
    $|=1;
    my $chromoLen = promptForParameter("chromoLen", "1E6");
    my $cloneLen = promptForParameter("cloneLen", "1E4");
    my $coverage = promptForParameter("coverage", "10");
    my $siteLen = promptForParameter("siteLen", "4");
    ## default restriction enzyme is a 4-cutter.
    my $sites = generateRestrictionSites($chromoLen, $siteLen);
    printf "\# There are %d sites total.\n", scalar @$sites;
    my $protoclones = generateProtoclones($chromoLen, $cloneLen, $coverage,$sites);
    my $probes = identifyProbes($protoclones);
    printf "\# There are %d probes total.\n", scalar @$probes;
    my $probesForClone = filterProtoclones($protoclones, $cloneLen,$probes);
    filterAndRenumberProbes($probesForClone); 
    printClones($probesForClone);
}

sub promptForParameter {
    my ($string, $default) = @_;
    print STDERR "Enter value for $string [$default]: ";
    my $value = <STDIN>; 
    chomp($value);
    $value = int($value || $default);
    print "\# Parameter value for $string: $value\n";
    return $value;
}

sub generateRestrictionSites {
    my ($chromoLen, $siteLen) = @_;
    my @sites;
    my $noCut = $siteLen;
    foreach (1..$chromoLen) {
	next if (rand() < 0.25 && --$noCut);
	push @sites, $_ if (!$noCut);
	$noCut = $siteLen;
    }
    return \@sites;
}

sub generateProtoclones {
    my ($chromoLen, $cloneLen, $coverage,$sites) = @_;
    my $avgFrag = $chromoLen / (1+@$sites);
    my $cutProb = $avgFrag / $cloneLen;
    my @protoclones;
    foreach (1..$coverage) {
	my @cuts = grep(rand() < $cutProb, @$sites);
	print "\# ", scalar @cuts, " cuts in copy $_\n";
	my $last = shift @cuts;
	foreach (@cuts) {
	    push @protoclones, [$last, $_];
	    $last = $_+1;
	}
    }
    ## Now shuffle the protoclones
    my $c = $#protoclones;
    while (--$c) {
	my $r = int(rand($c));
	@protoclones[$c,$r] = @protoclones[$r,$c];
    }
    return \@protoclones;
}

sub identifyProbes {
    my ($protoclones) = @_;
    return 
	[sort 
	 {$a<=>$b}
	 keys %{{
	     map((($$_[0]=>1),($$_[1]=>1)), @$protoclones)
	     }}
	 ];
}

sub filterProtoclones {
    my  ($protoclones, $cloneLen,$probes) = @_;
    my %probesForClone;
    my $savedProtoclone;
    my ($cCnt, $xCnt, $id);
    while (@$protoclones){
	my $protoclone = pop @$protoclones;
	my ($lo, $hi) = @$protoclone;
	next if ($hi-$lo) < (0.75*$cloneLen);  ## too small
	next if ($hi-$lo) > (1.75*$cloneLen);  ## too big
	if (($hi-$lo) < $cloneLen && !$savedProtoclone) {  ## save for chimera
	    $savedProtoclone = $protoclone;
	    next;
	}
	my $firstProbe = binarySearch($lo, $probes);
	my $lastProbe = binarySearch($hi, $probes);
	my @pList = @$probes[$firstProbe..$lastProbe];
	if (($hi-$lo) > $cloneLen) {  ## output is current clone
	    $cCnt++;
	    $id = "c$cCnt";
	} else {  ## output is chimera = current + saved
	    $xCnt++;
	    $id = "x$xCnt";
	    my ($lo, $hi) = @$savedProtoclone;
	    undef($savedProtoclone);
	    my $firstProbe = binarySearch($lo, $probes);
	    my $lastProbe = binarySearch($hi, $probes);
	    push @pList, @$probes[$firstProbe..$lastProbe];
	    ## Remove duplicates from pList.  They could occur if a chimera
	    ## is formed from protoclones that overlap.
	    @pList = sort keys %{{map(($_=>1), @pList)}};
	}
	$probesForClone{$id} = \@pList;
    }
    print "\# There are $cCnt simple clones and $xCnt chimerae.\n";
    return \%probesForClone;
}

sub binarySearch {
    my ($key, $arr) = @_;
    my ($lo, $hi) = (0, $#$arr);
    return $hi if $$arr[$hi] == $key;
    while ($lo < $hi) {
	my $mid = int(($lo+$hi)/2);
	if ($$arr[$mid] < $key) { $lo = $mid; }
	elsif ($$arr[$mid] > $key) { $hi = $mid; }
	else { return $mid; } 
    }
}

sub filterAndRenumberProbes {
    my ($probesForClone) = @_; ## Filter out indistinguishable probes.
    my %count; my %probeId;
    foreach (keys %$probesForClone) {  ## for each clone
	foreach my $p1 (@{$$probesForClone{$_}}) {  ## for each of its probes
	    $probeId{$p1}++;
	    foreach my $p2 (@{$$probesForClone{$_}}){ ## for each probe pair
		last if $p1 eq $p2;
		$count{$p1,$p2}++;
	    }
	}
    }
    
    my @probes = keys %probeId;
    foreach my $p1 (@probes) {
	foreach my $p2 (@probes) {
	    last if $p1 eq $p2;
	    my $common = $count{$p1,$p2} + $count{$p2,$p1};
	    if ($probeId{$p1} == $common && $probeId{$p2} == $common) {
		delete $probeId{$p1};
		last;
	    }
	}
    }
    @probes = sort {$a<=>$b} keys %probeId;
    my @breaks;
    foreach my $i (0..$#probes-1) {
	my ($p1, $p2) = @probes[$i, $i+1];
	push @breaks, $p1 if $count{$p1,$p2} + $count{$p2,$p1} == 0;
    }
    
    my $i;
    foreach (@probes) { $probeId{$_} = ++$i; }
    foreach my $clone (keys %$probesForClone) {
	$$probesForClone{$clone} = 
	    [sort {$a<=>$b}
	     map($probeId{$_}||(), @{$$probesForClone{$clone}})
	     ];
    }
    print "\# There are $i distinguishable probes.\n";
    warn "There are $i distinguishable probes.\n";
    my $breaks = join(' ',@probeId{@breaks});
    my $numBreaks = @breaks;
    print "\# There are $numBreaks breaks after probes $breaks.\n";
    warn "There are $numBreaks breaks after probes $breaks.\n";
}

sub printClones {
    my ($probesForClone) = @_;
    foreach my $cx ("c", "x") {
	for (my $i=1; $$probesForClone{$cx.$i}; $i++) {
	    my $cxi = $cx.$i;
	    print "\n$cxi:\n";
	    my $s = join(" ", @{$$probesForClone{$cxi}}) . " ";
	    while ($s =~ s/(.{1,70}) +//) {  print "  $1\n"; }
	    print "  $s\n" if $s;
	}
    }
}

