#!/usr/local/bin/perl

# An implementation of the fast, space-saving algorithm 
# for global alignment of similar strings, S&M p. 68.


# Except that I haven't done it yet...


sub max {
    local($m,@l) = @_;
	foreach $x (@l) { $m = $x if ($x > $m) }
	return $m;
}

sub compactScore {
	local($s,$t) = @_;
    local $i,$j,@a,$old, $temp;
	local $tlen = length($t);
	undef @a; $#a= $tlen;
    foreach $j (0..$tlen) { $a[$j] = $g * $j; }
    foreach $i (1..length($s)) {
        $old = $a[0];
		$a[0] = $g * $i;
		foreach $j (1..$tlen) {
            $temp = $a[$j];
            $m = (substr($s,$i-1,1) eq substr($t,$j-1,1))?1:-1;

			$a[$j] = &max($a[$j]+$g, $a[$j-1]+$g, $old+$m);
			$old = $temp;
		}
	}
    return(@a);
}

sub revString() {
	return join('',reverse(split(//,$_[0])));
}

sub recursiveAlign {
	local($s,$t) = @_;
    local $slen = length($s);
    local $tlen = length($t);
	if ($slen == 0) {
		return( ("-" x $tlen) . "\n$t\n" . (-2*$tlen) );
	} elsif ($tlen == 0) {
		return( "$s\n" . ("-" x $slen) . "\n" . (-2*$slen) );
	} elsif (($slen == 1) && ($tlen == 1)) {
		local $score = ($s eq $t)?1:-1;
		return ("$s\n$t\n$score");
	} else {
		local $mid = int($slen/2);
		local $sBeg = substr($s,0,$mid);
		local $sEnd = substr($s,$mid, $slen-$mid);
        local @aBeg = &compactScore($sBeg, $t);
        local @aEnd = reverse(&compactScore(&revString($sEnd), 
											  &revString($t)));
		local @a;
		local $k;
        foreach $k (0..$tlen) { $a[$k] = $aBeg[$k]+$aEnd[$k]; }
		local $kval = max(@a);
        for (local $kk=0; $a[$kk] != $kval; $kk++) {$k=$kk; }
        $k++;
		local ($ssBeg,$ttBeg,$scoreBeg) = 
			split("\n", &recursiveAlign($sBeg, substr($t,0,$k)));
		local ($ssEnd,$ttEnd,$scoreEnd) = 
			split("\n", &recursiveAlign($sEnd, substr($t,$k,$tlen-$k)));
        local $score = $scoreBeg + $scoreEnd;
		return( "$ssBeg$ssEnd\n$ttBeg$ttEnd\n$score" );
		}
}

#MAIN
{
	$s = <DATA>; chomp($s);
	$t = <DATA>; chomp($t);
	$g = (-2);
	print "\n", &recursiveAlign($s,$t), "\n";
}


__END__
RECONSTRUCTION
UNCONDITIONALLY

THANKSGIVING
CHANGELING
