#!/usr/bin/perl

use strict; 

my @structure;
my @bases;

my $INF = 2000000000;

sub max {
    my $ans = -$INF;
    foreach (@_) { $ans = $_ if $_>$ans; };
    $ans;
}

sub min {
    my $ans = $INF;
    foreach (@_) { $ans = $_ if $_<$ans; };
    $ans;
}


sub begin_PostScriptPicture {
    print "/rnaline { pt setlinewidth moveto lineto stroke} def\n";
    print "/rnabase { moveto -0.12 -0.12 rmoveto show} def\n";
    print "/rnapicture {\n\n";
}

my ($minx, $maxx, $miny, $maxy);

sub end_PostScriptPicture {
    print "\n} def\n\n";
    my $scale = 0.95 * (8.5*72) / max($maxx-$minx, $maxy-$miny);
    my $xorigin = ((8.5*72) * -$minx / ($maxx-$minx))  || "0";
    my $yorigin = ((8.5*72) * -$miny / ($maxy-$miny)) || "0";
    print "/pt {$scale div} def\n\n";
    print "/Helvetica findfont 8 pt scalefont setfont\n";
    print "$xorigin $yorigin translate\n";
    print "$scale dup scale\n";
    print "rnapicture showpage\n";
}


sub drawLine{
    my ($x1, $y1, $x2, $y2, $thick) = @_;
    $maxx = max($maxx, $x1, $x2);
    $minx = min($minx, $x1, $x2);
    $maxy = max($maxy, $y1, $y2);
    $miny = min($miny, $y1, $y2);
    ($x1, $y1, $x2, $y2) = (0.8 * $x1 + 0.2 * $x2,
			    0.8 * $y1 + 0.2 * $y2,
			    0.2 * $x1 + 0.8 * $x2,
			    0.2 * $y1 + 0.8 * $y2);
    print "$x1 $y1 $x2 $y2 $thick rnaline\n";
}

sub drawBase{
    my ($x, $y, $b) = @_;
    $b .= "'" if ("53" =~ $b);
    $b .= "---------------------------HERE" if ("53" =~ $b);
    print "($b) $x $y rnabase\n";
}

sub drawRna {
    my ($l, $r, $lx, $ly, $rx, $ry) = @_;
    my $level=0;
    my $stemcount=1;
    $stemcount=0 if $bases[$l]==5;
    my $maxdotrun=0;
    my $dotrun=0;
    for (my $i=$l+1; $i<$r; $i++) {
	$level-- if ($structure[$i] eq ")");
	if ($level==0) {
	    if ($structure[$i] eq ".") { 
		$maxdotrun = max($maxdotrun,++$dotrun);
	    } elsif ($structure[$i] eq "(") { 
		$stemcount++;
		$dotrun = 0;
	    }
 	}
	$level++ if ($structure[$i] eq "(");
    }

    my $count = (2+$maxdotrun)*$stemcount;
    my $theta = 2 * 3.14159 / $count;
    my $rad = 1 / (2*sin($theta/2));            ## radius 
    my $h = $rad * cos($theta/2);
    my ($cx, $cy) = ((($lx+$rx)/2.0)+$h*($ly-$ry),    ## center of circle
		     (($ly+$ry)/2.0)+$h*($rx-$lx)); 
    my $alphal = atan2($ly-$cy,$lx-$cx);
    my $alphar = atan2($ry-$cy,$rx-$cx);

    my ($ii,$xx,$yy) = ($l,$lx,$ly);
    
    my $stem=1;
    $dotrun=0;
    my $phi=0;
    my $alpha=$alphal;
    my ($x,$y);

    for (my $i=$l+1; $i<=$r; $i++) {
	$level-- if ($structure[$i] eq ")");
	if ($level==0) {
	    if ($structure[$i] eq "(") {
		$alpha = $alphar - (2*3.14159 * $stem/$stemcount);
		($x,$y) = ($cx+$rad*cos($alpha), $cy+$rad*sin($alpha));
	    } elsif ($structure[$i] eq ")") {
		$alpha = $alphal - (2*3.14159 * $stem/$stemcount);
		($x,$y) = ($cx+$rad*cos($alpha), $cy+$rad*sin($alpha));
		drawRna($ii,$i,$xx,$yy, $x, $y);
		$phi=0; $stem++;
	    } else {
		if ($phi == 0) {
		    my $j=0;
		    while ($structure[$i+$j] eq ".") { $j++ }
		    $phi = $theta * ($maxdotrun+1) / ($j+1);
		}
		$alpha -= $phi;
		($x,$y) = ($cx+$rad*cos($alpha), $cy+$rad*sin($alpha));
	    }
	    drawLine($xx,$yy,$x,$y,0);
	    drawBase($xx,$yy, $bases[$ii]);
	    ($xx,$yy)=($x,$y);
	    $ii = $i;
	}
	$level++ if ($structure[$i] eq "(");
    }
    drawLine($xx,$yy,$rx,$ry,0);
    drawBase($xx,$yy, $bases[$r-1]);
    drawBase($rx,$ry, $bases[$r]);
    my %bonds = (GU=>1,UG=>1,AU=>2,UA=>2,CG=>3,GC=>3);
    drawLine($lx,$ly,$rx,$ry,$bonds{$bases[$l].$bases[$r]}+0) if ($lx || $ly);
}

sub drawRnaStructure {
    my ($basestring,$structurestring) = @_;
    @bases = split(//, 5 . $basestring . 3);
    @structure = split(//, "($structurestring)");
    begin_PostScriptPicture();
    drawRna(0, $#structure, 0,0, 1,0);
    end_PostScriptPicture();
}


my $basestring = <STDIN>;
chomp($basestring);
my $structurestring = <STDIN>;
chomp($structurestring);
drawRnaStructure($basestring,$structurestring);
