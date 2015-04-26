#!/usr/local/bin/perl
use strict;

### Constants and tables:

my $INF = 1000000;
my $DEF = -50;   ## Default terminal mismatch, used for I 

my @hairpin = (
   $INF, $INF, $INF, 410, 490, 440, 470, 500, 510, 520, 531,
        542, 551, 560, 568, 575, 582, 589, 595, 601, 606,
        611, 616, 621, 626, 630, 634, 638, 642, 646, 650);

my @bulge = (
   $INF, 390, 310, 350, 420, 480, 500, 516, 531, 543, 555,
        565, 574, 583, 591, 598, 605, 612, 618, 624, 630,
        635, 640, 645, 649, 654, 658, 662, 666, 670, 673);

my @internalLoop = (
   $INF, $INF, 410, 510, 490, 530, 570, 587, 601, 614, 625,
        635, 645, 653, 661, 669, 676, 682, 688, 694, 700,
        705, 710, 715, 720, 724, 728, 732, 736, 740, 744);


my $MLclosing = 460;
my $MLintern =10, 
my $MLbase =40;

my $lxc = 107.856;

# These numbers are Ivo's "types"
my %pairable = ("CG"=>1, "GC"=>2, "GU"=>3, "UG"=>4, "AU"=>5, "UA"=>6);  



my %dangle3 = 
    ("CGA"=>-110, "CGC"=> -40, "CGG"=>-130, "CGU"=> -60,
     "GCA"=>-170, "GCC"=> -80, "GCG"=>-170, "GCU"=>-120,
     "GUA"=>-120, "GUC"=> -50, "GUG"=>-120, "GUU"=> -70,
     "UGA"=> -80, "UGC"=> -50, "UGG"=> -80, "UGU"=> -60,
     "AUA"=> -70, "AUC"=> -10, "AUG"=> -70, "AUU"=> -10,
     "UAA"=> -80, "UAC"=> -50, "UAG"=>- 80, "UAU"=> -60);

foreach (keys(%dangle3)) { $dangle3{$_} += $MLbase; }

sub xdangle3 {  
    ### $dangle3{shift} || $INF is NOT good enough due to 0 in table.
    my $x = shift;
    (defined($dangle3{$x})) ? $dangle3{$x} : $INF ;
}

my %dangle5 =
    ("CGA"=>-50, "CGC"=> -30, "CGG"=>-20, "CGU"=> -10,
     "GCA"=>-20, "GCC"=> -30, "GCG"=>  0, "GCU"=>   0,
     "GUA"=>-20, "GUC"=> -20, "GUG"=>-20, "GUU"=> -20,
     "UGA"=>-20, "UGC"=> -20, "UGG"=>-20, "UGU"=> -20,
     "AUA"=>-30, "AUC"=> -30, "AUG"=>-40, "AUU"=> -20,
     "UAA"=>-30, "UAC"=> -10, "UAG"=>-20, "UAU"=> -20);

foreach (keys(%dangle5)) { $dangle5{$_} += $MLbase; }

sub xdangle5 {  
    ### $dangle5{shift} || $INF is NOT good enough due to 0 in table.
    my $x = shift;
    (defined($dangle5{$x})) ? $dangle5{$x} : $INF ;
}

my %stack = (
 "CGCG"=>-290,"CGGC"=>-200,"CGGU"=>-120,"CGUG"=>-190,"CGAU"=>-180,"CGUA"=>-170,
 "GCCG"=>-340,"GCGC"=>-290,"GCGU"=>-140,"GCUG"=>-210,"GCAU"=>-230,"GCUA"=>-210,
 "GUCG"=>-210,"GUGC"=>-190,"GUGU"=> -40,"GUUG"=> 150,"GUAU"=>-110,"GUUA"=>-100,
 "UGCG"=>-140,"UGGC"=>-120,"UGGU"=> -20,"UGUG"=> -40,"UGAU"=> -80,"UGUA"=> -50,
 "AUCG"=>-210,"AUGC"=>-170,"AUGU"=> -50,"AUUG"=>-100,"AUAU"=> -90,"AUUA"=> -90,
 "UACG"=>-230,"UAGC"=>-180,"UAGU"=> -80,"UAUG"=>-110,"UAAU"=>-110,"UAUA"=> -90);


my %TETRA_ENERGY = 
    ("GAAA" => -200,
     "GCAA" => -200,
     "GAGA" => -200,
     "GUGA" => -200,
     "GGAA" => -200,
     "UUCG" => -200,
     "UACG" => -200,
     "GCGA" => -200,
     "UCCG" => -200,
     "GUAA" => -200,
     "CUUG" => -200,
     "AUUU" => -200,
     "UUUA" => -200);



my %mismatchH =
    ("CGAA"=>-140,"CGAC"=>-200,"CGAG"=>-210,"CGAU"=>-190,
     "CGCA"=>-100,"CGCC"=>-110,"CGCG"=>-100,"CGCU"=>- 80,
     "CGGA"=>-210,"CGGC"=>-190,"CGGG"=>-140,"CGGU"=>-190,
     "CGUA"=>-140,"CGUC"=>-150,"CGUG"=>-140,"CGUU"=>-120,
    
     "GCAA"=>-110,"GCAC"=>-130,"GCAG"=>-200,"GCAU"=>-130,
     "GCCA"=>-110,"GCCC"=> -60,"GCCG"=> -60,"GCCU"=> -50,
     "GCGA"=>-230,"GCGC"=>-150,"GCGG"=>-140,"GCGU"=>-150,
     "GCUA"=> -80,"GCUC"=> -80,"GCUG"=> -80,"GCUU"=> -70,
     
     "GUAA"=> -80,"GUAC"=>-100,"GUAG"=>-170,"GUAU"=>-100,
     "GUCA"=> -70,"GUCC"=> -70,"GUCG"=> -70,"GUCU"=> -70,
     "GUGA"=>-150,"GUGC"=>-100,"GUGG"=>-100,"GUGU"=>-100,
     "GUUA"=> -80,"GUUC"=> -80,"GUUG"=> -80,"GUUU"=> -80,
     
     "UGAA"=>-120,"UGAC"=>-140,"UGAG"=>-200,"UGAU"=>-140,
     "UGCA"=> -90,"UGCC"=> -90,"UGCG"=> -70,"UGCU"=> -70,
     "UGGA"=>-200,"UGGC"=>-140,"UGGG"=>-130,"UGGU"=>-140,
     "UGUA"=> -90,"UGUC"=>-110,"UGUG"=> -90,"UGUU"=> -90,
     
     "AUAA"=> -80,"AUAC"=>-100,"AUAG"=>-170,"AUAU"=>-100,
     "AUCA"=> -70,"AUCC"=> -70,"AUCG"=> -70,"AUCU"=> -70,
     "AUGA"=>-150,"AUGC"=>-100,"AUGG"=>-100,"AUGU"=>-100,
     "AUUA"=> -80,"AUUC"=> -80,"AUUG"=> -80,"AUUU"=> -80,
     
     "UAAA"=>-100,"UAAC"=> -80,"UAAG"=>-180,"UAAU"=> -90,
     "UACA"=> -70,"UACC"=> -60,"UACG"=> -30,"UACU"=> -50,
     "UAGA"=>-180,"UAGC"=> -90,"UAGG"=>-120,"UAGU"=> -90,
     "UAUA"=>- 30,"UAUC"=> -60,"UAUG"=> -30,"UAUU"=> -50);


my %mismatchI = 
    ("CGAA"=>-150,"CGAC"=>-150,"CGAG"=>-270,"CGAU"=>-190,
     "CGCA"=>-150,"CGCC"=>-150,"CGCG"=>-100,"CGCU"=>-150,
     "CGGA"=>-270,"CGGC"=>-190,"CGGG"=>-150,"CGGU"=>-190,
     "CGUA"=>-140,"CGUC"=>-150,"CGUG"=>-140,"CGUU"=>-250,
     
     "GCAA"=>-150,"GCAC"=>-150,"GCAG"=>-270,"GCAU"=>-130,
     "GCCA"=>-150,"GCCC"=>-150,"GCCG"=> -60,"GCCU"=>-150,
     "GCGA"=>-270,"GCGC"=>-150,"GCGG"=>-150,"GCGU"=>-150,
     "GCUA"=> -80,"GCUC"=>-150,"GCUG"=> -80,"GCUU"=>-250,
     
     "GUAA"=>-150,"GUAC"=>-150,"GUAG"=>-270,"GUAU"=>-130,
     "GUCA"=>-150,"GUCC"=>-150,"GUCG"=> -60,"GUCU"=>-150,
     "GUGA"=>-270,"GUGC"=>-150,"GUGG"=>-150,"GUGU"=>-150,
     "GUUA"=> -80,"GUUC"=>-150,"GUUG"=> -80,"GUUU"=>-250,
     
     "UGAA"=>-100,"UGAC"=>-100,"UGAG"=>-220,"UGAU"=> -40,
     "UGCA"=>-100,"UGCC"=>-100,"UGCG"=>  20,"UGCU"=>-100,
     "UGGA"=>-220,"UGGC"=> -40,"UGGG"=>-100,"UGGU"=> -40,
     "UGUA"=>  20,"UGUC"=>-100,"UGUG"=>  20,"UGUU"=>-200,
     
     "AUAA"=>-100,"AUAC"=>-100,"AUAG"=>-220,"AUAU"=> -50,
     "AUCA"=>-100,"AUCC"=>-100,"AUCG"=> -20,"AUCU"=>-100,
     "AUGA"=>-220,"AUGC"=> -50,"AUGG"=>-100,"AUGU"=> -50,
     "AUUA"=> -30,"AUUC"=>-100,"AUUG"=> -30,"AUUU"=>-200,
     
     "UAAA"=>-100,"UAAC"=>-100,"UAAG"=>-220,"UAAU"=> -40,
     "UACA"=>-100,"UACC"=>-100,"UACG"=>  20,"UACU"=>-100,
     "UAGA"=>-220,"UAGC"=> -40,"UAGG"=>-100,"UAGU"=> -40,
     "UAUA"=>  20,"UAUC"=>-100,"UAUG"=>  20,"UAUU"=>-200);


my $MAX_NINIO = 300;
my $FNinio2 = 30;


### Program variables

my @f5;
my @s;
my $jamesRule = 1;

my @c;
my @fML;
my @DML;
my $rna;


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


sub fold {
    my ($s) = @_;
    @s = ('X', split(//, $s));
    my ($i, $j, $p, $q);
    
    my $slen = length $s;
    
    for ($j = 1; $j<=$slen; $j++) {
	for ($i=max($j-3,1);  $i<=$j;  $i++) {
	    $c[$i][$j] = $fML[$i][$j] = $DML[$i][$j] = $INF;
	}       
    }
    
    for ($i = $slen-4; $i >= 1; $i--) { ### i,j in [1..$slen] 
	for ($j = $i+4; $j <= $slen; $j++) {
	    my $pair = $s[$i] . $s[$j];
	    my $rpair = $s[$j] . $s[$i];
	    if (!($pairable{$pair})) {
		$c[$i][$j] = $INF;
	    } else {
	    
		$c[$i][$j] = hairpin($i, $j, $pair, $s);
		
		### Try internal-loop decomposition:
		
		my $pmax = min($j-5,$i+31);  # limit length of internal loops
		for ($p = $i+1; $p <= $pmax ; $p++) {
		    my $qmin = max($p+4, $j - 32 + ($p-$i)); # limit length
		    for ($q = $qmin; $q < $j; $q++) {
			my $pair2 = $s[$p] . $s[$q];
			if ($pairable{$pair2}) {
			    $c[$i][$j] =
				min($c[$i][$j],
				     $c[$p][$q]+
				     internalLoop($i,$j,$p,$q,$pair,$pair2));
			}
		    } 
		} 
		### Try multi-loop decomposition:
		
		my $MLenergy = $MLclosing+$MLintern+
		    min($DML[$i+1][$j-1],
			 $DML[$i+2][$j-1] + $dangle3{$rpair.$s[$i+1]},
			 $DML[$i+1][$j-2] + $dangle5{$rpair.$s[$j-1]},
			 $DML[$i+2][$j-2] + $dangle5{$rpair.$s[$j-1]}
			 + $dangle3{$rpair.$s[$i+1]});
		$c[$i][$j] = min($MLenergy,$c[$i][$j]);
	    }  # end else
	    
	    
	    ### modular decomposition:
	    
	    $DML[$i][$j] = min( map {$fML[$i][$_]+$fML[$_+1][$j]} ($i+4..$j-5) );
	    
	    ### free ends ?
	    $fML[$i][$j] =   ### substring energy
		min($DML[$i][$j],
		     $fML[$i+1][$j]+$MLbase,
		     $fML[$i][$j-1]+$MLbase,
		     $c[$i][$j]+$MLintern,
		     $c[$i+1][$j]+xdangle5($s[$i+1].$s[$j].$s[$i])+$MLintern,
		     $c[$i][$j-1]+xdangle3($s[$i].$s[$j-1].$s[$j])+$MLintern,
		     $c[$i+1][$j-1]+xdangle5($s[$i+1].$s[$j-1].$s[$i])
		     +xdangle3($s[$i+1].$s[$j-1].$s[$j])+$MLintern);

#	    print "$i $s[$i] $j $s[$j] $c[$i][$j] $fML[$i][$j] $DML[$i][$j]\n";
	}  ### end for (j=...  
	print STDERR "finished i=$i...\n";
    }  ###  end for(i=...   
    
    ### calculate energies of 5' and 3' fragments 
    
    $f5[4]=0;
    for ($j=5; $j<=$slen; $j++) {
	$f5[$j] = $f5[$j-1];
	$f5[$j] = min($f5[$j], $c[1][$j])
	    if ($pairable{$s[1].$s[$j]});
	
	$f5[$j] = min($f5[$j], 
		       $c[1][$j-1]+$dangle3{$s[1].$s[$j-1].$s[$j]}-$MLbase)
	    if ($pairable{$s[1].$s[$j-1]});
	
	for ($i=$j-4; $i>1; $i--) {
	    my $pair = $s[$i].$s[$j];
	    if ($pairable{$pair}) {
		$f5[$j] = 
		    min($f5[$j], 
			 $f5[$i-1]+$c[$i][$j],
			 $f5[$i-2]+$c[$i][$j]+$dangle5{$pair.$s[$i-1]}-$MLbase);
	    }
	    $pair = $s[$i].$s[$j-1];
	    if ($pairable{$pair}) {
		$f5[$j] = 
		    min($f5[$j], 
			 $f5[$i-1]+$c[$i][$j-1]+$dangle3{$pair.$s[$j]}-$MLbase,
			 $f5[$i-2]+$c[$i][$j-1] - 2*$MLbase
			 + $dangle5{$pair.$s[$i-1]}+$dangle3{$pair.$s[$j]});
	    }
	}
#	print "f5[$j]=$f5[$j]\n";
    }
}

###---------------------------------------------------------------------------

sub hairpin {
    my ($i, $j, $pair, $s) = @_;
    my $energy = ($j-$i-1 <= 30) ?
	$hairpin[$j-$i-1] : $hairpin[30]+ int($lxc*log(($j-$i-1)/30.));    
    return $energy if ($j-$i-1 == 3);
    $energy += $TETRA_ENERGY{substr($s,$i,4)} if ($j-$i-1 == 4);
    return $energy + $mismatchH{$pair.$s[$i+1].$s[$j-1]};
}

###---------------------------------------------------------------------------

sub internalLoop {
    my ($i, $j, $p, $q, $pair, $pair2) = @_;
    ### compute energy of degree 2 loop (stack bulge or interior) 
    my $n1 = $p-$i-1;
    my $n2 = $j-$q-1;
    
    ($n1,$n2) = ($n2,$n1) if ($n1>$n2);
    
    return $stack{$pair.$pair2} if ($n2 == 0);
    
    return $bulge[$n2]+$stack{$pair.$pair2} if (($n2 == 1)&&($n1==0));

    return $bulge[$n2] if ($n1 == 0);

    ### special case for loop size 2 
    return 80 if (($n1+$n2==2) && $jamesRule);

    ### interior loop 

    my $energy = ($n1+$n2<=30)?($internalLoop[$n1+$n2]):
	($internalLoop[30]+int($lxc*log(($n1+$n2)/30.)));
    
    $energy += min($MAX_NINIO, ($n2-$n1)*$FNinio2);
    
    $energy += $mismatchI{$pair.$s[$i+1].$s[$j-1]}+
	$mismatchI{(scalar reverse $pair2).$s[$q+1].$s[$p-1]};
    
    return $energy;
}

sub f5Trace {
    my ($j) = @_;
    return ("." x $j) if ($j < 5);
    return f5Trace($j-1) . "."  if ($f5[$j] == $f5[$j-1]) ;
    my $jj;
    for (my $k=$j-4; $k>1; $k--) {
	$jj = $k-1;
	my $pair = $s[$k].$s[$j-1];
	if($pairable{$pair}) {
	    
	    return (f5Trace($k-1) . foundPair($k,$j-1) . ".")
		if ($f5[$j] == $f5[$k-1]+$c[$k][$j-1]+$dangle3{$pair.$s[$j]}-$MLbase);
	    
	    return (f5Trace($k-2) . "." . foundPair($k,$j-1) . ".")
		    if ($f5[$j] == $f5[$k-2]+$c[$k][$j-1]+ $dangle5{$pair.$s[$k-1]}+$dangle3{$pair.$s[$j]}-2*$MLbase);
	}
	$pair = $s[$k].$s[$j];
	if($pairable{$pair}) {
	    return (f5Trace($k-1) . foundPair($k,$j)) 
		if ($f5[$j] == $f5[$k-1]+$c[$k][$j]);
	    return (f5Trace($k-2) . "." . foundPair($k,$j))
		if ($f5[$j] == $f5[$k-2]+$c[$k][$j]+$dangle5{$pair.$s[$k-1]}-$MLbase);
	}
    }
    my $pair= $s[1].$s[$j-1];
    return (foundPair(1,$j) . ".")
	if ($pairable{$pair} && $f5[$j] == $c[1][$j-1]+$dangle3{$pair.$s[$j]}-$MLbase);
    return foundPair(1,$j) if ($f5[$j] == $c[1][$j]);
    die "f5Trace failed\n";
}


sub fMLTrace {
    my ($i,$j) = @_;
    return if ($j < $i+4);
    my $fij = $fML[$i][$j];

    return "." . fMLTrace($i+1,$j) 
	if ($fij == $fML[$i+1][$j]+$MLbase);
    return fMLTrace($i,$j-1) . "." 
	if ($fij == $fML[$i][$j-1]+$MLbase);
    return foundPair($i,$j) 
	if ($fij == $c[$i][$j] + $MLintern);
    return "." . foundPair($i+1,$j) 
	if ($fij == $c[$i+1][$j] + $dangle5{$s[$i+1].$s[$j].$s[$i]} + $MLintern);
    return foundPair($i,$j-1) . "."
	if ($fij == $c[$i][$j-1] + $dangle3{$s[$i].$s[$j-1].$s[$j]} + $MLintern);

    my $pair = $s[$i+1].$s[$j-1];
    return "." . foundPair($i+1,$j-1) . "." 
	if ($fij == $c[$i+1][$j-1] + $dangle5{$pair.$s[$i]} + $dangle3{$pair.$s[$j]}
	    +  $MLintern);

    for (my $k = $i+4; $k <= $j-5; $k++) {
	return fMLTrace($i,$k) . fMLTrace($k+1,$j)
	    if ($fij == ($fML[$i][$k]+$fML[$k+1][$j]));
    }
    die "fMLTrace failed\n";
}

sub foundPair {
    my ($i,$j) = @_;
	
#   print "basepair($i,$j)\n";
    
    my $pair = $s[$i].$s[$j];
    
    return ( "(" . ("." x ($j-$i-1)) . ")")
	if ($c[$i][$j] == hairpin($i, $j, $pair, $rna));
    
    my $pmax = min($j-5,$i+31);  # limit length of internal loops
    for (my $p = $i+1; $p <= $pmax ; $p++) {
	my $qmin = max($p+4, $j - 32 + ($p-$i)); # limit length
	for (my $q = $qmin; $q < $j; $q++) {
	    my $pair2 = $s[$p].$s[$q];
	    return ( "(" . ("."x($p-$i-1)) . foundPair($p,$q) . ("."x($j-$q-1)) . ")" )
		if ($pairable{$pair2} &&
		    ($c[$i][$j] == $c[$p][$q]+internalLoop($i,$j,$p,$q,$pair,$pair2)));
	}
    } 
    
    my $mm = $MLclosing+$MLintern;
    $pair = $s[$j].$s[$i];
    
    for (my $k = $i+5; $k <= $j-6; $k++) {
	return "(" . fMLTrace($i+1,$k) . fMLTrace($k+1,$j-1) . ")"
	    if ($c[$i][$j] == $fML[$i+1][$k]+$fML[$k+1][$j-1]+$mm);

	return "(." . fMLTrace($i+2,$k) . fMLTrace($k+1,$j-1) . ")"
	    if ($c[$i][$j] == ($fML[$i+2][$k]+$fML[$k+1][$j-1]+$mm+
			   $dangle3{$pair.$s[$i+1]}));

	return "(" . fMLTrace($i+1,$k) . fMLTrace($k+1,$j-2) . ".)"
	    if ($c[$i][$j] == ($fML[$i+1][$k]+$fML[$k+1][$j-2]+$mm+
			   $dangle5{$pair.$s[$j-1]})) ;

	return "(." . fMLTrace($i+2,$k) . fMLTrace($k+1,$j-2) . ".)"
	    if ($c[$i][$j] == ($fML[$i+2][$k]+$fML[$k+1][$j-2]+$mm+
			   $dangle3{$pair.$s[$i+1]}+$dangle5{$pair.$s[$j-1]}));
    }
    die "foundPair failed\n";
}


$rna = <STDIN>;
print $rna;
chomp($rna);
fold($rna);
print f5Trace(length $rna), "\n";
# print $f5[length $rna]/100., "\n";
