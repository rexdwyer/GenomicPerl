
use strict;
 
print <<EOF;
\\documentclass{book}
\\author{Rex}
\\title{pseudoknot}
\\begin{document}
\\setlength{\\unitlength}{1cm}
\\begin{picture}(8,6)
EOF


ns(0,5);
en(0,4);
ew(1,4);
ew(2,4);
ew(3,4);
nw(4,4);
se(4,5);
ew(5,5);
ew(6,5);
sw(7,5);
ns(7,4);
ns(7,3);
nw(7,2);
ew(6,2);
ew(5,2);
en(4,2);
sw(4,3);
ew(3,3);
ew(2,3);
se(1,3);
ns(1,2);
ns(1,1);
en(1,0);
ew(2,0);
ew(3,0);
nw(4,0);
se(4,1);
ew(5,1);
ew(6,1);
sw(7,1);
ns(7,0);
stitch(2,3.5);
stitch(3,3.5);
stitch(5,1.5);
stitch(6,1.5);

print <<EOF;
\\end{picture}
\\end{document}
EOF


sub ns {
    my ($x, $y) = @_;
    my ($y1,$y2) = ($y-0.5, $y+0.5);
    print "\\put($x,$y1){\\line(0,1){1}}\n";
    print "\\multiput($x,$y1)(0,0.25){5}{\\circle*{0.1}}\n";
}

sub ew {
    my ($x, $y) = @_;
    my ($x1,$x2) = ($x-0.5, $x+0.5);
    print "\\put($x1,$y){\\line(1,0){1}}\n";
    print "\\multiput($x1,$y)(0.25,0){5}{\\circle*{0.1}}\n";
}


sub en {
    my ($x, $y) = @_;
    my ($xc,$yc) = ($x+0.5, $y+0.5);
    my ($x0, $y0) = ($xc, $yc - 0.5);
    my ($x30, $y30) = ($xc - 0.25, $yc - sqrt(3.0)/4);
    my ($x60, $y60) = ($xc - sqrt(3.0)/4, $yc - 0.25);
    my ($x90, $y90) = ($xc - 0.5, $yc);
    print "\\put($xc,$yc){\\oval(1,1)[bl]}\n";
    print "\\put($x0,$y0){\\circle*{0.1}}\n";
    print "\\put($x30,$y30){\\circle*{0.1}}\n";
    print "\\put($x60,$y60){\\circle*{0.1}}\n";
    print "\\put($x90,$y90){\\circle*{0.1}}\n";
}

sub nw {
    my ($x, $y) = @_;
    my ($xc,$yc) = ($x-0.5, $y+0.5);
    my ($x0, $y0) = ($xc, $yc - 0.5);
    my ($x30, $y30) = ($xc + 0.25, $yc - sqrt(3.0)/4);
    my ($x60, $y60) = ($xc + sqrt(3.0)/4, $yc - 0.25);
    my ($x90, $y90) = ($xc + 0.5, $yc);
    print "\\put($xc,$yc){\\oval(1,1)[br]}\n";
    print "\\put($x0,$y0){\\circle*{0.1}}\n";
    print "\\put($x30,$y30){\\circle*{0.1}}\n";
    print "\\put($x60,$y60){\\circle*{0.1}}\n";
    print "\\put($x90,$y90){\\circle*{0.1}}\n";
}

sub se {
    my ($x, $y) = @_;
    my ($xc,$yc) = ($x+0.5, $y-0.5);
    my ($x0, $y0) = ($xc, $yc + 0.5);
    my ($x30, $y30) = ($xc - 0.25, $yc + sqrt(3.0)/4);
    my ($x60, $y60) = ($xc - sqrt(3.0)/4, $yc + 0.25);
    my ($x90, $y90) = ($xc - 0.5, $yc);
    print "\\put($xc,$yc){\\oval(1,1)[tl]}\n";
    print "\\put($x0,$y0){\\circle*{0.1}}\n";
    print "\\put($x30,$y30){\\circle*{0.1}}\n";
    print "\\put($x60,$y60){\\circle*{0.1}}\n";
    print "\\put($x90,$y90){\\circle*{0.1}}\n";
}


sub sw {
    my ($x, $y) = @_;
    my ($xc,$yc) = ($x-0.5, $y-0.5);
    my ($x0, $y0) = ($xc, $yc + 0.5);
    my ($x30, $y30) = ($xc + 0.25, $yc + sqrt(3.0)/4);
    my ($x60, $y60) = ($xc + sqrt(3.0)/4, $yc + 0.25);
    my ($x90, $y90) = ($xc + 0.5, $yc);
    print "\\put($xc,$yc){\\oval(1,1)[tr]}\n";
    print "\\put($x0,$y0){\\circle*{0.1}}\n";
    print "\\put($x30,$y30){\\circle*{0.1}}\n";
    print "\\put($x60,$y60){\\circle*{0.1}}\n";
    print "\\put($x90,$y90){\\circle*{0.1}}\n";
}

sub stitch {
    my ($x, $y) = @_;
    my ($x1,$y1) = ($x-0.5, $y-0.5);
    print "\\multiput($x1,$y1)(0.25,0){5}{\\line(0,1){1}}\n";
}
