package PrositeSuffixTree;
require SuffixTree;
@ISA = qw(SuffixTree);
use strict;

sub listMotifs {
    my ($this, $motif) = @_;
    my $pattern = lc $$motif{PA};
    return unless $pattern;
    my $description = $$motif{DE};
    $description =~ s/\n//gs;   ## remove newlines
    return "Can't handle $description\n" if $pattern =~ /\,/;
    foreach ($this->searchForPattern($pattern)) {
	print "$description seen at position $_\n";
    }
}

sub pairForm {
    my ($perlPatElem) = @_;
    my ($spec, $rep) = ($perlPatElem =~ /^([^\{]*)(\{.*\})?$/);
    $rep =~ s/[\{\}]//g;
    $rep ||= 1;
    return ($spec, $rep);
}

sub searchForPattern { 
    my ($this, $pattern) = @_;
    $pattern =~ s/\n//gs;   ## remove newlines
    $pattern =~ s/\.$//;
    $pattern =~ s/^</<-/;
    $pattern =~ s/>$/->/;
    my $patLinks = [];
    foreach (reverse split("-",$pattern)) {
	$patLinks = [pairForm(Prosite::perlizeElement($_)), $patLinks];
    }
    return ($this->searchHelp($this->{root}, $patLinks));
}

sub searchHelp {
    my ($this, $node, $patLinks) = @_;
    my $ssLen = $$node{len};
    my $ss = substr($this->{target},$$node{off},$ssLen);
#   print "At node $node representing $ss \n";
    my ($spec, $rep);
    while ($ss && @$patLinks) {
	($spec, $rep, $patLinks) = @$patLinks;
	my $pat = "$spec\{$rep\}";
#	print "Comparing $ss to $pat \n";
	if ($ss !~ s/^$pat//) {
	    return () if ($ssLen > $rep);
	    my $pat = "$spec\{$ssLen\}";
	    return () if ($ss !~ s/^$pat//);
	}
	$ssLen -= $rep;
    }
    ## Reached end of substring for this node.
    $patLinks = [$spec, -$ssLen, $patLinks] if $ssLen<0;
    if (@$patLinks==0) { return ($this->gatherStartPoints($node)); }
    $spec = $$patLinks[0];
    if ($spec =~ /[\^\.]/) {
#	my @goodKeys = grep(/^$spec$/, (keys %$node));
#	print "Spec $spec; recursive calls for keys @goodKeys\n";
	return map(($this->searchHelp($$node{$_}, $patLinks)),
		   grep(/^$spec$/, (keys %$node)));
    } 
    else {
#	my @goodKeys = grep($$node{$_},(split //,$spec));
#	print "Spec $spec; recursive calls for keys @goodKeys\n";
	return map(($this->searchHelp($$node{$_}, $patLinks)),
		   grep($$node{$_},(split //,$spec)));
    }
}
1;
