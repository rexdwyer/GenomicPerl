#########################################
package UnionFind;
#########################################
# Implements disjoint sets. 
# Uses inverted trees, union by rank, and path compression.
#########################################
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(new union find);

use strict;
use Util;

#########################################
sub new
# Creates a new union/find object.
#########################################
{
    my ($this) = @_;
    return bless [{ },{ }], (ref($this) || $this);
}

#########################################
sub union
#########################################
{
    my ($this,
	$key1, $key2) = @_;
    my ($p,$rank) = @$this;
    $key1 = $this->find($key1);
    $key2 = $this->find($key2);
    ($key1,$key2) = ($key2,$key1) if $$rank{$key1} > $$rank{$key2};
    $$rank{$key2} = max($$rank{$key2}, 1+$$rank{$key1});
    return ($$p{$key1} = $key2);
}

#########################################
sub find
#########################################
{
    my ($this,
	$key) = @_;
    my ($p,$rank) = @$this;
    my $pKey = $$p{$key};
    if (!defined($pKey)) {
	return ($$p{$key} = $key);
    } elsif ($pKey eq $key) {
	return ($key);
    } else {
        return ($$p{$key} = $this->find($pKey));
    }
}

1;
