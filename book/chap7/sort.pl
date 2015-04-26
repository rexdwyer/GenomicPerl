
my @L = qw(Five Four Three Two One five four three two one);
my %h = (Five=>5, Four=>4, Three=>3, Two=>2, One=>1,
	 five=>5, four=>4, three=>3, two=>2, one=>1);

print sort { $a cmp $b } @L;
print "\n";  ## sorts alphabetically
  ## prints: FiveFourOneThreeTwofivefouronethreetwo

print sort { lc($a) cmp lc($b) } @L;
print "\n";  ## sorts alphabetically, ignoring case
  ## prints: FivefiveFourfourOneoneThreethreetwoTwo

print sort { length($b) <=> length($a) } @L;
print "\n";  ## sorts longest string to shortest string
  ## prints: ThreethreeFiveFourfivefourOneTwotwoone

print sort { substr($a,-2) cmp substr($b,-2)} @L;
print "\n";  ## sorts alphabetically by last 2 letters
  ## prints: ThreethreeOneoneFourfourfiveFivetwoTwo
  ##            --   -- -- --  --  --  --  -- -- -- 

print sort { $h{$a} <=> $h{$b} || $b cmp $a } @L;
print "\n";  ## sorts first by value in hash, then lexicographically.
  ## prints: oneOnetwoTwothreeThreefourFourfiveFive
