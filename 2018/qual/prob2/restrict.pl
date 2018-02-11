#!/usr/bin/perl -l

use warnings;
use strict;

sub makeAmbig {
  my ($seq) = @_;
  $seq =~ s/\s+$//;
  my %ambigMap = (
    AAA => "AA#", AAC => "AA#", AAG => "AA#", AAU => "AA#",
    ACA => "AC#", ACC => "AC#", ACG => "AC#", ACU => "AC#",
    AGA => "#G#", AGC => "1!#", AGG => "#G#", AGU => "1!#",
    AUA => "AU#", AUC => "AU#", AUG => "AUG", AUU => "AU#",
    CAA => "CA#", CAC => "CA#", CAG => "CA#", CAU => "CA#",
    CCA => "CC#", CCC => "CC#", CCG => "CC#", CCU => "CC#",
    CGA => "#G#", CGC => "2G#", CGG => "#G#", CGU => "2G#",
    CUA => "#U#", CUC => "2U#", CUG => "#U#", CUU => "2U#",
    GAA => "GA#", GAC => "GA#", GAG => "GA#", GAU => "GA#",
    GCA => "GC#", GCC => "GC#", GCG => "GC#", GCU => "GC#",
    GGA => "GG#", GGC => "GG#", GGG => "GG#", GGU => "GG#",
    GUA => "GU#", GUC => "GU#", GUG => "GU#", GUU => "GU#",
    UAA => "U##", UAC => "UA#", UAG => "U1#", UAU => "UA#",
    UCA => "<:#", UCC => "1!#", UCG => "<:#", UCU => "1!#",
    UGA => "U#!", UGC => "UG#", UGG => "UGG", UGU => "UG#",
    UUA => "#U#", UUC => "UU#", UUG => "#U#", UUU => "UU#",
      );
  for(my $p = 0; $p < (length($seq)); $p += 3){
    if(exists($ambigMap{substr($seq, $p, 3)})){
      substr($seq, $p, 3) = $ambigMap{substr($seq, $p, 3)};
    } else {
      die("Not found: #", substr($seq, $p, 3),"#");
    }
  }
  return($seq);
}

my $debug = 0;
if($ARGV[0] eq "-v"){
  shift(@ARGV);
  if($ARGV[0] =~ /^([0-9]+)$/){
    $debug = $1;
    shift(@ARGV);
  } else {
    $debug = 1;
  }
}

my $tests = <>;

my %shiftMap = (
  "#" => 0,
  "1" => 1, "2" => 2, "!" => -1, "@" => -2,
  "<" => 2, ":" => 1,
    );

for(chomp $tests; $tests; $tests--){
  my $rna = <>; chomp($rna); chomp($rna); $rna = makeAmbig($rna);
  print(STDERR $rna) if ($debug);
  my $rangeCount = <>;
  #printf($rangeCount);
  my $valid = 1; #true
  my %changed = ();
  my $lastBasic = 0;
  my $lastEnd = 0;
  my %rangeEndStarts = ();
  my %rangeStartEnds = ();
  my %origRanges = ();
  my @openClose = (" ") x length($rna);
  ## reduce to unique start/end points for ranges
  for(chomp $rangeCount; $rangeCount; $rangeCount--){
    my $range = <>; chomp $range;
    $origRanges{$range} = 0;
    my ($rs, $re) = split(/\s+/, $range);
    if(($openClose[$rs-1] eq ")") || ($openClose[$rs-1] eq "O")){
      $openClose[$rs-1] = "O";
    } else {
      $openClose[$rs-1] = "(";
    }
    if($openClose[$re-1] ne " "){
      $openClose[$re-1] = "O";
    } else {
      $openClose[$re-1] = ")";
    }
    my $keep = 1; # true
    if(exists($rangeStartEnds{$rs})){
      my $rse = $rangeStartEnds{$rs};
      if($rse > $re){
        delete($rangeStartEnds{$rs});
        delete($rangeEndStarts{$rse});
      } else {
        $keep = 0; # false
      }
    }
    if(exists($rangeEndStarts{$re})){
      my $res = $rangeEndStarts{$re};
      if($res < $rs){
        delete($rangeStartEnds{$res});
        delete($rangeEndStarts{$re});
      } else {
        $keep = 0; # false
      }
    }
    if($keep){
      $rangeStartEnds{$rs} = $re;
      $rangeEndStarts{$re} = $rs;
    }
  }
  print(STDERR join("",@openClose)) if ($debug);
  @openClose = (" ") x length($rna);
  foreach my $rs (sort {$a <=> $b} (keys(%rangeStartEnds))){
    my $re = $rangeStartEnds{$rs};
    my $reb = $rangeEndStarts{$re};
    if(($openClose[$rs-1] eq ")") || ($openClose[$rs-1] eq "O")){
      $openClose[$rs-1] = "O";
    } else {
      $openClose[$rs-1] = "(";
    }
    if($openClose[$re-1] ne " "){
      $openClose[$re-1] = "O";
    } else {
      $openClose[$re-1] = ")";
    }
  }
  print(STDERR join("",@openClose)) if ($debug);
  ## post-start/end filtering
  if($debug > 1){
    print(STDERR " -- post-filtering --") if($debug);
    foreach my $rs (sort {$a <=> $b} (keys(%rangeStartEnds))){
      my $re = $rangeStartEnds{$rs};
      my $reb = $rangeEndStarts{$re};
      print(STDERR "$rs $re / $reb") if($debug);
    }
  }
  my @startList = sort {$a <=> $b} keys(%rangeStartEnds);
  my $lastRs = shift(@startList);
  my $lastRe = $rangeStartEnds{$lastRs};
  $valid = 1; # true
  ## Deal with overlapping regions iteratively, only concentrating on two regions at a time
  ## Assumptions: ranges have unique start and end points
  while($valid && @startList){
    print(STDERR "Start list: ", join(":", @startList)) if($debug);
    my $nextRs = shift(@startList);
    my $nextRe = $rangeStartEnds{$nextRs};
    printf(STDERR "Testing $nextRs / $nextRe pair (%d, %d)\n",
	   $rangeStartEnds{$nextRs}, $rangeEndStarts{$nextRe}) if ($debug);
    ## There are three possible situations
    ## 1) Region 2 is captured within Region 1, in which case Region 1 is removed
    if($nextRe <= $lastRe){
      printf(STDERR "nextRe <= lastRe; ".
	     "Deleting $lastRs / $lastRe pair (%d, %d) in preference for ".
	     "$nextRs / $nextRe pair (%d, %d)\n",
	     $rangeStartEnds{$lastRs}, $rangeEndStarts{$lastRe},
	     $rangeStartEnds{$nextRs}, $rangeEndStarts{$nextRe},
	  ) if ($debug);
      delete($rangeStartEnds{$lastRs});
      delete($rangeEndStarts{$lastRe});
      $lastRs = $nextRs;
      $lastRe = $nextRe;
    } elsif($lastRe > $nextRs) {
      ## 2) Region 1 overlaps Region 2. This is broken into five sub-cases
      my $intersectRs = $nextRs;
      my $intersectRe = $lastRe;
      my $intSst = substr($rna, $intersectRs-1, $intersectRe-$intersectRs+1);
      if($intSst =~ tr/#12!@//){
	## 2.a) the intersection contains a changeable portion of no
	##      greater than two changes, in which case the intersect is
	##      kept
	print(STDERR "Intersection $intersectRs..$intersectRe is 2-changeable (or better): $intSst") if ($debug);
	printf(STDERR "Deleting $lastRs / $lastRe pair (%d, %d)\n",
	       $rangeStartEnds{$lastRs}, $rangeEndStarts{$lastRe}) if ($debug);
	delete($rangeStartEnds{$lastRs});
	delete($rangeEndStarts{$lastRe});
	printf(STDERR "Deleting $nextRs / $nextRe pair (%d, %d)\n",
	       $rangeStartEnds{$nextRs}, $rangeEndStarts{$nextRe}) if ($debug);
	delete($rangeStartEnds{$nextRs});
	delete($rangeEndStarts{$nextRe});
	## add in the new region
	print(STDERR "Adding intersect range $intersectRs..$intersectRe") 
	    if ($debug);
	$rangeStartEnds{$intersectRs} = $intersectRe;
	$rangeEndStarts{$intersectRe} = $intersectRs;
	$lastRs = $intersectRs;
	$lastRe = $intersectRe;
      } else {
	my $leftSst = substr($rna, $lastRs-1, $intersectRs-1-$lastRs+1);
	my $rightSst = substr($rna, $intersectRe+1-1, $nextRe-$intersectRe+1-1+1);
	if((index($leftSst,"#") + 1) && (index($rightSst,"#") + 1)){
	  ## 2.b) non-intersecting regions both contain 1-changeable
	  ##      portions; they are chosen in preference.
	  print(STDERR "Non-intersects outside $intersectRs..$intersectRe ".
		"are 1-changeable: [$leftSst] ; [$rightSst]") if ($debug);
	  printf(STDERR "Deleting $lastRs / $lastRe pair (%d, %d)\n",
		 $rangeStartEnds{$lastRs}, $rangeEndStarts{$lastRe}) if ($debug);
	  delete($rangeStartEnds{$lastRs});
	  delete($rangeEndStarts{$lastRe});
	  printf(STDERR "Deleting $nextRs / $nextRe pair (%d, %d)\n",
		 $rangeStartEnds{$nextRs}, $rangeEndStarts{$nextRe}) if ($debug);
	  delete($rangeStartEnds{$nextRs});
	  delete($rangeEndStarts{$nextRe});
	  ## add in the new regions
	  print(STDERR "Adding left range $lastRs..".
		($intersectRs-1)) if ($debug);
	  $rangeStartEnds{$lastRs} = $intersectRs-1;
	  $rangeEndStarts{$intersectRs-1} = $lastRs;
	  if(exists($rangeStartEnds{$intersectRe+1})){
	    ## replacement already exists that starts in the same place, choose
	    ## the shortest
	    my $eRe = $rangeStartEnds{$intersectRe+1};
	    print(STDERR "Existing value exists for start at ".
		  ($intersectRe+1)."..$eRe") if ($debug);
	    if($eRe < $nextRe){
	      $nextRe = $eRe;
	    }
	  } else {
	    ## this situation ends a little differently; the right hand
	    ## side is pushed back on the stack, and the stack is re-sorted
	    unshift(@startList, $intersectRe+1);
	    @startList = sort {$a <=> $b} @startList;
	  }
	  print(STDERR "Adding right range ".
		($intersectRe+1)."..$nextRe") if ($debug);
	  $rangeStartEnds{$intersectRe+1} = $nextRe;
	  $rangeEndStarts{$nextRe} = $intersectRe+1;
	  $lastRe = $intersectRs - 1;
	} elsif($intSst =~ tr/<://){
	  ## 2.c) intersecting region contains a 3-changeable portion;
	  ##      it is kept
	  print(STDERR "Intersection $intersectRs..$intersectRe ".
		"is 3-changeable: $intSst") if ($debug);
	  printf(STDERR "Deleting $lastRs / $lastRe pair (%d, %d)\n",
		 $rangeStartEnds{$lastRs}, $rangeEndStarts{$lastRe}) if ($debug);
	  delete($rangeStartEnds{$lastRs});
	  delete($rangeEndStarts{$lastRe});
	  printf(STDERR "Deleting $nextRs / $nextRe pair (%d, %d)\n",
		 $rangeStartEnds{$nextRs}, $rangeEndStarts{$nextRe}) if ($debug);
	  delete($rangeStartEnds{$nextRs});
	  delete($rangeEndStarts{$nextRe});
	  ## add in the new region
	  print(STDERR "Adding range $intersectRs..$intersectRe") if ($debug);
	  $rangeStartEnds{$intersectRs} = $intersectRe;
	  $rangeEndStarts{$intersectRe} = $intersectRs;
	  $lastRs = $intersectRs;
	  $lastRe = $intersectRe;
	} elsif(($leftSst =~ tr/#12!@<://) && ($rightSst =~ tr/#12!@<://)){
	  ## 2.d) non-intersecting regions are both changeable [but
	  ##      not 1-changeable]
	  print(STDERR "Non-intersects outside $intersectRs..$intersectRe ".
		"are 2/3-changeable: [$leftSst] ; [$rightSst]") if ($debug);
	  printf(STDERR "Deleting $lastRs / $lastRe pair (%d, %d)\n",
		 $rangeStartEnds{$lastRs}, $rangeEndStarts{$lastRe}) if ($debug);
	  delete($rangeStartEnds{$lastRs});
	  delete($rangeEndStarts{$lastRe});
	  printf(STDERR "Deleting $nextRs / $nextRe pair (%d, %d)\n",
		 $rangeStartEnds{$nextRs}, $rangeEndStarts{$nextRe}) if ($debug);
	  delete($rangeStartEnds{$nextRs});
	  delete($rangeEndStarts{$nextRe});
	  ## add in the new regions
	  print(STDERR "Adding range $lastRs..".($intersectRs-1)) if ($debug);
	  $rangeStartEnds{$lastRs} = $intersectRs-1;
	  $rangeEndStarts{$intersectRs-1} = $lastRs;
	  print(STDERR "Adding range ".($intersectRe+1)."..$nextRe") if ($debug);
	  $rangeStartEnds{$intersectRe+1} = $nextRe;
	  $rangeEndStarts{$nextRe} = $intersectRe+1;
	  ## this situation ends a little differently; the right hand
	  ## side is pushed back on the stack, and the stack is re-sorted
	  unshift(@startList, $intersectRe+1);
	  @startList = sort {$a <=> $b} @startList;
	  $lastRe = $intersectRs - 1;
	} else {
	  ## 2.e) there are no changeable regions anywhere; a solution is impossible
	  print(STDERR "No changeable regions for $lastRs..$lastRe and $nextRs..$nextRe") if($debug);
	  print(STDERR substr($rna, $lastRs-1, $nextRe-$lastRs+1)) if($debug);
	  $valid = 0; # false
	}
      }
    } else {
      ## 3) regions don't overlap; move along, nothing to see here
      $lastRs = $nextRs;
      $lastRe = $nextRe;
    }
  }
  if(!$valid){
    print(STDERR "no possible solutions; clearing both range arrays");
    %rangeEndStarts = ();
    %rangeStartEnds = ();
  }
  @openClose = ("#") x length($rna);
  foreach my $rs (sort {$a <=> $b} (keys(%rangeStartEnds))){
    my $re = $rangeStartEnds{$rs};
    my $reb = $rangeEndStarts{$re};
    if(($openClose[$rs-1] eq ")") || ($openClose[$rs-1] eq "O")){
      $openClose[$rs-1] = "O";
    } else {
      $openClose[$rs-1] = "(";
    }
    if(($openClose[$re-1] eq "(") || ($openClose[$re-1] eq "O")){
      $openClose[$re-1] = "O";
    } else {
      $openClose[$re-1] = ")";
    }
  }
  print(STDERR join("",@openClose)) if ($debug);
  if($debug > 2){
    print(STDERR " -- post-overlap --") if($debug);
    foreach my $rs (sort {$a <=> $b} (keys(%rangeStartEnds))){
      my $re = $rangeStartEnds{$rs};
      print(STDERR "$rs $re") if($debug);
    }
  }
  foreach my $rs (sort {$a <=> $b} (keys(%rangeStartEnds))){
    my $re = $rangeStartEnds{$rs};
    my $sst = substr($rna, $rs-1, $re-$rs+1);
    if(($changed{$rs}) || (($rs != $re) && $changed{$rs+1})){
      ## already changed, nothing to do
    } elsif($sst =~ tr/#//){
      $lastBasic = $rs + index($sst, "#");
      $changed{$lastBasic} = 1;
      print(STDERR "Added basic at $lastBasic") if ($debug > 3);
    } elsif($sst =~ tr/12!@//) {
      my $ced = substr($sst, length($sst)-1, 1);
      my $cst = substr($sst, 0, 1);
      if($ced eq '1'){
        $changed{$re} = 1;
        $changed{$re + 1} = 1;
      } elsif($ced eq '2'){
        $changed{$re} = 1;
        $changed{$re + 2} = 1;
      } elsif($cst eq '!'){
        $changed{$rs} = 1;
        $changed{$rs - 1} = 1;
        if(($lastEnd == ($rs-1)) && $lastBasic && ($lastBasic != ($rs - 1))){
          print(STDERR "Deleting basic at $lastBasic") if ($debug > 3);
          delete($changed{$lastBasic});
        } else {
          print(STDERR "No basic to delete (was $lastBasic, last end was $lastEnd, current start is $rs)") if ($debug > 3);
        }
      } elsif($cst eq '@'){
        $changed{$rs} = 1;
        $changed{$rs - 2} = 1;
      } elsif($sst =~ /([12!@])/) {
        my $symbol = $1;
        my $pos = index($sst, $1);
        $changed{$rs + $pos} = 1;
        $changed{$rs + $pos + $shiftMap{$symbol}} = 1;
      }
      $lastBasic = 0;
    } elsif($sst =~ tr/<://) {
      my $ced = substr($sst, length($sst)-1, 1);
      my $cst = substr($sst, 0, 1);
      if($ced eq '<'){
        $changed{$re} = 1;
        $changed{$re + 1} = 1;
        $changed{$re + 2} = 1;
      } elsif($ced eq ':'){
        $changed{$re - 1} = 1;
        $changed{$re + 0} = 1;
        $changed{$re + 1} = 1;
        if(($lastEnd == ($re-1)) && $lastBasic && ($lastBasic != ($re - 1))){
          delete($changed{$lastBasic});
        }
      } elsif($cst eq ':'){
        $changed{$rs - 1} = 1;
        $changed{$rs + 0} = 1;
        $changed{$rs + 1} = 1;
        if(($lastEnd == ($rs-1)) && $lastBasic && ($lastBasic != ($rs - 1))){
          delete($changed{$lastBasic});
        }
      } else {
        my $pos = index($sst, "<");
        $changed{$rs + $pos + 0} = 1;
        $changed{$rs + $pos + 1} = 1;
        $changed{$rs + $pos + 2} = 1;
      }
      $lastBasic = 0;
    } else {
      $valid = 0;
    }
    print(STDERR "$rs..$re $valid ".substr($rna, $rs-1, $re-$rs+1).
	  " $lastBasic") if ($debug);
    $lastEnd = $re;
  }
  # my @changedPoss = sort {$a <=> $b} (keys(%changed));
  # print(STDERR "[".join(";", @changedPoss)."]") if ($debug);
  # my %captured = ();
  # foreach my $range (keys(%origRanges)){
  #   my ($rs, $re) = split(/\s+/, $range);
  #   my @capPoss = grep {($_ >= $rs) && ($_ <= $re)} (@changedPoss);
  #   $origRanges{$range} = scalar(@capPoss);
  #   if(scalar(@capPoss) == 1){
  #     $captured{$capPoss[0]} = 1;
  #   }
  # }
  # foreach my $pos (sort {$a <=> $b} keys(%captured)){
  #   my $sym = substr($rna, $pos-1, 1);
  #   my $offset = $shiftMap{$sym};
  #   if(!$captured{$pos+$offset}){
  #     print(STDERR "$pos => ".($pos + $offset)) if($debug);
  #     $captured{$pos+$offset} = 1;
  #     if($sym eq ":"){
  # 	$captured{$pos-$offset} = 1;
  #     }
  #   }
  # }
  # if($debug){
  #   if(scalar(keys(%captured)) != scalar(@changedPoss)){
  #     printf(STDERR "Only %d of %d locations are necessary for all ranges\n",
  # 	     scalar(keys(%captured)), scalar(@changedPoss));
  #     print(STDERR "Non-uniquely captured locations:") ;
  #     print(STDERR "  ".join(", ", (grep {!$captured{$_}} @changedPoss)));
  #   }
  # }
  # print($valid ? scalar(keys(%captured)) : -1);
  print($valid ? scalar(keys(%changed)) : -1); ## This seems to be more correct
}
