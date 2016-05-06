  ###################################################
  #Marten Boetzer 23-10-2010                        #
  #SSPACE perl subscript PairingAndScaffolding.pl   #
  #This script;                                     #
  #  -reads the contig sequences in a hash          #
  #  -stores Bowtie output in a hash                #
  #  -pairs the contigs                             #
  #  -generates scaffolds                           #
  ###################################################


  #THIS VERSION OF SCAFFOLDING FIRST ORDERS THE CONTIGS BASED ON THE NUMBER OF INGOING LINKS AND STARTS AT LOWEST LEVEL. AFTER ALL THESE CONTIGS ARE SCAFFOLDED, INGOING LINKS ARE RECALCULATED OF REMAINING CONTIGS, ITERATIVELY. 
  #ALSO, EACH CONTIG IS REPRESENTED ONCE IN THE SCAFFOLDS. 
  #METHOD OF SCAFFOLDING IS; IF MORE THAN ONE LINK, CHECK IF THOSE LINKS HAVE CONNECTION WITH EACH OTHER. IF SO, COMBINE THEM IN THE SCAFFOLD. IF NOT, ESTIMATE RATIO AND ONLY ALLOW EXTENSION OF SCAFFOLD IF IT'S BELOW THE RATIO THRESHOLD GIVEN BY THE USER.
  #FUTURE: INCLUDE NUMBER OF REPEATS THAT ARE POSSIBLY PRESENT
  use strict;
  use Storable;
  use File::Path;
  use File::Basename;

  my $seplines = ("-" x 60)."\n";

  my $contig = $ARGV[0];
  my $base_name = $ARGV[1];
  my $issues = $ARGV[2];
  my $distribution = $ARGV[3];
  my $verbose = $ARGV[4];
  my $library = $ARGV[5];
  my $insert_size = $ARGV[6];
  my $min_allowed = $ARGV[7];
  my $scaffold = $ARGV[8];
  my $min_links = $ARGV[9];
  my $max_link_ratio = $ARGV[10];
  my $outdir = $ARGV[11];

  my ($low_iz, $up_iz) = ($insert_size + $min_allowed, $insert_size - $min_allowed);
  my $bowtiefile = "$outdir/bowtieoutput/" . $base_name . ".$library.mapped";
  my $log = $outdir."/".$base_name . ".logfile.txt";
  my $summaryfile = $outdir."/".$base_name.".summaryfile.txt";

  open (LOG, ">>$log") || die "Can't write to $log -- fatal\n";
  open (SUMFILE, ">>$summaryfile") || die "Can't open $summaryfile -- fatal\n";

  my($contigstored, $tig_length) = &readFileContigHash($contig);
  my ($track_all, $fileread) = &readBowtie($bowtiefile, $tig_length);

#-------------------------------------------------READ CONTIGS INTO HASH AND STORE THEIR LENGTH. NEXT; PAIR THE CONTIGS
  my $pair = &pairContigs($fileread, $track_all, $tig_length, $issues, $distribution, $verbose, $library);
  ($track_all) =  (''); undef $track_all;

#-------------------------------------------------BUILDING SCAFFOLDS
  &buildScaffolds($pair, $tig_length, $verbose, $scaffold, $library);
  ($pair, $tig_length) = ('',''); undef $pair; undef $tig_length;

  close SUMFILE;
  close LOG;

#-------------------------------------------------

#READ THE CONTIG TO A HASH AND STORE THIS HASH
sub readFileContigHash{
  my ($file) = @_;

  &printMessage("\n=>".getDate().": Reading contig file($contig)\n");
    
  my ($contigs, $tig_length);
  my ($countContig, $seq, $prevhead) = (0, "", '');
  open(IN,$file) || die "Can't open $file -- fatal\n";
  while(<IN>){
     my $line = $_;
     $line =~ s/\r\n/\n/;
     chomp $line; 
     $seq.= uc($line) if(eof(IN));
     if (/\>(\S+)/ || eof(IN)){
       my $head=$1;
       if($prevhead ne ''){
         CounterPrint(++$countContig);
         $tig_length->{$countContig} = length($seq);
         $contigs->{$countContig}{'name'} = $prevhead;
         $contigs->{$countContig}{'seq'} = $seq;
       }
       $prevhead = $head;
       $seq='';
     }else{
        $seq.=$line;
     }
  }
  CounterPrint("                ");
  &FlushFiles();
  my $contigstore = "$outdir/tmp.$base_name/contigs.stored";
  store \%$contigs, "$contigstore";
  undef $contigs;
  return ($contigstore, $tig_length);
}

sub determineRepeats{
  my ($tig_length, $pair, $repeathash, $step) = @_;
 foreach my $tig (sort {$tig_length->{$b}<=>$tig_length->{$a}} keys %$tig_length){
    if(defined $repeathash->{$tig} || $step == 0){
        my $rtig = "r" . $tig;
        my $ftig = "f" . $tig;
        my $list = $pair->{$rtig};
        my ($countlinksforward, $countlinksbackward)= (0,0);
        my (@forward, @backward);
        my ($seenrepeatforw, $seenrepeatback) = ("", "");
        foreach my $match (sort {$list->{$b}{'links'}<=>$list->{$a}{'links'}} keys %$list){
          push @backward, $match if($list->{$match}{'links'} >= $min_links);
        }
        $countlinksbackward = $#backward + 1;
        my ($seenforward, $seenbackward);
        foreach my $backmatches (@backward){
          if(!defined $seenbackward->{$backmatches}){
            my $backlist = $pair->{$backmatches};
            my $duplicate = 0;
            foreach my $backmatches2 (@backward){
              if($backlist->{$backmatches2}{'links'} > 0 && $backmatches ne $backmatches2){
               $seenbackward->{$backmatches2}++;
               $countlinksbackward--;
              }
            }
          }
        }
        my $list2 = $pair->{$ftig};
        foreach my $match2 (sort {$list2->{$b}{'links'}<=>$list2->{$a}{'links'}} keys %$list2){
          push @forward, $match2 if($list2->{$match2}{'links'} >= $min_links);
        }

        foreach my $formatches (@forward){
          if(!defined $seenforward->{$formatches}){
            my $forlist = $pair->{$formatches};
            my $duplicate = 0;
            foreach my $formatches2 (@forward){
              if($forlist->{$formatches2}{'links'} > 0 && $formatches ne $formatches2){
               $seenforward->{$formatches2}++;
               $duplicate++;
              }
            }
            $countlinksforward++ if($duplicate == 0);
          }
        }
  
        my $repeatnum = 0;
        if($countlinksbackward > $countlinksforward){
           $repeatnum = $countlinksbackward;
        }else{
           $repeatnum = $countlinksforward;
        }
        if ($repeatnum > 2){
          if(!defined $repeathash->{$tig}){
            $repeathash->{$tig}->{'size'} = $tig_length->{$tig};
            $repeathash->{$tig}->{'links'} = $repeatnum;
          }
          if($countlinksbackward == 1 && $step == 1){
            foreach my $backmatches (@backward){
              my $tnum = $1 if($backmatches=~/[fr](\d+)/);
              if(!defined $seenbackward->{$backmatches} && !defined $repeathash->{$tnum}){
                 $repeathash->{$tnum}->{'size'} = $tig_length->{$tnum};
                 $repeathash->{$tnum}->{'links'} = $repeatnum;
              }
            }
          }
          if($countlinksforward == 1 && $step == 1){
            foreach my $formatches (@forward){
              my $tnum = $1 if($formatches=~/[fr](\d+)/);
              if(!defined $seenforward->{$formatches} && !defined $repeathash->{$tnum}){
               $repeathash->{$tnum}->{'size'} = $tig_length->{$tnum};
                $repeathash->{$tnum}->{'links'} = $repeatnum;
              }
            }
          }
        }
      }
   }
  return $repeathash;
}

sub buildScaffolds{

   my ($pair, $tig_length, $verbose, $scaffold, $lib) = @_;
   &printMessage("\n=>".getDate().": Building scaffolds file ($scaffold)\n");

   open (SC, ">$scaffold") || die "Can't write to $scaffold -- fatal\n";
   my ($sc_ct, $keyrep) = (0,0);
   my ($repeathash, $seen_start);

   $repeathash = determineRepeats($tig_length, $pair, $repeathash, 0);
   my $keyrep = 0;
   my $keyrepprev = keys (%$repeathash);

   while($keyrep ne $keyrepprev){
     $keyrep = keys (%$repeathash);
     $repeathash = determineRepeats($tig_length, $pair, $repeathash, 1);
     $keyrepprev = keys (%$repeathash);
   }
   if(keys (%$repeathash) > 0){
     my $totsize=0;
     open(REPEATS, ">$outdir/intermediate_results/$base_name.".$lib.".repeats") || die "Can't open $outdir/intermediate_results/$base_name.$lib.repeats -- fatal\n";
     foreach my $sizetig (keys %$repeathash){
       $totsize += ($repeathash->{$sizetig}->{'size'} * $repeathash->{$sizetig}->{'links'});
       print REPEATS "contig $sizetig is repeated $repeathash->{$sizetig}->{'links'} times. Size = $repeathash->{$sizetig}->{'size'}\n";
     }

     print SUMFILE "REPEATS: \n".$seplines;
     print SUMFILE "\tNumber of repeats = ".keys (%$repeathash)."\n";
     print SUMFILE "\tTotal size of repeats = $totsize\n$seplines.\n";

     close REPEATS;
   }

   SEED:
   foreach my $tig (sort {$tig_length->{$b}<=>$tig_length->{$a}} keys %$tig_length){
      my $ftig = "f" . $tig;
      my $rtig = "r" . $tig;

      if(! defined $seen_start->{$tig}){##should prevent re-using a contig as seed if it's already been incorporated into a scaffold
         CounterPrint(++$sc_ct);
         my $chainleft = "";
         my $ori_chainright = $ftig . "Z" . $tig_length->{$tig};
         my $chainright = $ori_chainright;
         my $total = $tig_length->{$tig};
         ($total, $chainright, $seen_start) = &computeLayout("R", $chainright, $ftig, $pair, $tig_length, $total, $seen_start, $tig);
         ($total, $chainleft, $seen_start) = &computeLayout("L", $chainleft, $rtig, $pair, $tig_length, $total, $seen_start, $tig);

         delete $pair->{$ftig};
         delete $pair->{$rtig};
         delete $tig_length->{$tig};
         $seen_start->{$tig}++;
         my $scaffold = $chainleft . $chainright;
         print SC "scaffold" . $sc_ct . ",$total,$scaffold\n";
      }
   }
   CounterPrint("                ");
   close SC;
   &FlushFiles();
}


sub computeLayout{
   my ($ext, $chain, $tig, $pair, $tig_length, $total, $seen_start, $orig_tig_number) = @_;
   my $orig_tig = $tig;
   my $extension = 1;
   EXTENSION:
   while($extension){
      my $tnum = $1 if($tig=~/[fr](\d+)/);
      my $tnumf = "f" . $tnum;
      my $tnumr = "r" . $tnum;
      my $ratio = 0.00;
      if(!defined $seen_start->{$tnum} || $tnumf eq $orig_tig || $tnumr eq $orig_tig){
         $seen_start->{$tnum}++;
         my $list = $pair->{$tig};
         my $matchhash;
         my ($match1,$link1,$gaps1,$match2,$link2,$gaps2,$cntloop, $countmatches)=("",0,0,"",0,0,0,0);

         LINK:
         foreach my $match (sort {$list->{$b}{'links'}<=>$list->{$a}{'links'}} keys %$list){
            my $matchnum = $1 if($match=~/[fr](\d+)/);
            if($list->{$match}{'links'} >= $min_links && !defined $seen_start->{$matchnum}){
              $matchhash->{$match}{'links'} = $list->{$match}{'links'};
              $matchhash->{$match}{'gaps'} = $list->{$match}{'gaps'};
              $matchhash->{$match}{'ratio'} = $list->{$match}{'gaps'}/$list->{$match}{'links'};
              $countmatches++;
            }else{
              last LINK;
            }
         }
         my $foundlinks = 0;
         if($countmatches > 1){
           my @arraymatch;
           foreach my $ratiosort (sort {$matchhash->{$a}{'ratio'}<=>$matchhash->{$b}{'ratio'}} keys %$matchhash){
             push @arraymatch, $ratiosort;
           }
           my $nummatch = $#arraymatch;
           for(my $i=0; $i <= $nummatch && $foundlinks < 1; $i++){
             my $listmatch = $pair->{$arraymatch[$i]};
              for(my $j=$i+1; $j <= $nummatch && $foundlinks < 1; $j++){
                 my $linkmatch = $listmatch->{$arraymatch[$j]}{'links'};
                 $foundlinks = 1 if(!($linkmatch >= $min_links));
              }
           }
           my $tignum = $1 if($arraymatch[$nummatch]=~/[fr](\d+)/);
           $countmatches=0 if(!$foundlinks && defined $seen_start->{$tignum});
         }
         if($foundlinks && $countmatches > 1){
             my @linkmatch;
             foreach my $linksort (sort {$matchhash->{$b}{'links'}<=>$matchhash->{$a}{'links'}} keys %$matchhash){
               push @linkmatch, $linksort;
             }
             my $link1 = $matchhash->{$linkmatch[1]}{'links'};
             my $link2 = $matchhash->{$linkmatch[0]}{'links'};
             $ratio = $link1 / $link2;        ## relative ratio of the two most abundant contig pairs
             if ($ratio =~ /(\d+\.\d{2})/){$ratio = $1;}
             if($ratio <= $max_link_ratio){
                 foreach my $mat (keys %$matchhash){
                   delete $matchhash->{$mat} if($mat ne $linkmatch[0]);
                 }
              $foundlinks = 0;
              $countmatches = 1;

             }
         }
         if((!$foundlinks) && $countmatches > 0){
           my $nummatch =0;
           my @chainlist;
           my @tiglist;
           foreach my $incl_matches (sort {$matchhash->{$a}{'ratio'}<=>$matchhash->{$b}{'ratio'}} keys %$matchhash){
             if($tig ne $incl_matches){
               $nummatch++;
               my $listmatch = $pair->{$tig};
               my $tempnum = $1 if($incl_matches =~ /[fr](\d+)/);
               my $link2 = $listmatch->{$incl_matches}{'links'};
               my $mean2 = $listmatch->{$incl_matches}{'gaps'}/$link2;

               $seen_start->{$tempnum}++if($nummatch < $countmatches);

               ($chain, $total, $tig) = &getChain($chain, $ext, $link2, $mean2, $incl_matches, $tempnum, $ratio, $tig_length, $total);
               delete $tig_length->{$tempnum};
             }
           }
           $extension = 1;

         }else{
           $extension = 0;
           last EXTENSION;
         }
      }else{
           $extension = 0;
           last EXTENSION;
      }
   }
   return $total, $chain, $seen_start;
}


sub getChain{
  my ($chain, $ext, $link, $mean, $match, $tempnum, $ratio, $tig_length, $total) = @_;
  my $tig = $match;
  if($ext eq "R"){
               $chain .= "k" . $link . "a" . $ratio . "m" . int($mean) . "_" . $match . "z" . $tig_length->{$tempnum};
  }else{
    my $temp_match = "";
    if($match =~ /^r(\d+)/){$temp_match = "f" . $1;}else{$temp_match = "r". $1;}
     $chain = $temp_match . "z" . $tig_length->{$tempnum} . "k" . $link . "a" . $ratio . "m" . int($mean) . "_" . $chain;
 }

  $total += $tig_length->{$tempnum};
  return ($chain, $total, $tig);
}


###GET THE DISTANCE BETWEEN TWO PAIRED READS
sub getDistance{

   my ($insert_size, $length_i, $start_i, $start_j) = @_;

   # L  ------  --------- R
   # i    ->        <-    j
   #      ....  ......    insert_span
   #      ============    insert_size

   my $insert_span = ($length_i - $start_i) + $start_j;
   my $gap_or_overlap = $insert_size - $insert_span;

   return $gap_or_overlap;
}

###PAIR CONTIGS BASED ON THE NUMBER OF LINKS
sub pairContigs{

   my ($readsfile,$track,$tig_length,$issues,$distribution,$verbose, $library) = @_;

   &printMessage("\n=>".getDate().": Reading paired reads and pairing contigs.\n");
   &FlushFiles();
   my ($step,$ct_illogical, $ct_ok_contig, $ct_ok_pairs, $ct_problem_pairs, $ct_iz_issues, $ct_single, $ct_both, $counter)= (100000,0,0,0,0,0,0,0,0 );
   my ($pair,$err,$track_insert);
   my %count = ();

   open(PET, ">$issues") || die "Can't open $issues for writing -- fatal\n";
   open(IN,$readsfile) || die "Can't open $readsfile -- fatal\n";
   while(<IN>){
     chomp;
     if(/^([^\>]*)$/i){
       my $sdna = $1;
       if($sdna =~/([ACGT]*)\:([ACGT]*)/i){
         ++$counter;
         if($counter == $step){
            CounterPrint($counter);
            $step = $step + 100000;
         }
         my $read_a = $1;
         my $read_b = $2;;
         if(defined $track->{$read_a} && defined $track->{$read_b}){
          my $combined = "$read_a:$read_b";
          my $revcombined = reverseComplement($combined);
          if(!defined $count{$combined} && !defined $count{$revcombined}){
            $count{$combined}++;

            ### both pairs assembled
             $ct_both++;
          
             my ($tig_a, $A_start, $A_end) = split(/\|/, $track->{$read_a});
             my ($tig_b, $B_start, $B_end) = split(/\|/, $track->{$read_b});
             
             my $ftig_a = "f" . $tig_a;
             my $ftig_b = "f" . $tig_b;
  
             my $rtig_a = "r" . $tig_a;
             my $rtig_b = "r" . $tig_b;
    
             my $A_length = $tig_length->{$tig_a};
             my $B_length = $tig_length->{$tig_b};
  
             if ($tig_a != $tig_b){####paired reads located on <> contigs
    
                ####Determine most likely possibility
                if ($A_start < $A_end){
    
                   if ($B_end < $B_start){####-> <- :::  A-> <-B  /  rB -> <- rA
                       my $d = &getDistance($insert_size, $A_length, $A_start, $B_start);
                       print "A-> <-B  WITH $tig_a -> <- $tig_b GAP $d A=$A_length ($A_start-$A_end) B=$B_length ($B_start-$B_end) Alen, Astart,Bstart\n" if($verbose);
                       if($d >= $min_allowed){
                          $pair->{$ftig_a}{$ftig_b}{'links'}++;
                          $pair->{$ftig_a}{$ftig_b}{'gaps'} += $d;
                          $pair->{$rtig_b}{$rtig_a}{'links'}++;
                          $pair->{$rtig_b}{$rtig_a}{'gaps'} += $d;
                          $ct_ok_pairs++;
                       }else{
                          my $err_pair = $ftig_a . "-". $ftig_b;
                          $err->{$err_pair}{'links'}++;
                          $err->{$err_pair}{'gaps'} += $d;
                          $ct_problem_pairs++;
                          print PET "Pairs unsatisfied in distance within a contig pair.  A-> <-B  WITH tig#$tig_a -> $d <- tig#$tig_b, A=$A_length nt (start:$A_start, end:$A_end) B=$B_length nt (start:$B_start, end:$B_end) CALCULATED DISTANCE APART: $d < $min_allowed\n";
                       }
                    }else{#### -> -> ::: A-> <-rB  / B-> <-rA 
                       my $rB_start = $B_length - $B_start;
                       my $d = &getDistance($insert_size, $A_length, $A_start, $rB_start);
                       print "A-> <-rB  WITH $tig_a -> <- r.$tig_b GAP $d A=$A_length ($A_start-$A_end) B=$B_length ($B_start-$B_end) Alen,Astart,rBstart\n"if($verbose);
                       if($d >= $min_allowed){
                          $pair->{$ftig_a}{$rtig_b}{'links'}++;
                          $pair->{$ftig_a}{$rtig_b}{'gaps'} += $d;
                          $pair->{$ftig_b}{$rtig_a}{'links'}++;
                          $pair->{$ftig_b}{$rtig_a}{'gaps'} += $d;
                          $ct_ok_pairs++;
                       }else{
                          my $err_pair = $ftig_a . "-". $rtig_b;
                          $err->{$err_pair}{'links'}++;
                          $err->{$err_pair}{'gaps'} += $d;
                          $ct_problem_pairs++;
                          print PET "Pairs unsatisfied in distance within a contig pair.  A-> <-rB  WITH tig#$tig_a -> $d <- tig#r.$tig_b, A=$A_length  nt (start:$A_start, end:$A_end) B=$B_length nt (start:$B_start, end:$B_end) CALCULATED DISTANCE APART: $d < $min_allowed\n";
                       }
                    }
                }else{
    
                   if ($B_end > $B_start){####<-  -> ::: B-> <-A / rA -> <- rB
                      my $d = &getDistance($insert_size, $B_length, $B_start, $A_start);
                      print "B-> <-A  WITH $tig_b -> <- $tig_a GAP $d A=$A_length ($A_start-$A_end) B=$B_length ($B_start-$B_end) Blen,Bstart,Astart\n" if($verbose);
                      if($d >= $min_allowed){
                         $pair->{$ftig_b}{$ftig_a}{'links'}++;
                         $pair->{$ftig_b}{$ftig_a}{'gaps'} += $d;
                         $pair->{$rtig_a}{$rtig_b}{'links'}++;
                         $pair->{$rtig_a}{$rtig_b}{'gaps'} += $d;
                         $ct_ok_pairs++;
                      }else{
                         my $err_pair = $ftig_b . "-". $ftig_a;
                         $err->{$err_pair}{'links'}++;
                         $err->{$err_pair}{'gaps'} += $d;
                         $ct_problem_pairs++;
                         print PET "Pairs unsatisfied in distance within a contig pair.  B-> <-A  WITH tig#$tig_b -> $d <- tig#$tig_a, B=$B_length nt (start:$B_start, end:$B_end) A=$A_length nt (start:$A_start, end:$A_end) CALCULATED DISTANCE APART: $d < $min_allowed\n";
                      }
                   }else{                          ####<- <-  :::  rB-> <-A / rA-> <-B
                      my $rB_start = $B_length - $B_start;
                      my $d = &getDistance($insert_size, $B_length, $rB_start, $A_start);
                      print "rB-> <-A WITH r.$tig_b -> <- $tig_a GAP $d A=$A_length ($A_start-$A_end) B=$B_length ($B_start-$B_end) Blen,rBstart,Astart\n" if($verbose);
                      if($d >= $min_allowed){
                         $pair->{$rtig_b}{$ftig_a}{'links'}++;
                         $pair->{$rtig_b}{$ftig_a}{'gaps'} += $d;
                         $pair->{$rtig_a}{$ftig_b}{'links'}++;
                         $pair->{$rtig_a}{$ftig_b}{'gaps'} += $d;
                         $ct_ok_pairs++;
                      }else{
                         my $err_pair = $rtig_b . "-". $ftig_a;
                         $err->{$err_pair}{'links'}++;
                         $err->{$err_pair}{'gaps'} += $d;
                         $ct_problem_pairs++;
                         print PET "Pairs unsatisfied in distance within a contig pair.  rB-> <-A WITH tig#r.$tig_b -> $d <- tig#$tig_a, B=$B_length nt (start:$B_start, end:$B_end) A=$A_length nt (start:$A_start, end:$A_end) CALCULATED DISTANCE APART: $d < $min_allowed\n";
                      }
                   }
                }
             }else{###Clone, paired reads located on the same contig -- could be used to investigate misassemblies
               
                print "Pair ($read_a and $read_b) located on same contig $tig_a ($A_length nt)\n" if ($verbose);
                my $pet_size = 0;
    
                if ($A_start > $B_start && ($B_start < $B_end) && ($A_start > $A_end)){    # B --> <-- A
                   $pet_size = $A_start - $B_start;
                   $track_insert->{$pet_size}++;
                   if($pet_size >= $low_iz && $pet_size <= $up_iz){
                      $ct_ok_contig++;
                   }else{
                      print PET "Pairs unsatisfied in distance within a contig.  Pair ($read_a - $read_b) on contig $tig_a ($A_length nt) Astart:$A_start Aend:$A_end Bstart:$B_start Bend:$B_end CALCULATED DISTANCE APART: $pet_size\n";
                      $ct_iz_issues++;
                   }
                }elsif($B_start > $A_start && ($B_start > $B_end) && ($A_start < $A_end)){ # A --> <-- B
                   $pet_size = $B_start - $A_start;
                   $track_insert->{$pet_size}++;
                   if($pet_size >= $low_iz && $pet_size <= $up_iz){
                      $ct_ok_contig++;
                   }else{
                      print PET "Pairs unsatisfied in distance within a contig.  Pair ($read_a - $read_b) on contig $tig_a ($A_length nt) Astart:$A_start Aend:$A_end Bstart:$B_start Bend:$B_end CALCULATED DISTANCE APART: $pet_size\n";
                      $ct_iz_issues++;
                   }
                }else{
                   $ct_illogical++;
                   print PET "Pairs unsatisfied in pairing logic within a contig.  Pair ($read_a - $read_b) on contig $tig_a ($A_length nt) Astart:$A_start Aend:$A_end Bstart:$B_start Bend:$B_end\n";
                }
             }
           }
        }else{###both pairs assembled
           $ct_single++;
        }
       }
     }
   }###each template
   CounterPrint("                ");

   my $read_number_message = "\tNumber of pairs found with pairing contigs / total pairs = $ct_both / $counter\n";
   printf SUMFILE $read_number_message.$seplines."\n";

   ### summary of contig pair issues
   print PET "------------- Putative issues with contig pairing - Summary  ----------------\n";
   foreach my $err_pair (sort {$err->{$b}{'links'}<=>$err->{$a}{'links'}} keys %$err){
      my $mean_iz = 0;
      $mean_iz = $err->{$err_pair}{'gaps'} / $err->{$err_pair}{'links'} if ($err->{$err_pair}{'links'});
      print PET "Pair $err_pair has $err->{$err_pair}{'links'} links and mean distance = $mean_iz\n";
   }
   close PET;

   my $satisfied = $ct_ok_pairs + $ct_ok_contig;
   my $unsatisfied = $ct_problem_pairs + $ct_iz_issues + $ct_illogical;
   my $ct_both_reads = $ct_both * 2;

   print SUMFILE "READ PAIRS STATS:\n";
   print SUMFILE "$seplines\tAt least one sequence/pair missing from contigs: $ct_single\n";
   print SUMFILE "\tAssembled pairs: $ct_both ($ct_both_reads sequences)\n";
   print SUMFILE "\t\tSatisfied in distance/logic within contigs (i.e. -> <-, distance on target: $insert_size +/$min_allowed): $ct_ok_contig\n";
   print SUMFILE "\t\tUnsatisfied in distance within contigs (i.e. distance out-of-bounds): $ct_iz_issues\n";
   print SUMFILE "\t\tUnsatisfied pairing logic within contigs (i.e. illogical pairing ->->, <-<- or <-->): $ct_illogical\n";
   print SUMFILE "\t\t---\n";
   print SUMFILE "\t\tSatisfied in distance/logic within a given contig pair (pre-scaffold): $ct_ok_pairs\n";
   print SUMFILE "\t\tUnsatisfied in distance within a given contig pair (i.e. calculated distances out-of-bounds): $ct_problem_pairs\n";
   print SUMFILE "\t\t---\n";
   print SUMFILE "\tTotal satisfied: $satisfied\tunsatisfied: $unsatisfied\n\n$seplines\n";

   #write distribution file
   open (CSV, ">$distribution") || die "Can't open $distribution for writing -- fatal";
   foreach my $is (sort {$a<=>$b} keys %$track_insert){
      print CSV "$is,$track_insert->{$is}\n";
   }
   close CSV;

   &FlushFiles();
   return $pair;
}

###READ BOWTIE OUTPUT AND STORE IT IN A HASH
sub readBowtie{
   my ($file) = @_;
   my $track_all;
   &printMessage("\n=>".getDate().": Storing Bowtie output ($file) to memory\n");

   open(IN,$file) || die "Can't open $file -- fatal\n";
   my $filefoundread = "$outdir/reads/$base_name.$library.foundpairedreads.fasta";
   open OUTFOUND, "> $filefoundread";
   my $lower = ($up_iz+200);
   my $sub = ($lower * 2) + 3;
   my ($prevline, $line, $prevread, $counter, $step) = ("","","",0, 100000);

   while(<IN>){
      chomp;
      $line = $_;
      ++$counter;
      if($counter == $step){
        CounterPrint($counter);
        $step = $step + 100000;
      }
      $line = $_;

      my ($read, $strand, $contig, $start, $seq, $other) = split(/\t/,$_);
      if($prevread eq $read){
        my ($seq1, $seq2);
        ($track_all, $seq1) = StoreResults($prevline, $track_all, $lower, $sub);
        ($track_all, $seq2) = StoreResults($line, $track_all, $lower, $sub);
        print OUTFOUND ">$counter\n$seq1:$seq2\n";
      }
      $prevread = $read;

      $prevline = $line;
   }
   close IN;
   close OUTFOUND;
   CounterPrint("                ");
   print SUMFILE "\nMAPPING READS TO CONTIGS:\n";
   print SUMFILE "$seplines\tNumber of single reads found on contigs = ".  keys( %$track_all )."\n";
   &FlushFiles();
   return $track_all, $filefoundread;
}

sub StoreResults{
  my ($input, $track_all, $lower, $sub) = @_;
  my ($read, $strand, $contig, $start, $seq, $other) = split(/\t/,$input);
  my ($tig,$startval, $endval);
      $tig = ($contig+1);
      if($start >  $lower && $tig_length->{$tig} > (($lower * 2)+100)){
        my $minsub = $sub - $start;
        $start = ($tig_length->{$tig} - $minsub);
      }
      my $upper;
      $upper = ($tig_length->{$tig} - ($up_iz + 100)) if($tig_length->{$tig} > ($up_iz +100));
      $upper = $lower if($tig_length->{$tig} <= ($up_iz +100));
      if($start <= ($lower) || $start >= ($upper)){

        $seq = reverseComplement($seq) if($strand eq "-");
        if($strand eq "+"){
           $startval = $start;
           $endval = $start + length($seq);
        }
        else{
           $startval = $start + length($seq);
           $endval = $start;
        }
        my $keyvalue =  "$tig"."|$startval"."|$endval\n";
        $track_all->{$seq} = $keyvalue;
      }
  return $track_all, $seq;
}

###FUNCTION TO REVERSE COMPLEMENT A SEQUENCE
sub reverseComplement{
   $_ = shift;
   tr/ATGC/TACG/;
   return (reverse());
}

###PRINTS A COUNTER ON THE SCREEN AND OVERWRITES PREVIOUS LINE
sub CounterPrint{
  my $countingMessager = shift;
  print "\r$countingMessager";
  $|++;
}

###FUNCTION TO PRINT MESSAGES TO THE SCREEN AND TO THE LOG FILE
sub printMessage{
  my $message = shift;
  print $message;
  print LOG $message;
}

###FUNCTION TO GET THE CURRENT DATE
sub getDate{
  my $date = scalar(localtime);
  return $date;
}

###FLUSHES THE SUMMARY AND LOG FILE
sub FlushFiles{
  select((select(SUMFILE), $| = 1)[0]);
  select((select(LOG), $| = 1)[0]);
  $|++;
}       

#########END PairingAndScaffolding.pl
