  ###############################################################
  #Marten Boetzer 23-10-2010                                    #
  #SSPACE perl subscript ExtendOrFormatContigs.pl               #
  #This script, based on the the -x parameter;                  #
  #  -Formats the contigs to appropriate format (-x 0)          #
  #  -Extends the contigs with available unmapped reads (-x 1)  #
  ###############################################################

  use strict;
  use File::Basename;

  my ($MAX, $MAX_TOP, $TRACK_COUNT) = (0, 100, 1);

  my $seplines = ("-" x 60)."\n";

  my $contig = $ARGV[0];
  my $base_name = $ARGV[1];
  my $extending = $ARGV[2];
  my $filecontig = $ARGV[3];
  my $MIN_READ_LENGTH = $ARGV[4];
  my $base_overlap = $ARGV[5];
  my $min_overlap = $ARGV[6];
  my $min_base_ratio = $ARGV[7];
  my $max_trim = $ARGV[8];
  my $verbose = $ARGV[9];
  my $Bin = $ARGV[10];
  my $outdir = $ARGV[11];
  my $seed;

  my $log = $outdir."/".$base_name . ".logfile.txt";
  my $summaryfile = $outdir."/".$base_name.".summaryfile.txt";

  open (SUMFILE, ">>$summaryfile") || die "Can't open $summaryfile -- fatal\n";
  open (LOG, ">>$log") || die "Can't write to logfile$log -- fatal\n";
  my $filenameOutExt = $base_name . ".singlereads.fasta";

  &ExtendContigs($base_name, $filecontig, $filenameOutExt, $Bin) if($extending == 1);
  &FormatContigs() if($extending == 0);

  close SUMFILE;
  close LOG;
#--------------------------------------------------

###EXTEND CONTIGS WITH UNMAPPED READS
sub ExtendContigs{
  my ($base_name, $filecontig, $filenameOutExt, $Bin) = @_;
  my ($set,$bin,$seq);
  
  my $encoded = &encodeBases();
  #-------------------------------------------------NOW MAP SINGLE READS TO INITIAL CONTIGS FILE.
  my $outfile = $outdir . "/reads/" . $filenameOutExt;
  my $lib = "start";
  system("perl $Bin/bin/mapWithBowtie.pl $base_name $filecontig $outfile $lib $Bin $outdir");
  #-------------------------------------------------READ UNMAPPED READS FOR EXTENSION
  ($set,$bin) = &readUnmappedReads($set,$bin,$base_name, $encoded);
  #-------------------------------------------------CONTIG EXTENSION USING UNMAPPED PAIRS STORED IN $SET
  &printMessage("\n=>".getDate().": Contig extension initiated\n");
  open (TIG, ">$contig") || die "Can't write to $contig -- fatal\n";
  #--------------------------------------------ASSEMBLY START

 ASSEMBLY:
   open(IN,$filecontig) || die "Can't open $filecontig -- fatal\n";
   my ($exttig_count, $counter, $NCount, $orig_mer, $prevhead) = (0,0,0,0,'');
   while(<IN>){
      s/\r\n/\n/;
      chomp;
      $seq.= uc($_) if(eof(IN));
      if (/\>(\S+)/ || eof(IN)){
         my $head=$1;
         $orig_mer = length($seq);
         if($seq ne ''){
             $NCount++ if($seq=~/([NX])/i);
             my $start_sequence = uc($seq);
             my $reads_needed = 1;                        #tracks coverage
             my $total_bases = $orig_mer * $reads_needed;

             ($seq, $set, $bin, $reads_needed, $total_bases) = doExtension("3", $orig_mer, $seq, $set, $bin, $reads_needed, $total_bases, $min_overlap, $base_overlap, $min_base_ratio, $verbose, $counter, $max_trim, $encoded) if($orig_mer >= $MIN_READ_LENGTH && $orig_mer >= $min_overlap);

             my $seqrc = reverseComplement($seq);
             ($seqrc, $set, $bin, $reads_needed, $total_bases) = doExtension("5", $orig_mer, $seqrc, $set, $bin, $reads_needed, $total_bases, $min_overlap, $base_overlap, $min_base_ratio, $verbose, $counter, $max_trim, $encoded) if($orig_mer >= $MIN_READ_LENGTH && $orig_mer >= $min_overlap);

             my $leng = length($seqrc);
             my $reversetig = reverseComplement($seqrc);                   ### return to sequence, as inputted
             my $trimmed_length = length($start_sequence) - 2*($max_trim);
             if($leng > $trimmed_length){ ### commented out: && $start_sequence ne $seqrc && $start_sequence ne $reversetig
               my $cov =  $total_bases / $leng;
               printf TIG ">extcontig%i|size%i|read%i|cov%.2f|seed:$prevhead\n%s\n", ($counter,$leng,$reads_needed,$cov,$reversetig);    #print contigs to file
               $exttig_count++;
             }else{
               my $cov = $reads_needed = 0;
               my $singlet_leng = length($start_sequence);
               printf TIG ">contig%i|size%i|read%i|cov%.2f|seed:$prevhead\n%s\n", ($counter,$leng,$reads_needed,$cov,$start_sequence);    #print singlets to file
             }
         }
         CounterPrint(++$counter);
         $prevhead = $head;
         $seq='';
      }else{
         $seq .= uc($_);
      }
   }
  CounterPrint("                ");
  print SUMFILE "\tNumber of contig sequences =".($counter-1). "\n";
  print SUMFILE "\t\tNumber of contigs containing N's (may prevent proper contig extension) = $NCount\n";

  print SUMFILE "\tNumber of contigs extended = $exttig_count\n".$seplines;
  close IN;

  if($@){
     my $message = $@;
     &printMessage("\nSomething went wrong running $0 ".getDate()."\n$message\n");
  }
  close TIG;

}

###STORE CONTIGS TO APPROPRIATE FORMAT WHEN CONTIGS WILL NOT BE EXTENDED
sub FormatContigs{
   &printMessage("\n=>".getDate().": Storing contigs to format for scaffolding\n");
   open (TIG, ">$contig") || die "Can't write to $contig -- fatal\n";
   open(IN,$filecontig) || die "Can't open $filecontig -- fatal\n";
   my ($counter, $seq, $prevhead) = (0,'','');
   while(<IN>){
      s/\r\n/\n/;
      chomp;
      $seq.= uc($_) if(eof(IN));
      if (/\>(\S+)/ || eof(IN)){
        my $head=$1;
        my $length_seq = length($seq);
        printf TIG ">contig%i|size%i|read%i|cov%.2f|seed:$prevhead\n%s\n", ($counter,$length_seq,0,0.00,$seq) if($seq ne '');
        CounterPrint(++$counter);
        $prevhead = $head;
        $seq = '';
      }else{
         $seq .= uc($_);
      }
   }
   CounterPrint("                ");
   print SUMFILE "CONTIG SUMMARY:\n".$seplines;
   print SUMFILE "\tNumber of unique contig sequences =".($counter-1). "\n".$seplines;
   close IN;
   close TIG;
}

###EXTEND CONTIGS
sub doExtension{

   my ($direction, $orig_mer, $seq, $set, $bin, $reads_needed, $total_bases, $min_overlap, $base_overlap, $min_base_ratio, $verbose, $tig_count, $max_trim, $e) = @_;

   my $previous = $seq;
   my ($extended, $trim_ct) = (1,0);

   if($orig_mer > $MAX){$orig_mer=$MAX;}  ### Deals with special cases where the seed sequences are different from the read set (and possibly very large) - goal here is not to increase sequence coverage of seed, but rather to extend it.

   TRIM:
   while($trim_ct <= $max_trim){
      while($extended){

         my ($pos,$current_reads,$current_bases,$span) = (0,0,0,"");

         ### Added 19March08
         if(length($seq) >= $MAX){   # $seq is length of contig being extended -- if larger than largest read, make sure the largest read could align and all subsequent rds.
            $span = $MAX - $TRACK_COUNT;
         }else{
            $span = length($seq) - $TRACK_COUNT;
         }
         my $startspan = $span;
         my $overhang = {};
         my @overlapping_reads = ();
         for (my $x=1;$x <= ($orig_mer * 2);$x++){
            ($overhang->{$x}{'A'},$overhang->{$x}{'C'},$overhang->{$x}{'G'},$overhang->{$x}{'T'}) = (0,0,0,0);
         }

         ### COLLECT SEQUENCES 
         while ($span >= $min_overlap){  # will slide the subseq, until the user-defined min overlap size

            $pos = length($seq) - $span;
            print "MAX:$MAX, SPAN:$span, POS:$pos" if ($verbose);

            my $subseq = substr($seq, $pos, $span);              #make a sub-sequence of length l-(1..i) for searching
            my $sub = substr($subseq, 0, 10);                    #grab first 10 nucleotides and get all reads having this subset stored in $bin
            my $subset = $bin->{$sub};                           #Will grab everything even the reverse complement ones
            print "####$direction' SEARCH Position:$pos Span:$span - Subseq:$subseq Previous:$previous\n" if ($verbose);
            ### SEARCH -- this cycles through limited k-mer space
            foreach my $pass (sort {$subset->{$b} <=> $subset->{$a}} keys %$subset){
               if($pass =~ /^$subseq([ACGT]*)/){
                  #can we align perfectly that subseq to another rd start?
                  my $dangle = $1;
                  print "\n", "=" x 80, "\n$direction'- FOUND sequence: $pass -> subset: $subseq -> overhang: $dangle\n", "=" x 80, "\n\n" if ($verbose);

                  # Collect all overhangs
                  push @overlapping_reads, $pass;                  ### all overlapping reads
                  my @over = split(//,$dangle);
                  my $ct_oh = 0;
 
                  foreach my $bz(@over){
                     $ct_oh++;                                     ### tracks overhang position passed the seed  
                     if(defined $set->{$pass}){
                        $overhang->{$ct_oh}{$bz} += $set->{$pass};      ### reflects read coverage (often real duplicates)
                     }else{
                        my $pass_rc = reverseComplement($pass);
                        $overhang->{$ct_oh}{$bz} += $set->{$pass_rc};
                     }
                     print "$ct_oh - $bz = $overhang->{$ct_oh}{$bz}\n" if($verbose);
                  }
               }
            }
            $span--;
         }#while overlap >= user-defined -m minimum

         my $consensus = "";
         print "Finished Collecting Overlapping Reads - BUILDING CONSENSUS...\n" if ($verbose);
         print Dumper(@overlapping_reads) if ($verbose);

         ### Build consensus
         CONSENSUS:
         foreach my $ohpos (sort {$a<=>$b} keys %$overhang){
            if($ohpos){

               my $coverage = $overhang->{$ohpos}{'A'}+$overhang->{$ohpos}{'C'}+$overhang->{$ohpos}{'G'}+$overhang->{$ohpos}{'T'};
               print "pos:$ohpos cov:$coverage A:$overhang->{$ohpos}{'A'} C:$overhang->{$ohpos}{'C'} G:$overhang->{$ohpos}{'G'} T:$overhang->{$ohpos}{'T'}\n" if($verbose);
               if ($coverage < $base_overlap){
                  print "COVERAGE BELOW THRESHOLD: $coverage < -o $base_overlap @ $ohpos :: will extend by: $consensus\n" if ($verbose);
                  last CONSENSUS;
               }

               my $baselist = $overhang->{$ohpos};

               my ($ct_dna, $previous_bz) = (0, "");

               BASE:
               foreach my $bz (sort {$baselist->{$b}<=>$baselist->{$a}} keys %$baselist){
                 if($ct_dna){## the two most abundant bases at that position
                     if($previous_bz ne "" && ($baselist->{$previous_bz} / $coverage) >= $min_base_ratio && $baselist->{$previous_bz} > $baselist->{$bz}){### a simple consensus btw top 2
                        $consensus .= $previous_bz;                                         ### build consensus
                        print "Added base $previous_bz (cov = $baselist->{$previous_bz}) to $consensus **\n" if ($verbose);
                        last BASE;
                     }else{
                        print "ISSUES EXTENDING: best base = $previous_bz (cov=$baselist->{$previous_bz}) at $ohpos.  Second-Best: $bz (cov=$baselist->{$bz}) (ratio best=$baselist->{$previous_bz} / total=$coverage) >= $min_base_ratio (-r) -- will terminate with $consensus\n" if($verbose);
                        last CONSENSUS;
                     }
                  }
                  $previous_bz = $bz;                 
                  $ct_dna++;
               }
            }
         }

         ### deal with sequence reads making up the consensus/newly formed contig
         if($consensus ne ""){

            print "Will extend $seq\nwith: $consensus\n\n" if($verbose);
            my $temp_sequence = $seq . $consensus;  ## this is the contig extension  
            my $integral = 0;
            my $position = length($temp_sequence) - ($startspan + length($consensus));
            my $temp_sequence_end = substr($temp_sequence, $position);
            foreach my $ro (@overlapping_reads){
             my $result = index($temp_sequence_end, $ro);

               if($temp_sequence_end =~ /$ro/){
                  my $complement_ro = reverseComplement($ro);
                  $integral=1;

                 if(defined $set->{$ro}){
                     $current_reads = $set->{$ro};
                     $current_bases = length($ro) * $current_reads;
                     $reads_needed += $current_reads;
                     $total_bases += $current_bases;
                     ($bin,$set,$seed) = deleteData($bin,$set,$seed,$ro,$e)
                  }

                  if(defined $set->{$complement_ro}){
                     $current_reads = $set->{$complement_ro};
                     $current_bases = length($complement_ro) * $current_reads;
                     $reads_needed += $current_reads;
                     $total_bases += $current_bases;
                     ($bin,$set,$seed) = deleteData($bin,$set,$seed,$complement_ro,$e);
                  }
               }
            }
            if(! $integral){### no reads are found overlapping with the consensus might be indicative of low complexity regions -- Stop the extension
               print "No overlapping reads agree with the consensus sequence.   Stopping extension" if ($verbose);
               $extended = 0;
            }else{
               $seq = $temp_sequence;
               $temp_sequence = "";
               print "New Contig is: $seq\n" if ($verbose);
               $extended = 1;
            }
            $previous = $seq;
         }else{### no consensus built, will stop the extension
            $extended = 0;
         }

      }###while get the OK for extension

      $trim_ct++;
      if ($trim_ct <= $max_trim){
         last TRIM if (length($seq) <= $MIN_READ_LENGTH); #terminate assembly if trimming becomes too agressive
         $seq = substr($seq, 0, -1);
         $extended = 1;
         print "\n$direction prime EXTENSION ROUND $trim_ct COMPLETE UNTIL $max_trim nt TRIMMED OFF => TRIMMED SEQUENCE:$seq\n\n" if ($verbose);
      }
      
   }### while trimming within bounds

   print "\n*** NOTHING ELSE TO BE DONE IN $direction prime- PERHAPS YOU COULD DECREASE THE MINIMUM OVERLAP -m (currently set to -m $min_overlap) ***\n\n" if ($verbose);

   return $seq, $set, $bin, $reads_needed, $total_bases;
}


###DELETE READ DATA IF IT HAS BEEN USED FOR EXTENDING A CONTIG
sub deleteData {
   my ($bin,$set,$seed,$sequence,$e) = @_;
   
   my $subnor = substr($sequence, 0, 10);
   my $comp_seq = reverseComplement($sequence);
   my $subrv = substr($comp_seq, 0, 10);

   #remove k-mer from hash table and prefix tree
   delete $bin->{$subrv}{$comp_seq};
   delete $bin->{$subnor}{$sequence};
   delete $set->{$sequence};
   return $bin, $set, $seed;
}

sub readUnmappedReads{
  my ($set,$bin,$base_name, $encoded) = @_;

  &printMessage("\n=>".getDate().": Reading unmapped reads to memory \n");

  my $outfileExt =  $base_name . ".start.unmapped";
  my $file = $outdir."/"."bowtieoutput/$outfileExt";
  my ($counter, $step) = (0, 100000);
  my $e = $encoded;
  if(-e $file){
    open(INUNMAPPED,$file);
    while(<INUNMAPPED>){
      s/\r\n/\n/, chomp;
      if(!(/\>(\S+)/)){
        ++$counter;
        if($counter == $step){
          CounterPrint($counter);
          $step = $step + 100000;
        }
        my $orig = $_;
        my $orig_mer = length($orig);
        if ($orig ne '' && $orig_mer >= $MIN_READ_LENGTH && $orig_mer >= $min_overlap){
          $_ = $orig;
          tr/ACTG/TGAC/;
          my $rc=reverse();

          $MAX=$orig_mer if ($orig_mer > $MAX);

          my $subrv = substr($rc, 0, 10);
          my $subnor = substr($orig, 0, 10);

          $set->{$orig}++;
          $bin->{$subnor}{$orig}++;
          $bin->{$subrv}{$rc}++;
        }
        $MAX = $MAX_TOP if ($MAX > $MAX_TOP);
      }
    }
  }
  else{
    print "WARNING: NO SEQUENCES USED FOR EXTENSION\n\n";
  }
  close INUNMAPPED;
  CounterPrint("                ");

  print SUMFILE "CONTIG EXTENSION:\n".$seplines;
  print SUMFILE "\tNumber of unmapped reads used for contig extension = $counter\n";
  &FlushFiles();
  return ($set,$bin);
}


#CONVERT EACH THREE NUCLEOTIDE CODON TO COMPACT FORMAT
sub encodeBases{
   my $encoded;
   my @pos1= ('A','C','G','T');
   my @pos2 = @pos1;
   my @pos3 = @pos1;
 
   my $ascii_num = 33;

   foreach my $p1 (@pos1){
      foreach my $p2 (@pos2){
         foreach my $p3 (@pos3){
            my $codon = $p1.$p2.$p3;
            my $char = chr($ascii_num);
            $encoded->{$codon}=$char;
            $ascii_num++;
         } 
      }
   }
   return $encoded;
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

#########END ExtendOrFormatContigs.pl
