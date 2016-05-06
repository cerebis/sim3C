  #############################################################
  #Marten Boetzer 23-11-2010                                  #
  #SSPACE perl subscript readLibFiles.pl                      #
  #This script;                                               #
  #  -reads, converts and filters original input sequences    #
  #############################################################

  use strict;
  use Storable;
  use File::Path;
  use File::Basename;

  my $seplines = ("-" x 60)."\n";

  my $libraryfile = $ARGV[0];
  my $base_name = $ARGV[1];
  my $extending = $ARGV[2];
  my $unpaired_file = $ARGV[3];
  my $outdir = $ARGV[4];
  
  my $log = $outdir."/".$base_name . ".logfile.txt";
  my $summaryfile = $outdir."/".$base_name.".summaryfile.txt";

  open (SUMFILE, ">>$summaryfile") || die "Can't open $summaryfile -- fatal\n";
  open (LOG, ">>$log") || die "Can't write to $log -- fatal\n";
  
  my $filenameOutFilt = "filtered.readpairs.fasta";
  my $filenameOutExt = $base_name . ".singlereads.fasta";

  open OUTFILEExt, "> $outdir/reads/$filenameOutExt" if($extending == 1);

#-------------------------------------------------READ UNPAIRED FILE CONTAINING SINGLE READS
  &readUnpairedFile($unpaired_file) if ($unpaired_file);
#-------------------------------------------------LOOP THROUGH EACH LIBRARY IN LIBRARYFILE AND STORE AND FILTER READS
  open(FILELIB, "< $libraryfile");
  
  my ($library, $fileA, $fileB, $insert_size, $insert_stdev, $reverse);
  my ($totcounter, $totNcount, $prevlibrary) = (0,0, "");
  my ($filesA, $filesB, $inserts, $stdevs, $reverses);
  while(<FILELIB>){
    chomp;
    ($library, $fileA, $fileB, $insert_size, $insert_stdev, $reverse) = split(/\s+/, $_);

    next if($library eq "");
    my ($fileBaseName1, $dirName1, $fileExtension1) = fileparse($fileA);
    my ($fileBaseName2, $dirName2, $fileExtension2) = fileparse($fileB);

    if($library ne $prevlibrary){
      if($prevlibrary ne ""){
        my $filt = $totcounter-$totNcount;
        print SUMFILE "READING READS $prevlibrary:\n";
        print SUMFILE "$seplines\tTotal inserted pairs = $totcounter \n";
        print SUMFILE "\tNumber of pairs containing N's = $totNcount \n\tRemaining pairs = $filt\n$seplines\n";

        &printMessage("\n=>".getDate().": Reading, filtering and converting input sequences of library '$library' initiated\n");
        ($filesA, $filesB, $inserts, $stdevs, $reverses) = ("", "","","","");
        &FlushFiles();
        ($totNcount, $totcounter) = (0,0);
        close OUTSINGLEFILE;
      }elsif($prevlibrary eq ""){
        &printMessage("\n=>".getDate().": Reading, filtering and converting input sequences of library '$library' initiated\n");
      }
     open (OUTSINGLEFILE, ">$outdir/reads/$base_name.$library.filtered.readpairs.singles.fasta") || die "Can't write to single file -- fatal\n";

      CounterPrint("                ");
    }
    print "Reading: $fileBaseName1 and $fileBaseName2...\n";

    $filesA .= "$fileA |  ";
    $filesB .= "$fileB |  ";
    $inserts .= "$insert_size |  ";
    $stdevs .= "$insert_stdev |  ";
    $reverses .= "$reverse |  ";
    my ($counter2, $Ncount2) = &generateInputFiles($library, $fileA, $fileB, $extending, $reverse);
    $totcounter += $counter2;
    $totNcount += $Ncount2;

    $prevlibrary = $library;
  }

  if($prevlibrary ne ""){
    my $filt = $totcounter-$totNcount;
    print SUMFILE "READING READS $prevlibrary:\n";
    print SUMFILE "$seplines\tTotal inserted pairs = $totcounter \n";
    print SUMFILE "\tNumber of pairs containing N's = $totNcount \n\tRemaining pairs = $filt\n$seplines\n";
  }

  &printMessage("\n$seplines");

  close OUTSINGLEFILE;
  close FILELIB;
  close OUTFILEExt if($extending == 1);
  close SUMFILE;
  close LOG;

#--------------------------------------------------

###CONVERT INPUT SEQUENCES BY REMOVING PAIRED READS HAVING AN 'N' OR DUPLICATE PAIRS
sub generateInputFiles{
  my ($lib, $fileA, $fileB, $extension, $reverse) = @_;

  my ($combined, $name,$seq1,$seq2);
  my ($counterext, $Ncount, $counter, $countsinglet, $fastq, $step) = (0,0,0,0,0,100000);

  open(TEST, "< $fileA");
  $name = <TEST>;
  close TEST;
  $fastq = 1 if ($name =~ /^[@]/);

  open(FILEA, "< $fileA");
  open(FILEB, "< $fileB");

  while(<FILEA>) {
    <FILEB>;
    $seq1 = uc(<FILEA>), $seq1 =~ s/\r\n/\n/, chomp $seq1;
    $seq2 = uc(<FILEB>), $seq2 =~ s/\r\n/\n/, chomp $seq2;
    #FASTQ FORMAT
    <FILEA>,<FILEA>,<FILEB>,<FILEB> if ($fastq);
    # ELSE FASTA FORMAT
    $seq1 = reverseComplement($seq1), $seq2 = reverseComplement($seq2) if($reverse);

    if ($extension == 1){
        if($seq1 =~ /^([ACGT]*)$/i){
         print OUTFILEExt ">$counterext\n$seq1\n";
         $counterext++;
      }
      if ($seq2 =~ /^([ACGT]*)$/i){
         print OUTFILEExt ">$counterext\n$seq2\n";
         $counterext++;
      }
    }
    $combined = "$seq1:$seq2";
    if(!($combined =~ /[N]/)){
      ++$countsinglet;
       print OUTSINGLEFILE ">read$countsinglet\n$seq1\n";
       print OUTSINGLEFILE ">read$countsinglet\n$seq2\n";
    }
    else{
      $Ncount++;
    }
    ++$counter;
    if($counter == $step){
      CounterPrint($counter);
      $step = $step + 100000;
    }
  }
  CounterPrint("                ");
  close FILEA;
  close FILEB;
  return $counter, $Ncount;
}

#------------------READ UNPAIRED SINGLE READS FILE WHEN -u IS SET

sub readUnpairedFile{
  my ($file) = @_;
  open(INUNPAIRED, "< $file") || die "Can't open $file -- fatal\n";

  &printMessage("\n=>Reading, filtering and converting unpaired input sequences initiated ".getDate()."\n");

  my ($seq1, $name);
  my ($counterext, $counter, $step) = (0,0, 100000);
  while(<INUNPAIRED>) {
    chomp;
    $name = $_;
    $seq1 = uc(<INUNPAIRED>); $seq1 =~ s/\r\n/\n/; chomp $seq1;

    #FASTQ FORMAT
    if ($name =~ /^[@]/){
      <INUNPAIRED>; <INUNPAIRED>;
    }
    # ELSE FASTA FORMAT
    
    if ($seq1 =~ /^([ACGT]*)$/i){
       print OUTFILEExt ">$counterext\n$seq1\n";
       $counterext++;
    }
    ++$counter;
    if($counter == $step){
      CounterPrint($counter);
      $step = $step + 100000;
    }
  }
  CounterPrint("                ");

  print SUMFILE "READING UNPAIRED READS:\n";
  print SUMFILE "$seplines\tTotal inserted reads = $counter \n";
  print SUMFILE "\tNumber of reads containing N's = ".($counter-$counterext)."\n\tRemaining reads = $counterext\n";

  close INUNPAIRED;
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

#########END readLibFiles.pl
