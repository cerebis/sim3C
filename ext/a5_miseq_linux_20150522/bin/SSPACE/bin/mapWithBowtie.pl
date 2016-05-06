  ####################################################################
  #Marten Boetzer 23-10-2010                                         #
  #SSPACE perl subscript mapWithBowtie.pl                            #
  #This script;                                                      #
  #  -Calls the external program Bowtie to map the reads to contigs  #
  ####################################################################

  use strict;
  use File::Basename;

  my $base_name = $ARGV[0];
  my $contigFile = $ARGV[1];
  my $singlereads = $ARGV[2];
  my $library = $ARGV[3];
  my $Bin = $ARGV[4];
  my $outdir = $ARGV[5];

  my $log = "$outdir/$base_name.logfile.txt";
  open (LOG, ">>$log") || die "Can't write to $log -- fatal\n";
  &MapReadsToContigs($base_name,$contigFile, $singlereads, $library);
  close LOG;
#--------------------------------------------------

###MAP SINGLE READS TO CONTIGS WITH BOWTIE
sub MapReadsToContigs{
    my ($base_name, $contigFile, $singlereads, $library) = @_;

    #my $bowtieout = $outdir."/".$base_name . ".$library.bowtieIndex";
	# Andrew Tritt 8/10/2011: bowtieout was constructed incorrectly
    my $bowtieout = $base_name . ".$library.bowtieIndex";
    my $bowbuildpath = "$Bin"."/bowtie/bowtie-build";
    my $bowtiepath = "$Bin"."/bowtie/bowtie";
    $bowtiepath =~ s/ /\\ /g;
    $bowbuildpath  =~ s/ /\\ /g;
    #my $outfileExt =  $outdir."/".$base_name . ".$library.unmapped";
    #my $outfileNotExt =  $outdir."/".$base_name . ".$library.mapped";
	# Andrew Tritt 8/10/2011: outfileExt and outfileNotExt were constructed incorrectly. Same as bowtieout
    my $outfileExt =  $base_name . ".$library.unmapped";
    my $outfileNotExt =  $base_name . ".$library.mapped";
    
    die "Contig file ($contigFile) not found. Exiting...\n" if(!(-e $contigFile));
    &printMessage("\n=>".getDate().": Building Bowtie index for contigs ($contigFile)\n");
    system("$bowbuildpath $contigFile $outdir/bowtieoutput/$bowtieout --quiet --noref") == 0 || die "\nBowtie-build error; $?"; # returns exit status values
    die "Single read file ($singlereads) not found. Exiting...\n" if(!(-e $singlereads));
    &printMessage("\n=>".getDate().": Mapping reads ($singlereads) to Bowtie index\n");
    system("$bowtiepath -v 0 $outdir/bowtieoutput/$bowtieout --suppress 2,3,4,6,7 -f $singlereads $outdir/bowtieoutput/$outfileNotExt --un $outdir/bowtieoutput/$outfileExt --quiet --refidx") == 0 || die "\nBowtie error; $?" if($library eq "start"); # returns exit status values
    system("$bowtiepath -v 0 -m 1 $outdir/bowtieoutput/$bowtieout --suppress 6,7 -f $singlereads $outdir/bowtieoutput/$outfileNotExt --quiet --refidx") == 0 || die "\nBowtie error; $?" if($library ne "start"); # returns exit status values
    &FlushFiles();
}

###FLUSHES THE SUMMARY AND LOG FILE
sub FlushFiles{
  select((select(LOG), $| = 1)[0]);
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

#########END MapWithBowtie.pl
