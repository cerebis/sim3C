#!/usr/bin/env perl
#
# a5 pipeline
# (c) 2011-2014 Andrew Tritt and Aaron Darling
# This is free software licensed under the GPL, v3. Please see the file LICENSE for details.
#
# Usage: a5_pipeline.pl <library file> <output directory or basename>
#   Or:  a5_pipeline.pl <fastq 1> <fastq 2> <output directory or basename>
#
use strict;
use warnings;
use File::Basename;
use Cwd 'abs_path';
use Getopt::Long;

my $VERSION="20160825";

=pod

=head1 NAME

a5_pipeline.pl -- Assemble isolate genomes from Illumina data with ease

=head1 SYNOPSIS

a5_pipeline.pl takes FastQ format sequence reads directly from the Illumina base-calling software and cleans, filters, and assembles them.
For example the command:

a5_pipeline.pl read1.fastq read2.fastq my_assembly

Will assemble the paired reads in read1.fastq and read2.fastq and store the result in files whose names begin with "my_assembly".
The final scaffolded assembly will be named my_assembly.final.scaffolds.fasta

=head1 DESCRIPTION

The flow of execution is roughly:
   1) Clean up the sequence data
      - Error correct, trim adapter, quality trim
   2) Build contigs on error-corrected reads
   3) Scaffold contigs using paired-read and mate-pair information
   4) Detect any misassemblies and break contigs/scaffolds
   5) Rescaffold the broken contigs/scaffolds

Throughout the pipeline there are various steps used to estimate
alignment parameters

=head1 EXAMPLES

With a single paired-end illumina library:
a5_pipeline.pl <read1.fastq> <read2.fastq> <output directory or basename>

With two or more libraries:
a5_pipeline.pl <library file> <output directory or basename>

=head1 AUTHORS

Andrew Tritt <ajtritt@lbl.gov>
Aaron Darling <aaron.darling@uts.edu.au>

=head1 AVAILABILITY

http://sourceforge.net/projects/ngopt

=head1 COPYRIGHT

Copyright 2011-2014
This software is licensed under the terms of the GPLv3

=cut

use constant {
	IDBA_MIN_K => 35,
	SGA_Q_TRIM => 25,
	SGA_Q_FILTER => 20,
	SGA_DUST => 8,
	SGA_MIN_READ_LENGTH => 35,
	MAX_PAIREDEND_INSERT_MEAN => 1500,
	JAVA_HEAP => -Xmx512m,
};

my $DIR = shell_safe(dirname(abs_path($0)));
my $AVAILMEM = 4000000;
my $def_up_id="upReads";
my @KEYS = ("id","p1","p2","shuf","up","rc","ins","err","nlibs","libfile");
my $pname = basename($0);
my $usage= qq{
A5-miseq version $VERSION
Usage: $pname [--begin=1-5] [--end=1-5] [--preprocessed] [--threads=4] [--debug] [--metagenome] <lib_file> <out_base>

Or: $pname <Read 1 FastQ> <Read 2 FastQ> <out_base>

Or: $pname <Read 1,2 Interleaved FastQ> <out_base>

<out_base> is the base file name for all output files. When assembling from 
a single library, the fastq files may be given directly on the command line.
If using more than one library, a library file must be given as <lib_file>.
The library file must contain the filenames of all read files.

If --preprocessed is used, <lib_file> is expected to be the library file
created during step 2 of the pipeline, named <out_base>.preproc.libs. Note 
that this flag only applies if beginning pipeline after step 2.
};
die $usage if ! @ARGV;

Getopt::Long::Configure(qw{no_auto_abbrev no_ignore_case_always pass_through});
my $start = 1;
my $end = 2;
my $preproc = 0;
my $debug = 0;
my $metagenome = 0;
my $threads = 4;
my $adapter = dirname(abs_path($0))."/../adapter.fasta";
GetOptions( 'begin=i' => \$start,
            'end=i' => \$end,
			'preprocessed' => \$preproc,
			'debug' => \$debug,
			'metagenome' => \$metagenome,
			'threads=i' => \$threads,
			'adapter=s' => \$adapter);

die $usage if (@ARGV < 2);

$AVAILMEM = get_availmem();

#
# Check that java is available and in the path
#
my $java = `which java`;
die "Unable to find a java runtime environment. The A5 pipeline requires java 6 or later. Please ensure that java is installed and in the \$PATH" unless length($java)>1;
# OS compatibility checks
if($^O !~ /darwin/ && $^O !~ /inux/){
	die "Sorry, A5-miseq only works on Linux and OS X platforms\n";
}
if($^O =~ /darwin/){
	my $java_version = `java -version 2>&1`;
	if($java_version !~ /java version/){
		die "Sorry, in order to run A5 you must first install command-line Java for Mac OS X. It is available from http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html\n";
	}
	my $sga_test = system("$DIR/sga 2> /dev/null > /dev/null");
	if($sga_test){
		die "Sorry, A5-miseq does not work on your system.\nThe SGA component of A5-miseq requires OS X 10.6 or later\n";
	}
}

# check for spaces in file and path names -- these are not supported
my $message = "contains unsupported characters such as quotes, spaces, or &:%?*><\$";

my $libfile = $ARGV[0];
my $OUTBASE = $ARGV[1];

# check whether command-line input was FastQ files or a library file
if(@ARGV==3){
	# assume the user provided two FastQ files instead of a library file
	# first check the files exist
	validate_sequence_file($ARGV[0], $message);
	validate_sequence_file($ARGV[1], $message);
	$OUTBASE = $ARGV[2];
	open(TMPLIBFILE, ">$OUTBASE.tmplibs");
	print TMPLIBFILE "[LIB]\n";
	print TMPLIBFILE "p1=$ARGV[0]\n";
	print TMPLIBFILE "p2=$ARGV[1]\n";
	close TMPLIBFILE;
	$libfile = "$OUTBASE.tmplibs";
} else { 
	$OUTBASE = $ARGV[1];
	my $file = $ARGV[0];
	my $file_type = `file $file`;
	my $first_line = "";
	if ($file_type =~ /gzip/){
		$first_line = `gunzip -c $file | head -n 1`;
	} elsif ($file_type =~ /bzip2/) {
		$first_line = `bunzip2 -c $file | head -n 1`;
	} else {
		$first_line = `head -n 1 $file`;
	}
	if ($first_line =~ /^@/){ # assume interleaved
		open(TMPLIBFILE, ">$OUTBASE.tmplibs");
		print TMPLIBFILE "[LIB]\n";
		print TMPLIBFILE "shuf=$ARGV[0]\n";
		close TMPLIBFILE;
		$libfile = "$OUTBASE.tmplibs";
	} elsif ($first_line =~ /^\[LIB\]/) {
		$libfile = $ARGV[0];
	} else {
		print STDERR "$file is neither a library file nor a fastq file.\n";
		exit;
	}
}

die "\nERROR: use of an absolute path for the output base name is not supported.\n\"$OUTBASE\" is an absolute path.\nPlease specify a relative path instead.\n" if $OUTBASE =~ /^\//;
die "ERROR: The output location provided \"$OUTBASE\" $message.\nPlease provide an alternate output location\n" unless validate_filename($OUTBASE);
die "ERROR: The input file $libfile $message.\nPlease provide an alternate output location\n" unless validate_filename($libfile);

if(-f $OUTBASE && !-d $OUTBASE){
	die "$OUTBASE already exists and is not a directory.\nPlease specify a new file path for a base name after the library file or FastQ files.\n";
}

# Format for output of time resource statistics
my $TIME_FMT = "###PROCESS\n".
               "%C\n".
               "###\n".
               "%Uuser %Ssystem %Eelapsed %PCPU (%Xtext+%Ddata %Mmax)k\n".
               "%Iinputs+%Ooutputs (%Fmajor+%Rminor)pagefaults %Wswaps";
my $STAT_FILE = "$OUTBASE.res_stat";
my $TIME = "/usr/bin/time -o $STAT_FILE -a -f $TIME_FMT ";

#variables used for stats generation
my ($total_reads,$total_bases,$EC_reads,$EC_bases) = (0,0,0,0);

# An extension to give unzipped files.
my $UNZIP_EXT = "a5unzipped";
my %RAW_LIBS = read_lib_file($libfile);
print STDERR "[a5] Found the following libraries:\n";
print_lib_info(\%RAW_LIBS);
for my $lib (keys %RAW_LIBS){
	for my $att (keys %{$RAW_LIBS{$lib}}){
		$RAW_LIBS{$lib}{$att} = check_and_unzip($RAW_LIBS{$lib}{$att}) if ($att eq "p1" || $att eq "p2" || $att eq "shuf");
		die "ERROR: The file \"".$RAW_LIBS{$lib}{$att}."\" $message.\nPlease ensure that all sequence input file names are free from unsupported characters\n" unless validate_filename($RAW_LIBS{$lib}{$att});
		if ($att eq "shuf"){
			my ($fq1, $fq2) = split_shuf($RAW_LIBS{$lib}{$att},"$OUTBASE.".$RAW_LIBS{$lib}{"id"});
			$RAW_LIBS{$lib}{"p1"} = $fq1;
			$RAW_LIBS{$lib}{"p2"} = $fq2;
		}
	}
}

die "[a5] No libraries found in $libfile\n" unless(keys %RAW_LIBS);
print "[a5] Found ".scalar(keys %RAW_LIBS)." libraries\n";
my $maxrdlen = -1;
my $scafs;
my $ctgs;
my $reads;
my $WD="$OUTBASE";

print "[a5] Starting pipeline at step $start\n";

# if we are assembling metagenomic data do not run the scaffolding
$end = ($end < 2 ? $end : 2) if $metagenome;

if ($start <= 1) {
	print "[a5] Cleaning reads with SGA\n";
	print STDERR "[a5] Cleaning reads with SGA\n";
	$WD="$OUTBASE.s1";
	mkdir($WD) if ! -d $WD;
	$reads = sga_clean($OUTBASE, \%RAW_LIBS);
	`mv $reads $OUTBASE.ec.fastq`;
	`rm -f $WD/*.fastq`;
	`rm -f $WD/*.ec.fa`;
} 
if ($end == 1){
	print STDERR "[a5] Done cleaning reads. Results at $OUTBASE.ec.fastq\n";
	exit;
}
if ($start <= 2) {
	my $fq_reads = "$OUTBASE.ec.fastq";
	if (-f "$fq_reads.gz" && ! -f "$fq_reads") { 
		`gunzip $fq_reads.gz`;
	}
	$reads = $fq_reads;
	die "[a5_s2] Can't find error corrected reads $reads" unless -f $reads && (-s $reads > 0);
	print "[a5_s2] Building contigs from $reads with IDBA\n";
	print STDERR "[a5_s2] Building contigs from $reads with IDBA\n";
	$WD="$OUTBASE.s2";
	mkdir($WD) if ! -d $WD;
	($reads, $maxrdlen) = fastq_to_fasta($reads,"$WD/$OUTBASE.ec.fasta");
	`gzip -f $fq_reads`;

	$ctgs = idba_pe_assemble($OUTBASE, \%RAW_LIBS, $maxrdlen);
	`rm -f $reads`;
	open (my $filtered, ">", "$OUTBASE.contigs.fasta");
	open (IN,"<",$ctgs);
	while (my $hdr = <IN>){
		my $seq = <IN>;
		chomp $seq;
		next if (length($seq) < 2*$maxrdlen); 
		print $filtered $hdr.$seq."\n";
	}
	close $filtered;
	`rm $ctgs`;
	#`mv $ctgs $OUTBASE.contigs.fasta`;
	
} 
if ($end == 2 && !$metagenome){
	my $stats = map_libs_and_polish($OUTBASE, \%RAW_LIBS, ".contigs");
	finalize($stats, "$OUTBASE.contigs.fasta");
	print STDERR "[a5] Done building assembly. Results at $OUTBASE.contigs.fasta\n";
	print STDERR "[a5] Note, large insert scaffolding skipped (can be enabled with --end=5). Contig and scaffold statistics will be identical.\n";
	exit;
}
$WD="$OUTBASE.s3";
mkdir($WD) if ! -d $WD;
my %PAIR_LIBS; 
if (!$preproc || $start <= 2){
	print STDERR "[a5] Preprocess libraries for scaffolding with SSPACE\n";
	%PAIR_LIBS = preprocess_libs(\%RAW_LIBS,"$OUTBASE.contigs.fasta");
	print STDERR "[a5] Processed libraries:\n";
	print_lib_info(\%PAIR_LIBS);
} else {
	print STDERR "[a5] Libraries already preprocessed\n" if ($preproc);
	%PAIR_LIBS = %RAW_LIBS;
	delete($PAIR_LIBS{$def_up_id}) if (defined($PAIR_LIBS{$def_up_id}));
}
if ($start <= 3) {
	$ctgs = "$OUTBASE.contigs.fasta";
	die "[a5_s3] Can't find starting contigs $ctgs.\n" unless -f $ctgs;
	print "[a5_s3] Scaffolding contigs from $ctgs with SSPACE\n";
	print STDERR "[a5_s3] Scaffolding contigs from $ctgs with SSPACE\n";
	$scafs = scaffold_sspace($OUTBASE, \%PAIR_LIBS, $ctgs);
	`mv $scafs $OUTBASE.crude.scaffolds.fasta`; 
} 
if ($end == 3){
	print STDERR "[a5] Done building crude scaffolds. Results at $OUTBASE.crude.scaffolds.fasta\n";
	exit;
}
my $need_qc = 1;
if ($start <= 4) {
	my $prev_scafs = "$OUTBASE.crude.scaffolds.fasta";
	$scafs = $prev_scafs;
	die "[a5_s4] Can't find starting crude scaffolds $scafs\n" unless -f $scafs;
	print STDERR "[a5_s4] Detecting and breaking misassemblies in $scafs with A5QC\n";
	$WD="$OUTBASE.s4";
	mkdir($WD) if ! -d $WD;
	$scafs = break_all_misasms($prev_scafs,\%PAIR_LIBS,"$OUTBASE.qc"); 
	if ($scafs eq $prev_scafs){
		$need_qc = 0;
		`cp $scafs $OUTBASE.final.scaffolds.fasta`;
		`cp $scafs $OUTBASE.broken.scaffolds.fasta`; # so that pipeline can begin at step 5
		print "[a5_s5] No misassemblies found.\n";
	} else {
		`mv $scafs $OUTBASE.broken.scaffolds.fasta`; 
	}
} 
if ($end == 4){
	print STDERR "[a5] Done running QC. ";
	if ($need_qc){
		print STDERR " Results at $OUTBASE.broken.scaffolds.fasta\n";
	} else {
		print STDERR " No misassemblies detected.";
	}
	exit;
}
if ($start <= 5 && $end >= 5) {
	if($need_qc) {
		$scafs = "$OUTBASE.broken.scaffolds.fasta";
		die "[a5_s5] Can't find starting broken scaffolds $scafs\n" unless -f $scafs;
		print "[a5_s5] Scaffolding broken contigs with SSPACE\n";
		print STDERR "[a5_s5] Scaffolding broken contigs with SSPACE\n";
		$WD="$OUTBASE.s5";
		mkdir($WD) if ! -d $WD;
		$scafs = scaffold_sspace("$OUTBASE.rescaf",\%PAIR_LIBS,$scafs,1);
		`mv $scafs $OUTBASE.final.scaffolds.fasta`;
	}
	# get rid of |size... in FastA header because NCBI doesn't like it
	`perl -p -i -e "s/\\|size.+//g" $OUTBASE.final.scaffolds.fasta`;
	my $stats = map_libs_and_polish($OUTBASE, \%RAW_LIBS, ".final.scaffolds");
	finalize($stats, "$OUTBASE.final.scaffolds.fasta");
	print "[a5] Final assembly in $OUTBASE.final.scaffolds.fasta\n"; 
} 

# clean up any files we unzipped files if the raw input was compressed
my @unzipped = glob("$OUTBASE/*.$UNZIP_EXT");
for my $file (@unzipped){
	print STDERR "[a5] Removing unzipped file $file\n";
	`rm $file`;
}

for my $lib(keys %RAW_LIBS){
	if(defined($RAW_LIBS{$lib}{"p1_cor"})){
		`rm -f $RAW_LIBS{$lib}{"p1_cor"}`;
		`rm -f $RAW_LIBS{$lib}{"p2_cor"}`;
	}
}

#
# Begin subroutines
#
#

=item get_availmem 
Reads total installed memory in a platform-dependent way
=cut 
sub get_availmem {
	my $mem = 4000000;
	if ($^O =~ m/darwin/) {
		my $inf = `/usr/sbin/system_profiler SPHardwareDataType | grep Memory`;
		if ($inf =~ /      Memory: (\d+) GB/){
			$mem = $1 * 1048576;
		}
	} else {
		my $inf = `cat /proc/meminfo | grep MemTotal`;
		if ($inf =~ /MemTotal:\s+(\d+) kB/){
			$mem = $1;
		}
	}
	return $mem
}

sub shuffle_fastq {
	my $r1file = shift;
	my $r2file = shift;
	my $outfile = shift;
	open(R1, $r1file);
	open(R2, $r2file);
	open(OUTTER, ">$outfile");
	while(1){
		my $r1_entry;
		my $r2_entry;
		for(my $i=0; $i<4; $i++){
			my $r1l = <R1>;
			my $r2l = <R2>;
			last unless defined $r1l && $r2l;
			$r1_entry .= $r1l;
			$r2_entry .= $r2l;
		}
		last unless defined $r1_entry && $r2_entry;
		print OUTTER $r1_entry.$r2_entry;
	}
}

=item qfilter_trim_paired

Preprocesses paired reads with Trimmomatic and SGA to remove adapters, filtering low quality reads, 
do quality trimming, and separate out unpaired reads.

Takes 4 parameters: $r1file, $r2file, $t, and $lib
where 
	$r1file = read 1 in paired reads
	$r2file = read 2 in paired reads.
	$t = the number of threads to use when running SGA
	$lib = the library name to use for output generating throughout this subroutine 

=cut
sub qfilter_trim_paired {
	my $r1file = shift;
	my $r2file = shift;
	my $lib = shift;
	
	# shuffle reads
	my $p1_name = basename($r1file);
	shuffle_fastq($r1file, $r2file, "$WD/$p1_name.both.fastq");

	# run trimmomatic on the reads
	my $phred64 = get_phred64($r1file);
	my $trim_phred = $phred64 ? "-phred64 " : "-phred33 ";
	my $trim_cmd = "java ".JAVA_HEAP." -jar $DIR/trimmomatic.jar SE -threads $threads $trim_phred $WD/$p1_name.both.fastq $WD/$p1_name.trim.fastq ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36";
	print STDERR "[a5] $trim_cmd\n";
	$trim_cmd = $TIME.$trim_cmd if $debug;	
	system($trim_cmd);

	# quality filter. unsure whether to use --permute-ambiguous
	my $pp_cmd = "$DIR/sga preprocess -q ".SGA_Q_TRIM." -f ".SGA_Q_FILTER." -m ".SGA_MIN_READ_LENGTH."  --pe-mode=0 ";
	$pp_cmd .= "--phred64 " if $phred64;
	$pp_cmd .= " $WD/$p1_name.trim.fastq > $WD/$p1_name.both.pp 2> /dev/null";
	print STDERR "[a5] $pp_cmd\n";
	$pp_cmd = $TIME.$pp_cmd if $debug;	
	system($pp_cmd);

	# dust filter (use a very high dust threshold to distinguish seq adapter readthrough error from natural low complexity)
#	my $dust_cmd = "$DIR/sga preprocess --dust-threshold ".SGA_DUST." -q ".SGA_Q_TRIM." -f ".SGA_Q_FILTER." -m ".SGA_MIN_READ_LENGTH."  --pe-mode=0 ";
#	$dust_cmd .= "--phred64 " if (get_phred64($r1file));
#	$dust_cmd .= " $WD/$p1_name.both.pp > $WD/$p1_name.both.pp.dust 2> /dev/null";
#	print STDERR "[a5] $dust_cmd\n";
#	$dust_cmd = $TIME.$dust_cmd if $debug;	
#	system($dust_cmd);

}

sub correct_paired {
	my $r1file = shift;
	my $r2file = shift;
	my $lib = shift;
	
	# OS X is missing some required multithreading gadgetry for SGA
	my $t = $^O =~ m/darwin/ ? 1 : $threads;
	# shuffle reads
	my $p1_name = basename($r1file);
	# error correct
	my $ec_cmd = "$DIR/sga correct -t $t -p $OUTBASE.pp -o $WD/$p1_name.both.pp.ec.fastq $WD/$p1_name.both.pp > $WD/$lib.correct.out";
	print STDERR "[a5] $ec_cmd\n";
	$ec_cmd = $TIME.$ec_cmd if $debug;	
	system($ec_cmd);
	system("rm -rf $WD/$p1_name.both.pp"); # clear up some disk space
	
}

sub repair_fastq{
	my $input = shift;
	my $output = shift;
	my $unpair = shift;
	my $pair_r1 = shift;
	my $pair_r2 = shift;
	open(TDPIPE, $input);
	open(REPAIR, ">".$output);
	open(UNPAIR, ">".$unpair);
	open(READ1, ">$pair_r1");
	open(READ2, ">$pair_r2");
	my $lc = -1;
	my $prev_read = "";
	my $prev_record = "";
	my $faq_modulo = 4;
	my $base_count = 0;
	my $read_count = 0;
	while( my $line = <TDPIPE> ){
		$lc++;
		if($lc%$faq_modulo == 1){
		    # We are reading the sequence line
		    # -1 to not count the \n character
		    $base_count += length($line) - 1;
		    $read_count++;
		}
		if($lc%$faq_modulo != 0){
			$prev_record .= $line;
			next;
		}
		if($prev_read eq ""){
			$prev_read = $line;
			$prev_record = $line;
			next;
		}
		# we're at a header line and we have a previous FastQ record loaded
		chomp $prev_read;
		$prev_read =~ s/\/\d$//g;
		$prev_read =~ s/ .+//g;  # illumina 1.8+ tags
		my $cur_read = $line;
		chomp $cur_read;
		$cur_read =~ s/\/\d$//g;
		$cur_read =~ s/ .+//g; # illumina 1.8+ tags
		# read the rest of the current FastQ entry
		my $cur_record = $line;
		for(my $i=1; $i<$faq_modulo; $i++){
		    my $newline = <TDPIPE>;
		    $cur_record .= $newline;
		    $lc++;
		    $base_count += length($newline) - 1 if $i == 1;
		    $read_count++ if $i ==1;
		}
		# if records are paired, write them out
		# otherwise ignore the previous entry
		if($cur_read eq $prev_read){
			print REPAIR $prev_record.$cur_record;
			print READ1 $prev_record;
			print READ2 $cur_record;
			$prev_record = "";
			$prev_read = "";
		}else{
			print UNPAIR $prev_record;
			$prev_record = $cur_record;
			$prev_read = $line;
		}
	}
	close REPAIR;
	close UNPAIR;
	close READ1;
	close READ2;
	return ($base_count,$read_count);
}

=item get_phred64
Returns whether or not the given FastQ file has phred64 quality scores
=cut 
sub get_phred64{
	my $tail_file = shift;
	my $qline = `head -n 1000 $tail_file | perl -ne 'print if \$. % 4 == 0'`;
	chomp $qline;
	my $phred64 = 1;
	for my $q (split(/\n*/,$qline)){
		if (ord($q) < 59) {
			$phred64 = 0;
			last;
		}
	}
	return $phred64;
}

=item sga_clean
Clean reads using SGA. Multi-thread if running in linux environment.

Takes two parameters: an basename for output, $outbase, and a reference to a hash storing library information, $libsref.

First, error correct all reads together to prepare for building contigs using IDBA. To do so, this function passes all 
fastq files into the following command:

	sga preprocess -q SGA_Q_TRIM -f SGA_Q_FILTER -m SGA_MIN_READ_LENGTH --permute-ambiguous <fastq1> <fastq2> . . . <fastqN>

applying the --phred64 flag to sga if phred-64 quality scores are detected, and writing preprocessed reads to $outbase.pp.fastq.

$outbase.pp.fastq is indexed using the following command:

	sga index -d 0.5*AVAL_SYSTEM_MEMORY -t <n_threads> $outbase.pp.fastq


$outbase.pp.fastq is error corrected using the following command:

	sga correct -t <n_threads> -o $outbase.pp.ec.fa $outbase.pp.fastq 
 
For the indexing and correcting step, n_threads is to 1 if running in a Macintosh environment, 4 otherwise. 

Second, for each paired library in $libsref, error correct reads from each library individually for scaffolding using
the qfilter_trim_paired function.

=cut 
sub sga_clean {
	my $outbase = shift;
	my $libsref = shift;
	my %libs = %$libsref;
	# figure out which files we need to pass to SGA
	my $files = "";
	my $tail_file = "";
	for my $lib (keys %libs) {
		if (defined($libs{$lib}{"p1"})) {
			my $fq1 = $libs{$lib}{"p1"}; 
			my $fq2 = $libs{$lib}{"p2"}; 
			$tail_file = $fq1 unless length($tail_file);
			$files .= "$fq1 $fq2 ";
		}
		if (defined($libs{$lib}{"up"})){
			my $up = $libs{$lib}{"up"}; 
			$tail_file = $up unless length($tail_file);
			$files .= "$up ";
		}
	}
	print STDERR "cleaning each PE lib\n\n";
	`rm -f $WD/$outbase.pp.fastq`;
	#
	# qfilter individual paired-end libraries
	for my $lib (keys %libs) {
		print STDERR "p1 is ".$libs{$lib}{"p1"}."\n";
		next unless defined($libs{$lib}{"p1"});
		my $ec1 = $libs{$lib}{"p1"};
		my $ec2 = $libs{$lib}{"p2"};
		# capturing raw reads counts and raw base counts
		my ($reads_c, $bases_c) = get_counts( $ec1 );
		$total_reads += $reads_c;
		$total_bases += $bases_c;
		($reads_c, $bases_c) = get_counts( $ec2 );
		$total_reads += $reads_c;
		$total_bases += $bases_c;

		# clean these reads!
		print STDERR "Cleaning reads\n\n";
		qfilter_trim_paired($ec1, $ec2, $lib);
		my $p1_name = basename($ec1);
		`cat $WD/$p1_name.both.pp >> $WD/$outbase.pp.fastq`;
	}

	# build a bwt index for all of the reads
	if(!$metagenome){
		my $sga_ind = "";
		my $sga_ind_kb = $AVAILMEM/4;
		my $err_file = "$WD/index.err";
		do{
			# OS X is missing some required multithreading gadgetry for SGA
			my $t = $^O =~ m/darwin/ ? 1 : $threads;
			my $cmd = "sga index -d ".($sga_ind_kb)." -t $t $WD/$outbase.pp.fastq > $WD/index.out 2> $err_file";
			print STDERR "[a5] $cmd\n";
			$cmd = "$DIR/$cmd";
			$cmd = $TIME.$cmd if $debug;	
			$sga_ind = `$cmd`;
			$sga_ind = read_file($err_file);
			$sga_ind_kb = int($sga_ind_kb/2);
		}while(($sga_ind =~ /bad_alloc/ || $? != 0) && $sga_ind_kb > 500000);
	
		`rm -f $WD/$outbase.pp.fastq`;
	}

	#
	# error correct individual paired-end libraries
	# and construct a merged fasta
	my $merged_pe_fa = "$WD/$OUTBASE.ec.fastq";
	`rm -f $merged_pe_fa`; #ensure we are starting fresh
	for my $lib (keys %libs) {
		next unless defined($libs{$lib}{"p1"});
		my $p1_name = basename($libs{$lib}{"p1"});
		if(!$metagenome){
			correct_paired($libs{$lib}{"p1"}, $libs{$lib}{"p2"}, $lib);
		}else{
			`mv $WD/$p1_name.both.pp $WD/$p1_name.both.pp.ec.fastq`;
		}
		# re-pair the reads
		$libs{$lib}{"p1_cor"} = "$WD/$p1_name.split.r1.fq";
		$libs{$lib}{"p2_cor"} = "$WD/$p1_name.split.r2.fq";
		($EC_bases,$EC_reads)=repair_fastq("$WD/$p1_name.both.pp.ec.fastq", "$WD/$p1_name.both.repair.fastq", "$WD/$p1_name.both.unpaired.fastq", $libs{$lib}{"p1_cor"}, $libs{$lib}{"p2_cor"});
		# don't merge if this is a mate pair lib
		next if(defined($libs{$lib}{"ins"}) && $libs{$lib}{"ins"} > MAX_PAIREDEND_INSERT_MEAN);
		my $both = "$OUTBASE.s1/$p1_name.both.repair.fastq";
		print STDERR "Running cat $both >> $merged_pe_fa\n";
		`cat $both >> $merged_pe_fa`;
		if( -f "$OUTBASE.s1/$p1_name.both.unpaired.fastq" ){
			`cat $OUTBASE.s1/$p1_name.both.unpaired.fastq >> $WD/$OUTBASE.merged.clean.unpaired.fa`;
		}
	}
	print STDERR "Done merging libraries\n";
	die "[a5] ERROR: Unable to identify any properly paired reads\n[a5] ERROR: Please check that you have provided at least one library of paired-end reads with names conforming to the Illumina read pair convention.\n" if -s $merged_pe_fa < 5;
	
	die "[a5] Error correcting reads with SGA\n" if( $? != 0 );
	`rm $OUTBASE.pp.*`;	# clean up temporary files used in error correction 
	return $merged_pe_fa;
}

=item get_counts
 
 Reads through a fastq file and returns the number of reads and the number of bases
 
 Input : 1 fastQ files 
 
=cut
 
sub get_counts {
    my $file = shift;
    open(INFILE, $file);
    my ($r_count, $b_count) = 0;
    my @read = ();
    while(<INFILE>){
	$read[0] = $_;
	$read[1] = <INFILE>;
	$read[2] = <INFILE>;
	$read[3] = <INFILE>;
	$r_count++;
	$b_count += length($read[1])-1;
	@read=();
    }
    close(INFILE);
    return ($r_count, $b_count);
}

=item check_and_unzip
Check if the given file is zipped up. If it is, unzip and return the file it was unzipped to.
=cut
sub check_and_unzip {
	my $file = shift;
	my $file_type = `file $file`;
	my $ret = $file;
	validate_sequence_file($file, $message);
	my $bname = basename($file);
	if ($file_type =~ /gzip/){
		$ret = "$OUTBASE/$bname.$UNZIP_EXT";
		`mkdir -p $OUTBASE ; gunzip -c $file > $ret`;
	} elsif ($file_type =~ /bzip2/) {
		$ret = "$OUTBASE/$bname.$UNZIP_EXT";
		`mkdir -p $OUTBASE ; bunzip2 -c $file > $ret`;
	}
	return $ret;
}

=item map_libs_and_polish
Maps all libraries back to the final scaffolds generated by the A5 pipeline. 
Takes three parameters: $outbase, $libsref, and $suffix

where 
	$outbase = the basename for output files generated throughout this subroutine
	$libsref = a hash storing information on all the libraries passed in.
	$suffix = usually either .contigs or .final.scaffolds to indicate which assembly to use
=cut
sub map_libs_and_polish {
	my $outbase = shift;
	my $libsref = shift;
	my $suffix = shift;  # either .contigs or .final.scaffolds
	my %libs = %$libsref;
	my $index_cmd = "bwa index $OUTBASE$suffix.fasta";
	print STDERR "[a5] $index_cmd\n";
	$index_cmd = "$DIR/$index_cmd";
	$index_cmd = $TIME.$index_cmd if $debug;	
	system($index_cmd);
	my %coverage;
	my $lib_list = "";
	for my $lib (keys %libs) {
		if (defined($libs{$lib}{"p1"})) {
			my $seqs = $libs{$lib}{"p1"}." ".$libs{$lib}{"p2"};
			map_lib("$lib.pe", $seqs, \%coverage, $suffix);
			# report insert size
			open(VIEW, "$DIR/samtools view $OUTBASE.$lib.pe.sort.bam |");
			$lib_list .= " $OUTBASE.$lib.pe.sort.bam";
			my $count = 0;
			my $ins = 0;
			my @sdlist;
			while(my $line = <VIEW>){
				my @d = split(/\t/, $line);
				$ins += $d[8] if $d[8] > 0;
				$count++ if $d[8] > 0;
				$sdlist[$count % 1000] = $d[8] if $d[8] > 0;
			}
			my $sd = int(stddev($ins/$count, \@sdlist));
			print "[a5] Library $lib avg. insert size ".int($ins/$count)." +/- $sd\n";
		}
		if (defined($libs{$lib}{"up"})){
			map_lib("$lib.up", $libs{$lib}{"up"}, \%coverage, $suffix);
		}
	}
	# "polish" the assembly
	`$DIR/samtools faidx $OUTBASE$suffix.fasta`;
	`$DIR/samtools mpileup -uf $OUTBASE$suffix.fasta $lib_list | $DIR/bcftools view -cg - | $DIR/vcfutils.pl vcf2fq -d 1 > $OUTBASE$suffix.fastq`;
	# convert fastq to qvl file
	my $above_q40 = fastq_to_qvl("$OUTBASE$suffix.fastq", "$OUTBASE$suffix.qvl");
	my $corrected = correct_basecalls("$OUTBASE$suffix.fastq", "$OUTBASE$suffix.fasta");
	`mv $corrected $OUTBASE$suffix.fasta`;
	# calculate median cov
	my @c;
	foreach my $scaf(keys(%coverage)){
		foreach my $i(keys(%{$coverage{$scaf}})){
			push(@c, $coverage{$scaf}{$i});
		}
	}
	my @cc = sort {$a <=> $b} @c;
	print "[a5] Median depth of coverage: ".$cc[@cc/2].", 10th percentile: ".$cc[@cc/10]."\n";
	print "[a5] $above_q40 / ".scalar(@cc)." bases at or above Q40\n";
	# clean up temporary files used for read mapping
	`rm $OUTBASE$suffix.fasta.*`;
	return {median_cov => $cc[@cc/2], cov_10 => $cc[@cc/10], above_q40 => $above_q40};
}

=item correct_basecalls

Use a FastQ file containing revised base calls (e.g. as produced by bcftools/vcfutils)
to correct the base calls present in the FastA scaffold assembly. Where differences exist
between the FastA and FastQ, the basecall in the FastA is replaced with that contained in 
the FastQ. This is done in all cases except where bcftools has changed the sequence 
length, possibly due to some bug in that software. 

=cut
sub correct_basecalls {
	my $fastq = shift;
	my $fasta = shift;
	open(FASTA, $fasta);
	open(FASTQ, $fastq);
	my @fastq_lines;
	my $in_sequence = 0;
	while(my $line = <FASTQ>){
		if($line =~ /^\@scaf/){
			$line =~ s/^@/>/g;
			push(@fastq_lines, $line);
			$in_sequence = 1;
		}elsif($line =~ /^\+/){
			$fastq_lines[-1] .= "\n";
			$in_sequence = 0;
		}elsif($in_sequence > 0){
			chomp $line;
			push(@fastq_lines, $line) if $in_sequence == 1;
			$fastq_lines[-1] .= $line if $in_sequence > 1; # append to last element
			$in_sequence++;
		}
	}

	my @fasta_lines = <FASTA>;
	my $corrected_bases = 0;
	for(my $lineI=0; $lineI < @fastq_lines; $lineI++){
		next if $lineI % 2 == 0; # skip header lines (even numbered)
		my @fa_line = split( //, $fasta_lines[$lineI] );
		my @fq_line = split( //, $fastq_lines[$lineI] );
		next if @fa_line != @fq_line;
		print STDERR "ERROR: different lengths between consensus and IDBA-UD assembly!\nSequence $fasta_lines[$lineI-1], lengths ".scalar(@fa_line)." and ".scalar(@fq_line)."\n.qvl file may be unusable." if @fa_line != @fq_line;
		for(my $i=0; $i<@fq_line; $i++){
			if($fq_line[$i] ne $fa_line[$i] && $fq_line[$i] ne "n" && $fq_line[$i] ne "N"){
				$fa_line[$i] = $fq_line[$i];
				$corrected_bases++;
			}
		}
		$fasta_lines[$lineI] = join("", @fa_line);
	}
	close FASTA;
	open(FASTA, ">$fasta.corrected.fa");
	print FASTA @fasta_lines;
	print "[a5] Corrected $corrected_bases base(s) in consensus sequence\n";
	return "$fasta.corrected.fa";
}

=item fastq_to_qvl

Create a qvl file using the quality scores contained in a phred33 encoded fastq file

=cut
sub fastq_to_qvl {
	my $fastq = shift;
	my $qvl = shift;
	open(FQ, $fastq);
	open(QVL, ">$qvl");
	my $in_sequence = 0;
	my $above_q40 = 0;
	while(my $line = <FQ>){
		chomp $line;
		if($line =~ /^@(scaff.+)/){
			print QVL "\n" if $in_sequence == 2;
			print QVL ">$1\n";
			$in_sequence = 1;
		}elsif($line =~ /^\+/ && $in_sequence == 1){
			$in_sequence = 2;
		}elsif($in_sequence == 2){
			my @qual = split( //, $line );
			foreach my $q(@qual){
				$above_q40++ if (ord($q) - 33) >= 40;
				print QVL (ord($q) - 33);
				print QVL " ";
			}
		}
	}
	return $above_q40;
}

=item map_lib

Uses bwa mem to map reads from a library back to the assembly

=cut
sub map_lib {
	my $lib = shift;
	my $seqs = shift;
	my $coverage = shift;
	my $suffix = shift;
	my $sampe_cmd = "bwa mem -t $threads $OUTBASE$suffix.fasta $seqs 2> /dev/null | $DIR/samtools view -b -S - | $DIR/samtools sort - $OUTBASE.$lib.sort";
	my $bamindex_cmd = "samtools index $OUTBASE.$lib.sort.bam";
	print STDERR "[a5] $sampe_cmd\n";
	$sampe_cmd = "$DIR/$sampe_cmd";
	$sampe_cmd = $TIME.$sampe_cmd if $debug;	
	system($sampe_cmd);
	print STDERR "[a5] $bamindex_cmd\n";
	$bamindex_cmd = "$DIR/$bamindex_cmd";
	$bamindex_cmd = $TIME.$bamindex_cmd if $debug;	
	system($bamindex_cmd);
	open(MPILEUP, "$DIR/samtools mpileup -d1000000 -BQ0 $OUTBASE.$lib.sort.bam |");
	my $i=0;
	while(my $line = <MPILEUP>){
		my @d = split(/\t/, $line);
		$coverage->{$d[0]}{$d[1]} = 0 unless defined $coverage->{$d[0]}{$d[1]};
		$coverage->{$d[0]}{$d[1]} += $d[3];
	}	
}

=item stddev

Calculates standard deviation of a collection of numbers

=cut
sub stddev{
	my $average = shift;
	my($data) = @_;
	return 0 if(@$data == 1);
	my $sqtotal = 0;
	foreach(@$data) {
		$sqtotal += ($average-$_) ** 2;
	}
	my $std = ($sqtotal / (@$data-1)) ** 0.5;
	return $std;
}

=item read_lib_file

Parses the library file, returning the given information stored in a hash. 
Takes a singled parameter: $libfile, the path to the library file to be parsed.

Returns a hash of hashes index as such: $hash{$library}{$component}
where $component = 'up' || 'p1' || 'p2' || 'ins' || 'id'

=cut
sub read_lib_file {
	my $libfile = shift;
	my %libs = ();
	open(LIB,"<",$libfile);
	my $lib_count = 0;
	my $id = "";
	my %hash = ();
	while(<LIB>){
		chomp;
		if ($_ =~ m/\[LIB\]/){
			if ($lib_count > 0) {
				$libs{$id}{"id"} = $id;
				for my $key (keys %hash){
					$libs{$id}{$key} = $hash{$key};
					delete($hash{$key});
				}
			} 
			$lib_count++;
			$id = "raw$lib_count";
		} elsif ($_ =~ m/id=([\w\/\-\.]+)/) { 
			$id = $1;
		} elsif ($_ =~ m/(\w+?)=(.+)/) { 
			$hash{$1} = $2;
		} else {
			die "[a5] Unrecognizable line in library file: >$_<\n";
		}
	} 
	$libs{$id}{"id"} = $id;
	for my $key (keys %hash){
		$libs{$id}{$key} = $hash{$key};
	}
	return %libs;
}

=item read_file
Read in entire file, returning its contents as a single string.
=cut
sub read_file {
	my $file = shift;
	local $/=undef;
	open(FILE,"<",$file);
	my $str = <FILE>;
	return $str;
}

=item fastq_to_fasta
Convert a fastq file to fasta.
=cut
sub fastq_to_fasta {
	my $fastq = shift;
	my $fasta = shift;
	my $maxrdlen = 0;
	open(FQ,"<$fastq");
	open(FA,">$fasta");
	while(my $hdr = <FQ>) {
		my $seq = <FQ>;
		my $qhdr = <FQ>;
		my $qseq = <FQ>;
		$hdr =~ s/^@/>/g;
		print FA $hdr.$seq;
		$maxrdlen = length($seq) if length($seq) > $maxrdlen;
	}
	close FA;
	# why are we removing 2 here?
	return ($fasta, $maxrdlen - 2); # -1 removes newline char
}

=item split_shuf
Splits a shuffled (a.k.a. interleaved) fastq file into two separate fastq files. 

Takes two parameters: $shuf, $outbase
where
	$shuf = the shuffled fastq file to split
	$outbase = the basename for the two fastq files resulting from splitting $shuf. The two files output
	           are $outbase_p1.fastq and $outbase_p2.fastq


=cut
sub split_shuf {
	my $shuf = shift;
	my $outbase = shift;
	my $fq1 = "$outbase\_p1.fastq";
	my $fq2 = "$outbase\_p2.fastq";
	print STDERR "[a5] Splitting shuffled file $shuf into $fq1 and $fq2\n"; 
	open(FQ1,">$fq1");
	open(FQ2,">$fq2");
	open(IN,"<$shuf");
	while(my $hdr = <IN>){
		my $seq = <IN>;
		my $qhdr = <IN>;
		my $qual = <IN>;
		print FQ1 "$hdr"."$seq"."$qhdr"."$qual";
		$hdr = <IN>;
		$seq = <IN>;
		$qhdr = <IN>;
		$qual = <IN>;
		print FQ2 "$hdr"."$seq"."$qhdr"."$qual";
	}
	close FQ1;
	close FQ2;
	return ($fq1, $fq2);
}


=item idba_pe_assemble
Assemble reads using IDBA. 
Takes 3 parameters: $outbase, $reads, $maxrdlen
where 
	$outbase = the basename for output generated throughout this subroutine
	$reads = the path to the file containing the reads to assemble into contigs
	$maxrdlen = the maximum kmer to iterate up to. 
expects a file called $outbase.clean.fa in the current working directory
=cut
sub idba_pe_assemble {
	my $outbase = shift;
	my $libs = shift;
	my $maxrdlen = shift;
	
	my $merged_pe_fa = "$WD/$outbase.ec.fasta";
	my $unpaired_option = "";
	$unpaired_option = " -l $WD/$outbase.merged.clean.unpaired.fa " if -f "$WD/$outbase.merged.clean.unpaired.fa";
	my $idba_bin = "idba_ud";
	$idba_bin = "idba_ud500" if $maxrdlen > 120;
	die "[a5] Reads are too long to be assembled with A5. This version of A5 supports read lengths up to 500nt." if $maxrdlen > 500;
	my $idba_cmd = "$DIR/$idba_bin --num_threads $threads -r $merged_pe_fa $unpaired_option -o $WD/$outbase --mink ".IDBA_MIN_K." --maxk $maxrdlen --min_pairs 2 ";
	my $pe_size = -s $merged_pe_fa;
	print STDERR "[a5] $merged_pe_fa has $pe_size bytes of FastA sequence data\n";
	print STDERR "[a5] $idba_cmd\n";
	`$idba_cmd > $WD/idba.out`;
	# idba_ud has a race condition that appears to crop up only on Macs when running multi-threaded.
	# if we're on a Mac try running it single-threaded in the hope that it works.
	# if it still fails then there may be a problem with the sequence input, or another bug.
	if($? != 0 && $^O =~ m/darwin/){
		$idba_cmd .= " --num_threads 1 ";
		print STDERR "[a5] Previous IDBA-UD run failed. Trying again single-threaded.\n";
		`$idba_cmd > $WD/idba.out`;
	}
	if ($? != 0){
		# idba_ud has a bug in estimating the insert size which causes a crash before scaffolding.
		# when that happens, the contigs can still be used. If contigs were not generated then
		# idba failed early and there is no point continuing.
		die "[a5] Error building contigs with IDBA\n" unless -e "$WD/$outbase/contig.fa";
		`cp $WD/$outbase/contig.fa $WD/$outbase-contig.fa`;		
	}else{
		`cp $WD/$outbase/scaffold.fa $WD/$outbase-contig.fa`;
	}
	`rm -rf $WD/$outbase`; # free up disk space (the idba_ud output can be very large)
	return "$WD/$outbase-contig.fa";
}


=item scaffold_sspace
Scaffold contigs using SSPACE. 

Takes 3 or 4 parameters: $outbase, $libsref, $curr_ctgs, $rescaffold (optional)
where
	$outbase = the basename for output generated throughout this subroutine
	$libsref = a hash storing information on all the libraries passed in.
	$curr_ctgs = the contigs to scaffold
	$rescaffold = a binary value indicating whether or not we are scaffolding our assembly 
                  for the first time (i.e. scaffolding contigs generated from IDBA), indicated 
                  by a 0, or if we are scaffolding for the second time (i.e. scaffolding 
                  scaffolds broken by the A5qc algorithm), indicated by a 1. 

Returns the contigs resulting from scaffolding.

=cut 
sub scaffold_sspace {
	my $outbase = shift;
	# build library file
	my $libsref = shift;
	my %libs = %$libsref;
	my $curr_ctgs = shift;
	my $rescaffold = shift;

	my @library_file = ();
	my @lib_files;

	my $genome_size = get_genome_size($curr_ctgs);
	print STDERR "[a5] Total contig length $genome_size\n";
	my $libraryI=1;
	my @curr_lib_file = ();
	my $curr_ins = -1;
	my %run_lib;
	$rescaffold=0 unless defined($rescaffold);
	# sort library.txt to so that we scaffold with smaller insert libraries first
	for my $lib (sort { $libs{$a}{"ins"} <=> $libs{$b}{"ins"} } keys %libs) {
		my ($exp_link, $cov) = calc_explinks( $genome_size, $libs{$lib}{"ins"}, $libs{$lib}{"p1"} ); 
		printf STDERR "[a5] %s\: Insert %.0f, coverage %.2f, expected links %.0f\n", $libs{$lib}{"id"}, $libs{$lib}{"ins"}, $cov, $exp_link;
		if (-f "$WD/$OUTBASE.ec.fa") { # run sspace with unpaired library if we have one
			$curr_ctgs = run_sspace($genome_size, $libs{$lib}{"ins"}, $exp_link, $outbase.".".$libs{$lib}{"id"},
                                                  $libs{$lib}{"libfile"}, $curr_ctgs, $rescaffold, "$WD/$OUTBASE.ec.fa");
		} else {
			$curr_ctgs = run_sspace($genome_size, $libs{$lib}{"ins"}, $exp_link, $outbase.".".$libs{$lib}{"id"},
                                                  $libs{$lib}{"libfile"}, $curr_ctgs, $rescaffold);
		}
	}

	return $curr_ctgs;
}

=item preprocess_libs
Process a library hash generated with the read_lib_file subroutine. For each library, check to see
if an error corrected version exists, and swap that into the hash if so. Then, do initial insert size
estimates for each library. After estimating insert sizes, merge libraries with similar insert size distributions
using the subroutine aggregate_libs. Insert size distributions are deemed similiar if their calculated ranges overlap. 
If two libraries are merged, the insert size of the resulting merged libarary is recalculated. Lastly, for each library,
create a library file for SSPACE>

Takes 2 parameters: $libsref, $ctgs
where
	$libsref = a hash storing information on all the libraries passed in.
	$ctgs = contigs to use to estimate insert sizes

Returns a new library hash containing insert size calculations and new library information for merged libraries.
=cut
sub preprocess_libs {
	my $libsref = shift;
	my %libs = %$libsref;
	my $ctgs = shift;
	if (-f "$OUTBASE.unpaired.fastq"){
		print STDERR "[a5] Removing unpaired reads $OUTBASE.unpaired.fastq\n";
		`rm $OUTBASE.unpaired.fastq`;
	}
	# for each lib, check whether an error corrected version exists and use that instead if possible
	for my $lib (keys %libs) {
		next unless defined($libs{$lib}{"p1_cor"});
		if(!-e $libs{$lib}{"p1_cor"}){
			die "Trimmed, error corrected reads unavailable for library $lib\nShould be located at ".$libs{$lib}{"p1_cor"}."\nUnable to proceed with scaffolding\n";
		}
		$libs{$lib}{"p1"} = $libs{$lib}{"p1_cor"};
		$libs{$lib}{"p2"} = $libs{$lib}{"p2_cor"};
	}

#
#       MAKE OUR INITIAL ESTIMATES OF INSERT SIZE HERE AND GET ERROR ESTIMATES.
#       ALSO CONCATENATE UNPAIRED LIBRARIES AND REMOVE FROM HASH
#
	my $have_up = 0;
	print STDERR "[a5] Making initial estimates of insert size\n";
	for my $libid (sort { $libs{$a}{"id"} cmp $libs{$b}{"id"} }keys %libs) {
		if (defined($libs{$libid}{"p1"})) {
			my $fq1 = $libs{$libid}{"p1"};
			my $fq2 = $libs{$libid}{"p2"};
			my ($ins_mean, $ins_err, $outtie) = get_insert($fq1,$fq2,"$OUTBASE.$libid",$ctgs);
			$libs{$libid}{"rc"} = $outtie;
			if ($ins_mean > 0) {
				$libs{$libid}{"ins"} = $ins_mean;
				$libs{$libid}{"err"} = $ins_err;
				$libs{$libid}{"rc"} = $outtie;
			} else {
				print STDERR "[a5_preproc] Insert estimate for $libid not reliable.";
				if (defined($libs{$libid}{"ins"})) {
					print STDERR "[a5_preproc] Using given insert estimate instead.\n";
					$libs{$libid}{"err"} = 0.95;
				} else {
					print STDERR "[a5_preproc] No insert estimate for $libid. Please provide one in a library file. Exiting\n";
					exit -1;
				}
		
			}
			if (defined($libs{$libid}{"up"})){
				my $up = $libs{$libid}{"up"};
				`cat $up >> $OUTBASE.unpaired.fastq`;
				 # don't aggregate unpaired along with paired reads, just dump them all into one file 
				delete($libs{$libid}{"up"});
				$have_up = 1;
			}
		} elsif (defined($libs{$libid}{"up"})){
			my $up = $libs{$libid}{"up"};
			`cat $up >> $OUTBASE.unpaired.fastq`;
			 # isn't paired, don't include in aggregated libraries because we'll just dump these into one file
			delete($libs{$libid});
			$have_up = 1;
		}
	}
#
#       NOW MERGE SIMILIAR LIBRARIES
#
	my @curr_lib_file = ();
	my $curr_lib = "";
	my $prev_lib;
	my $libraryI = 1;
	my %processed = ();
	my $run_lib;
	print STDERR "[a5] Will merge libraries if similar enough\n";
	# sort libraries so that we scaffold with smaller insert libraries first
	for my $lib (sort { $libs{$a}{"ins"} <=> $libs{$b}{"ins"} } keys %libs) {
		if (defined($prev_lib)) {
			my $curr_ins = $libs{$lib}{"ins"};
			my $prev_ins = $libs{$prev_lib}{"ins"};
			my $curr_min = $libs{$lib}{"ins"}*(1-$libs{$lib}{"err"});	
			my $curr_max = $libs{$lib}{"ins"}*(1+$libs{$lib}{"err"});	
			my $prev_min = $libs{$prev_lib}{"ins"}*(1-$libs{$prev_lib}{"err"});	
			my $prev_max = $libs{$prev_lib}{"ins"}*(1+$libs{$prev_lib}{"err"});	
			#print STDERR "[a5] \$prev_ins = $prev_ins, \$curr_ins = $curr_ins\n";
			
			# if we've hit a substantially different insert size (i.e. min and max 
			# inserts don't overlap or currrent insert is twice the size of the 
			# previous insert), combine are current aggregate, and move one with
			# the current library
			if (!($curr_min <= $prev_max && $prev_min <= $curr_max) || $curr_ins/$prev_ins > 2){
				# scaffold with the previous insert....
				# combine libraries if necessary, and return just one library hash
				$run_lib = aggregate_libs(\@curr_lib_file,$curr_lib,$ctgs);
				$run_lib->{"libfile"} = print_libfile("$OUTBASE.library_$libraryI.txt", $run_lib, $run_lib->{"err"});
				print_libfile("$OUTBASE.library_$libraryI.txt.strict", $run_lib, $run_lib->{"err"});
				for my $key (keys %$run_lib){
					$processed{$run_lib->{"id"}}{$key} = $run_lib->{$key};
				}
				# now move on to the next library...
				$libraryI++;
				@curr_lib_file = ();
				$curr_lib = "";
			}
		}

		push (@curr_lib_file,$libs{$lib});
		$curr_lib .= $lib;
		$prev_lib = $lib;
	}
	$run_lib = aggregate_libs(\@curr_lib_file,$curr_lib,$ctgs);
	$run_lib->{"libfile"} = print_libfile("$OUTBASE.library_$libraryI.txt", $run_lib, $run_lib->{"err"});
	print_libfile("$OUTBASE.library_$libraryI.txt.strict", $run_lib, $run_lib->{"err"});
	open(LIB,">$OUTBASE.preproc.libs");
	print STDERR "[a5] Printing preprocessed library file to $OUTBASE.preproc.libs\n";
	for my $key (keys %$run_lib){
		$processed{$run_lib->{"id"}}{$key} = $run_lib->{$key};
	}
	for my $lib (sort keys %processed){
		print LIB "[LIB]\n";
		for my $att (@KEYS) {
			next if !defined($processed{$lib}{$att});
			print LIB "$att=".$processed{$lib}{$att}."\n";
		}
	}
	print LIB "[LIB]\nid=$def_up_id\nup=$OUTBASE.unpaired.fastq\n" if $have_up;
	close LIB;
	return %processed;
}

=item aggregate_libs
Recreate the library hash for the given library/libraries. Merge files if an aggregate of libraries is passed in.
Calculate the insert size of the final library. 

Takes 3 parameters: $curr_lib_file, $curr_lib, $curr_ctgs
where 
	$curr_lib_file = an array of single-library hashes. If the array has multiple entries (i.e. if 
	                 an aggregate has been passed in) the corresponding libraries are merged using the
	                 subroutine merge_libraries
	$curr_lib = the id of the library being processed here
	$curr_ctgs = contigs to use for insert size calculations
=cut
sub aggregate_libs {
	my $curr_lib_file = shift;
	my $curr_lib = shift;
	my $curr_ctgs = shift;
	my ($fq1, $fq2,$up);
	my %fin_lib = ();
	print STDERR "[a5] aggregating libraries. n = ".scalar(@$curr_lib_file)."\n";
	if (scalar(@$curr_lib_file) > 1) { # if this is an aggregate of libraries, combine them into one file
		($fq1, $fq2, $up) = merge_libraries($curr_lib_file,$OUTBASE.".".$curr_lib);
	} else {
		($fq1, $fq2) = ($curr_lib_file->[0]{"p1"}, $curr_lib_file->[0]{"p2"});
		$up = $curr_lib_file->[0]{"up"} if (defined($curr_lib_file->[0]{"up"}));
	}

	$fin_lib{"id"} = $curr_lib;
	$fin_lib{"p1"} = $fq1;
	$fin_lib{"p2"} = $fq2;
	$fin_lib{"rc"} = 0;
	# re-estimate insert sizes
	my ($ins_mean, $ins_err, $outtie) = get_insert($fq1,$fq2,"$OUTBASE.$curr_lib",$curr_ctgs);
	$fin_lib{"ins"} = abs($ins_mean);
	$fin_lib{"err"} = abs($ins_err);	
	$fin_lib{"rc"} = $outtie;
	$fin_lib{"nlibs"} = scalar(@$curr_lib_file);
	return \%fin_lib;		
}

=item print_libfile
Print an SSPACE library file. 

Takes 3 parameters: $file, $libref, $err_estimate
where
	$file = the file to write the SSPACE library file to
	$libref = library hash to use to create the SSPACE library file.
	$err_estimate = the error estimate for this library to give to SSPACE
=cut
sub print_libfile {
	my $file = shift;
	my $libref = shift;
	my $err_estimate = shift;
	open( LIBRARY, ">$file" );
	print LIBRARY $libref->{"id"}." ".$libref->{"p1"}." ".$libref->{"p2"}." ".
                  $libref->{"ins"}." $err_estimate ".$libref->{"rc"}."\n";
	close LIBRARY;
	return $file;
}

=item print_lib_info
Print library information from the given library hash to standard error. 
=cut
sub print_lib_info {
	my $LIBS = shift;
	for my $lib (sort keys %$LIBS) {
		print STDERR "     $lib:\n";
		for my $att (@KEYS){
			#print STDERR "$att is undefined for $lib\n" if !defined($LIBS->{$lib}{$att});
			next if !defined($LIBS->{$lib}{$att});
			print STDERR "      $att=".$LIBS->{$lib}{$att}."\n";
		}
	}
}
=item merge_libraries
Merge library files corresponding to the libraries given by the array of library hashes. 

Takes 2 parameters: $curr_lib_file, $curr_lib
where
	$curr_lib_file = an array of single-library hashes. 
	$curr_lib = the id of the library being processed here

Returns an array with the resulting files.
=cut
sub merge_libraries {
	my $curr_lib_file = shift; # array of hashes
	my $curr_lib = shift;
	my $fq1 = "$curr_lib\_p1.fastq";
	my $fq2 = "$curr_lib\_p2.fastq";
	my $up = "$curr_lib\_up.fastq";
	open(my $fq1h,">","$fq1");
	open(my $fq2h,">","$fq2");
	my $uph;
	# merge files for scaffolding, reverse complementing as necessary
	print STDERR "[a5] Merging libraries $curr_lib\n";
	for my $sublib (@$curr_lib_file){
		my @ar  = split(' ',$sublib);
	#	print STDERR " ".$sublib->{"id"};
		print STDERR "[a5] Piping ".$sublib->{"p1"}."\n"; 
		open(my $fq,"<",$sublib->{"p1"});
		pipe_fastq($fq,$fq1h,$sublib->{"rc"});
		print STDERR "[a5] Piping ".$sublib->{"p2"}."\n"; 
		open($fq,"<",$sublib->{"p2"});
		pipe_fastq($fq,$fq2h,$sublib->{"rc"});
		if (defined($sublib->{"up"})){  # build an unpaired file if there are unpaired reads
			if (defined($uph)){
				open($uph,">","$up");
			}
			open($fq,"<",$sublib->{"up"});
			pipe_fastq($fq,$uph,$sublib->{"rc"});
		}
	}	
	close $fq1h;
	close $fq2h;
	close $uph if defined($uph);
	return ($fq1, $fq2, $up);
}
=item break_all_misasms
Using the A5QC algorithm, remove any misassembled regions from the given contigs using the given libraries. 

Takes 3 parameters: $ctgs, $libsref, $outbase
where
	$ctgs = the contigs to remove misassemblies from
	$libsref = a hash storing information on all the libraries passed in.
	$outbase = the basename for output generated throughout this subroutine

Return the contigs resulting from removing misassemblies.
=cut
sub break_all_misasms {
	my $ctgs = shift;
	my $libsref = shift;
	my %libs = %$libsref;
	my $outbase = shift;
	my @lib_files;
	my $genome_size = get_genome_size($ctgs);
	my $broken = 0;
	# sort libraries by insert so we break with larger inserts first.
	for my $lib (sort { $libs{$b}{"ins"} <=> $libs{$a}{"ins"} } keys %libs) {
		#AED: only use large inserts for now, or the largest of the small insert libraries
		next unless ($libs{$lib}{"ins"} > 2000 || $broken==0);
		$broken = 1;
		my $ins = $libs{$lib}{"ins"};
		my $fq1 = $libs{$lib}{"p1"};
		my $fq2 = $libs{$lib}{"p2"};
		my ($exp_link, $cov) = calc_explinks( $genome_size, $libs{$lib}{"ins"}, $libs{$lib}{"p1"} ); 
		my $min_len = int($libs{$lib}{"ins"}*(1-$libs{$lib}{"err"}));	
		my $max_len = int($libs{$lib}{"ins"}*(1+$libs{$lib}{"err"}));	
		my $sspace_k = int(log($exp_link)/log(1.4)-9.5);
		my $nlibs = $libs{$lib}{"nlibs"};
		if ($libs{$lib}{"ins"} > MAX_PAIREDEND_INSERT_MEAN) {
			$nlibs++;
			print STDERR "[a5_break_misasm] Expecting shadow library in library $lib\n"; 
		}
		my $prev_ctgs = $ctgs;
		$ctgs = break_misasms($ctgs,$fq1,$fq2,"$outbase.lib$lib",$nlibs);
	}	
	return $ctgs;
}

=item break_misasms
Remove any misassmblies from the given contigs using the given paired library. 

Takes 5 parameters: $ctgs, $fq1, $fq2, $outbase, $nlibs
where
	$ctgs = the contigs to remove misassemblies from
	$fq1 = read 1 in the paired library
	$fq2 = read 2 in the paired library
	$outbase = the basename for output generated throughout this subroutine
	$nlibs = the number of libraries in the paired library. Use to indicate if this
	         library is the result of merging similiar libraries or if there is an 
	         expected shadow library.
=cut
sub break_misasms {
	my $ctgs = shift;
	my $fq1 = shift;
	my $fq2 = shift;
	my $outbase = shift;
	my $nlibs = shift;
	print STDERR "[a5] Identifying misassemblies in $ctgs with $outbase\n";
	my $sai = "$WD/$outbase.sai";
	my $sam = "$WD/$outbase.sam";
	`$DIR/bwa index $ctgs > $WD/$outbase.index.out`;
	`cat $fq1 $fq2 | $DIR/bwa mem -M -t $threads $ctgs - > $sam`;
	`$DIR/samtools view -bhS $sam | $DIR/samtools sort -o -n - $outbase.sort | $DIR/samtools view -F 0x100 -h - > $outbase.sorted.sam`;
	`mv $outbase.sorted.sam $sam`;
	`rm $outbase.sort*` if -f "$outbase.sort*";
	`rm $ctgs.*`;	# remove temp files required by mapper
	# Let Java use 2/3 of available memory
	my $mem = (int((($AVAILMEM)*0.66)/1024))."m";
	my $cmd = "A5qc.jar $sam $ctgs $WD/$outbase.broken.fasta $nlibs > $WD/$outbase.qc.out";
	print STDERR "[a5] java -Xmx$mem -jar $cmd\n"; 
	if ($debug) {
		`$TIME java -Xmx$mem -jar $DIR/$cmd`;
	} else {
		`java -Xmx$mem -jar $DIR/$cmd`;
	}
	die "[a5] Error in detecting misassemblies.\n" if ($? != 0);
	`rm -f $sam`;
	my $qc_stdout = read_file("$WD/$outbase.qc.out");
	if ($qc_stdout =~ /No blocks were found/){
		return $ctgs;
	} else {
		return "$WD/$outbase.broken.fasta";
	}
}

=item pipe_fastq
Write fastq data from an input stream to an output stream, reverse completmenting 
fastq entries if indicated.

Takes 3 parameters: $from, $to, $rc
where
	$from = the input stream to read from
	$to = the output stream to write to
	$rc = a binary value indicating whethere or not fastq entries should be reverse complemented
=cut
sub pipe_fastq {
    my $from = shift;
    my $to = shift;
    my $rc = shift;
    while (!eof($from)){
        my $line = readline($from);
        print $to $line;
        $line = readline($from);
        if ($rc){
            chomp $line;
            $line =~ tr/ACGTacgt/TGCAtgca/; # complement
            $line = reverse($line);         # reverse
            print $to $line."\n";
			$line = readline($from);
			chomp $line;
            print $to $line."\n";
            $line = readline($from);
            chomp $line;
            $line = reverse($line);         # reverse
            print $to $line."\n";
        } else {
			chomp $line;
            print $to $line."\n";
			$line = readline($from);
			chomp $line;
            print $to $line."\n";
			$line = readline($from);
			chomp $line;
            print $to $line."\n";
        }   
    }   
}
=item calc_explinks
Calculate the expected number of read pairs to span a point in the assembly.
=cut
sub calc_explinks {
	my $genome_size = shift;
	my $ins_len = shift;
	my $read_file = shift;
	my $read_count = 0;
	my $maxrdlen = -1;
	open( READFILE, $read_file );
	while( my $line = <READFILE> ){
		$read_count++;
		$maxrdlen = length($line) if $read_count % 4 == 2 && $maxrdlen < length($line);
	}
	$read_count /= 4;
	my $cov = $maxrdlen * $read_count / $genome_size;
	my $exp_link = $cov * $ins_len / $maxrdlen;
	return ($exp_link, $cov);
}

=item run_sspace
Run SSPACE using the given parameters and data.
Takes 7 parametes: $genome_size $insert_size $exp_links $outbase $libfile $input_fa $rescaffold
	$genome_size = the expected size of the genome (used to calculate SSPACE parameters)
	$insert_size = the insert size of the library (used to calculate SSPACE parameters)
	$exp_links = the expected links to span a point in the contigs we are scaffold
	$outbase = the basename for output generated throughout this subroutine
	$libfile = the SSPACE library file to give SSPACE
	$input_fa = the contigs to scaffold
	$rescaffold = a binary value indicating whether or not we are scaffolding our assembly 
                  for the first time (i.e. scaffolding contigs generated from IDBA), indicated 
                  by a 0, or if we are scaffolding for the second time (i.e. scaffolding 
                  scaffolds broken by the A5qc algorithm), indicated by a 1. 

Return the scaffolds resulting from running SSPACE
=cut
sub run_sspace {
	my $genome_size = shift;
	my $insert_size = shift;
	my $exp_links = shift;
	my $outbase = shift;
	my $libfile = shift;
	my $input_fa = shift;
	my $rescaffold = shift;

	my $sspace_x = 0;
	my $sspace_m = int(log2($genome_size)+3.99);
	$sspace_m = 15 if $sspace_m < 15; # min overlap of read during extension. as genome size goes down, k-mers of a particular size are more likely to be unique, thus overlaps can be trusted at smaller sizes in smaller genomes
	my $sspace_n = int(log2($insert_size)*1.25+.99); # min overlap required to merge contigs. smaller the insert size, the more certain we can be about a short overlap being correct, so scale this up/down according to insert size
	my $sspace_k = int(log($exp_links)/log(1.4)-11.5);
	# rationale: paired-end: the chimerism rate in small paired end libraries is low. in cases where there is ambiguous pairing, 
	#            the -a parameter will recognize this and resolve it
	#            mate-pair: the chimerism rate in mate-pair libraries can be high, up to 25% or more, but using the soup 
	#            strategy lowers the rate to 1% or less. At this low rate, even single links for scaffolds are almost always reliable

#	$sspace_k = 5;	# gives best results on mediterranei

	my $sspace_a = 0.4;	# be less stringent about ambiguity by default, since these will be fixed by the misassembly detector

	if(defined($rescaffold)&&$rescaffold==1){
		# when rescaffolding we want to pick up the low coverage
		# and small pieces that were cut out from misassembled regions
		$sspace_a = 0.3; # be stringent here -- do not want more misassembly
		$sspace_k -= 2; # if( $insert_size > 1500 );
		$sspace_x = 0; # do not extend contigs -- risks further misassembly
		$libfile .= ".strict";
	}

	# require at least 1 link
	$sspace_k = $sspace_k < 1 ? 1 : $sspace_k;

	my $sspace_cmd = "SSPACE -m $sspace_m -n $sspace_n -k $sspace_k -a $sspace_a -o 1 -x $sspace_x ".
                                        "-l $libfile -s $input_fa -b $outbase -d $WD";
#	if (@_) {
#		print STDERR "[a5] Running SSPACE with unpaired reads\n";
#		my $up = shift;
#		$sspace_cmd .= " -u $up";
#	}
	print STDERR "[a5] $sspace_cmd > $WD/$outbase.out\n";
	`$DIR/SSPACE/$sspace_cmd > $WD/$outbase.out`;
	`rm -rf $WD/bowtieoutput/ $WD/reads/`;
	return "$WD/$outbase.final.scaffolds.fasta";
}

=item validate_filename
Check for characters that are not shell-safe
=cut
sub validate_filename {
	my $fname = shift;
	return $fname !~ /[\ \"\'\&\:\%\?\*\$\<\>]/;
}

sub validate_sequence_file{
	my $seq_file = shift;
	my $message = shift;
	die "ERROR: The file \"$seq_file\" $message.\nPlease ensure that all sequence input file names are free from unsupported characters\n" unless validate_filename($seq_file);
	die "Unable to read input FastQ file $seq_file" unless -e $seq_file && -s $seq_file > 0;
}

=item shell_safe
Escape a file name for the shell.
=cut
sub shell_safe {
	return $_[0] if $_[0] =~ /^[a-zA-Z0-9_\-]+\z/;
	my $s = $_[0];
	$s =~ s/'/'\\''/g;
	return "'$s'";
}


sub log2 {
	my $n = shift;
	return (log($n)/ log(2));
}
=item get_genome_size
Read the given fasta file, calculating the size the genome in said file.
=cut
sub get_genome_size {
	my $fasta = shift;
	open( FASTA, $fasta );
	my $len = 0;
	while( my $line = <FASTA> ){
		next if $line =~ /^>/;
		chomp $line;
		$len += length($line);
	}
	return $len;
}

=item get_insert
Calculate the insert size of the given library using the A5 GetInsertSize.jar code.

Takes 4 parameters: $r1fq, $r2fq, $outbase, $ctgs
where
	$r1fq = read 1 in the paired library
	$r2fq = read 2 in the paired library
	$outbase = the basename for output generated throughout this subroutine
	$ctgs = the contigs to use to calculate the insert size
=cut
sub get_insert($$$$) {
	my $r1fq = shift;
	my $r2fq = shift;
	my $outbase = shift;
	my $ctgs = shift;
	my $estimate_pair_count = 40000;
	my $require_reads = 4000;
	my $fq_linecount = $estimate_pair_count*4;
	# estimate the library insert size with bwa
	# just use a subsample of $estimate_pair_count reads
	`$DIR/bwa index $ctgs`;
	# check that there are enough reads in the file, don't want to duplicate reads...
	my $headcheck = `head -n $fq_linecount $r1fq | wc -l`;
	chomp $headcheck;
	print STDERR "[a5] Found only ".($headcheck/4)." read pairs in library\n" if($headcheck < $fq_linecount);
	$fq_linecount = int($headcheck / 8)*8;
	# allow insert size estimates with fewer than 5000 reads if enough have mapped
	# but never use fewer than 1000.
	$require_reads = ($fq_linecount/4)*0.25 < $require_reads ? ($fq_linecount/4)*0.25 : $require_reads;
	$require_reads = 1000 > $require_reads ? 1000 : $require_reads;
	print STDERR "[a5] Using ".($fq_linecount/4)." read pairs for mapping\n";
	my $half_lines = $fq_linecount / 2;
	`head -n $half_lines $r1fq > $r1fq.sub`;
	`tail -n $half_lines $r1fq >> $r1fq.sub`;
	`head -n $half_lines $r2fq > $r2fq.sub`;
	`tail -n $half_lines $r2fq >> $r2fq.sub`;
	# use bwa to map the reads
	my $cmd = "$DIR/bwa mem -t $threads $ctgs $r1fq.sub $r2fq.sub ".
			"> $outbase.sub.pe.sam 2> $outbase.sampe.out";
	system($cmd);
	$cmd = "java -jar ".JAVA_HEAP." $DIR/GetInsertSize.jar $outbase.sub.pe.sam";
	print STDERR "[a5] $cmd\n"; 
	my $cmdout = `$cmd 2> /dev/null`;
	chomp $cmdout;
	my ($ins_mean,$ins_sd,$ins_n, $ori) = split(/,/,$cmdout);
	my $ins_error = 0;
	`rm $r1fq.sub* $r2fq.sub* $ctgs.* $outbase.sub.pe.sam $outbase.sampe.out`;
	my $n_sd = 6;
	if ($ins_n > $require_reads) {		
		$ins_error = $ins_sd*$n_sd < $ins_mean ? $ins_sd*$n_sd / $ins_mean : 0.95;
		$ins_mean = sprintf("%.0f",$ins_mean);
		$ins_error = sprintf("%.3f",$ins_error);
		$ins_error =~ s/0+$//g;
	} else {
		print STDERR "[a5] Discarding estimate. Not enough data points: $ins_n but requiring $require_reads\n";
		$ins_mean *= -1; 
		$ins_error *= -1;
	}
	return ($ins_mean, $ins_error, $ori);
}

=item finalize
Calcualtes and reports summary statistics on the final scaffolded genome assembly

One parameter: stats, which is a hash containing some previously calculated statistics to be reported
=cut
sub finalize {
	my $stats = shift;
	my $assembly = shift;
	# check whether we need to generate an AGP for NCBI submission
	# any scaffold containing more than 10 N's needs to be reported with an AGP
	# calculate some basic assembly stats while we're at it
	my $agp = 0;
	open(SCAFIN, $assembly);
	my @lens;
	my $total = 0;
	my ($gc,$at) = 0;
	while(my $line = <SCAFIN>){
		next if $line =~ /^>/;
		$agp = 1 if $line =~ /NNNNNNNNNN/;
		push(@lens, length($line)-1);
		$total += length($line)-1;
		$gc += $line =~ tr/GCgc//;
		$at += $line =~ tr/ATat//;
	}
	@lens = sort{$b<=>$a} @lens;
	my ($i, $cumsum) = (0,0);
	for(; $i < @lens; $i++){
		$cumsum += $lens[$i];
		last if $cumsum > $total / 2;
	}
	my $pct_gc = substr(100*$gc/($gc+$at), 0, 4);
	print "[a5] Total size: $total nt, scaffolds: ".@lens.", max $lens[0], N50 $lens[$i], $pct_gc% GC\n";
	my $ec_reads_pct = $total_reads > 0 ? sprintf("%.2f", 100*$EC_reads/$total_reads) : "unknown";
	my $ec_bases_pct = $total_reads > 0 ? sprintf("%.2f", 100*$EC_bases/$total_bases) : "unknown";
	my $x_cov = $total_reads > 0 ? sprintf("%.2f", $EC_bases/$total) : "unknown";
	my $raw_cov = $total_reads > 0 ? sprintf("%.2f", $total_bases/$total): "unknown";
	my $contig_num = `grep -c '>' $OUTBASE.contigs.fasta`;
	chomp($contig_num);
	print "[a5] Total reads used: $total_reads, Error Corrected reads used: $EC_reads, Pct of reads that passed EC step: $ec_reads_pct\n";
	print "[a5] Total bases used: $total_bases, Error Corrected bases used: $EC_bases, Pct of bases that passed EC step: $ec_bases_pct\n";
	print "[a5] Raw coverage: $raw_cov\n";
	print "[a5] EC coverage: $x_cov\n";
	open(STATS, ">$OUTBASE.assembly_stats.csv");
	print STATS "File Name\tContigs\tScaffolds\tGenome Size\tLongest Scaffold\tN50\tRaw reads\tEC Reads\t% reads passing EC\tRaw nt\tEC nt\t% nt passing EC\tRaw cov\tEC cov\tMedian cov\t10th percentile cov\tbases >= Q40\t% GC\tassembler version\n";
	print STATS "$OUTBASE\t$contig_num\t".@lens."\t$total\t$lens[0]\t$lens[$i]\t$total_reads\t$EC_reads\t$ec_reads_pct\t$total_bases\t$EC_bases\t$ec_bases_pct\t$raw_cov\t$x_cov\t$stats->{median_cov}\t$stats->{cov_10}\t$stats->{above_q40}\t$pct_gc\tA5-miseq $VERSION\n";
	close STATS;
	if($agp){
		print "[a5] Large scaffold gaps present in assembly...generating an AGP for NCBI submission\n";
		`$DIR/fasta2agp.pl $OUTBASE.final.scaffolds.fasta > $OUTBASE.agp`;
	}
}


