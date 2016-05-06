#!/usr/bin/perl
### david.studholme@tsl.ac.uk
### modified for input from A5 assemblies and updated NCBI requirements for contigs: aaron.darling@uts.edu.au
### Bio::SeqIO requirement removed by Guillaume Jospin: gjospin@ucdavis.edu
### Generates contigs (in FastA) and scaffolding information (in AGP) from Velvet 'contigs.fa' supercontigs file

### Use entirely at you own risk!! There may be bugs!

use strict;
use warnings;

my $sequence_file = shift or die "Usage: $0 <sequence file>\n" ;


### Output file for contigs in Fasta format
my $fasta_outfile = "$sequence_file.contigs.fsa";
open (FILE, ">$fasta_outfile") and
    warn "Will write contigs to file '$fasta_outfile' and AGP to STDOUT\n" or
    die "Failed to write to file '$fasta_outfile'\n";
print "# Generated from assembly file $sequence_file using script $0\n";


my $i = 0;# a counter, used for generating unique contig names

my %sequences = ();
my %desc = ();
open(INSEQ, $sequence_file);
my $seq = "";
my $seqID = "";
my $seqDesc = "";
while (<INSEQ>) {
    chomp($_);
    if($_ =~ m/^>(\S+)\s?(.*)$/){
	if ($seq ne ""){
	    $sequences{$seqID} = $seq;
	    $desc{$seqID}= $seqDesc;
	}
	$seqID = $1;
	$seqDesc = $2;
	$seq = "";
    }
    else{
	$seq .= $_;
    }
}
#adding the last sequence read.
$sequences{$seqID}=$seq;
$desc{$seqID}=$seqDesc;
close(INSEQ);
#print "Found : ".scalar(keys(%sequences))."\n";

foreach my $sID(keys %sequences){


    my $supercontig_id = $sID ;
    my $supercontig_seq = $sequences{$sID} ;
    my $supercontig_desc = $desc{$sID} ;
    my $supercontig_length = length($supercontig_seq);
    
   ### NCBI do not allow coverage and length information in the FastA identifier
   ### e.g. NODE_1160_length_397673_cov_14.469489 is an illegal FastA ID
   ### So we will replace these with simple numbers
   if ($supercontig_id =~ m/NODE_(\d+)_length_\d+_cov_\d+/ or
   $supercontig_id =~ m/^(\d+)$/) {
   $supercontig_id = "scf_$1";
   }
   $supercontig_id =~ s/\|size\d+//g;

   my $j = 1;# another counter, to generate part IDs
   my $start_pos = 1; # keep track of whereabouts in this supercontig we are
    my %substring_sequences;

   ## contigs < 200nt must be eliminated as NCBI does not accept these. Replace with N and merge to neighboring contig
   my $subseqs_join;
   foreach my $substring_sequence ( split /(N{10,})/i, $supercontig_seq ) {
     if (  $substring_sequence =~ m/^N+$/i ) {
        $subseqs_join .= $substring_sequence;
     }elsif(length($substring_sequence) >= 200){
        $subseqs_join .= $substring_sequence;
     }else{
        $subseqs_join .= "N" x length($substring_sequence);
     }
   }
   
   # can't include this if the total amount of resolved sequence is < 200nt
   my $nfree_seq = $subseqs_join;
   $nfree_seq =~ s/N//g;
   next if length($nfree_seq) < 200;
   $subseqs_join =~ s/^N+//g;
   $subseqs_join =~ s/N+$//g;

   foreach my $substring_sequence ( split /(N{10,})/i, $subseqs_join ) {
   ### NB that NCBI do not allow gaps of fewer than 10 nucleotides between contigs.
   ### Gaps of fewer than 10 nucleotides are treated as ambiguities rather than gaps.
   ### So this split is a bit of a fudge.

   #warn "\n$substring_sequence\n" if $supercontig_id eq '1160'; for #debugging only

   ### Define the AGP column contents
   my $object1 = $supercontig_id;
   my $object_beg2 = $start_pos;
   my $object_end3 = $start_pos + length($substring_sequence) - 1;
   my $part_number4 = $j;
   my $component_type5;
   my $component_id6a;
   my $gap_length6b;
   my $component_beg7a;
   my $gap_type7b;
   my $component_end8a;
   my $linkage8b;
   my $orientation9a;
   my $filler9b;
   $j++;
     if (  $substring_sequence =~ m/^N+$/i ) {
       ### This is poly-N gap between contigs
       $component_type5 = 'N';
       $gap_length6b = length($substring_sequence);
       $gap_type7b = 'scaffold';
       $linkage8b = 'yes';
       $filler9b = 'paired-ends';
         } elsif ( $substring_sequence =~ m/^[ACGTRYSWKMBDHVN]+$/i ) {
       ### This is a contig
       $i++; # a counter, used for generating unique contig names
       $component_type5 = 'W';
       $component_id6a = "contig_$i";
       $component_beg7a = 1;
       $component_end8a = length($substring_sequence);
       $orientation9a = '+';
             ### Print FastA formatted contig
       print FILE ">$component_id6a\n$substring_sequence\n";
   } else {
       die "Illegal characters in sequence $sID\n";#$substring_sequence\n";
   }
     $start_pos += length ($substring_sequence);
     if ($component_type5 eq 'N') {
       ### print AGP line for gap
       print "$object1\t$object_beg2\t$object_end3\t$part_number4\t$component_type5\t$gap_length6b\t$gap_type7b\t$linkage8b\t$filler9b\n";
   } else {
       ### print AGP line for contig
       print "$object1\t$object_beg2\t$object_end3\t$part_number4\t$component_type5\t$component_id6a\t$component_beg7a\t$component_end8a\t$orientation9a\n";
         }
   }
} 
