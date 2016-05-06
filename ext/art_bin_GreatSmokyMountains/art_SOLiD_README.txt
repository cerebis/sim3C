ART_SOLiD (version 1.3.3) README (last updated on 03/19/2015) Weichun Huang at whduke@gmail.com                              

DESCRIPTION:

	ART_SOLiD  is a simulation program to generate sequence read data of SOLiD sequencing reads.
       	ART generates reads according to a SOLiD read error profile. The built-in error profile is an
       	empiricial error profile summarized from large SOLiD sequencing data. ART has been using for 
	testing or benchmarking a variety of method or tools for next-generation sequencing data analysis,
       	including read alignment, de novo assembly, detection of SNP, CNV, or other structure variation.
       	
	art_SOLiD can generate both single-end, matepair, and paired-end of SOLiD sequencing platform.
       	art_SOLiD also support amplicon sequencing simulation with RNA references. Its outputs include 
	FASTQ read, MAP alignment, and optional SAM alignment files. The map2bed.pl tool included can
       	convert ART MAP alignment files to a UCSC BED file 	

COMPILATION AND INSTALLATION:

	PREREQUISITES: 
			1) GNU g++ 4.0 or above (http://gcc.gnu.org/install) 
			2) GNU gsl library (http://www.gnu.org/s/gsl/) 

	COMPILATION & INSTALLATION

		1) add gsl library installation directory to compiler's search path

	       	If your GSL library is not installed in system standard search path, your GSL installation 
		directory need to be added to gcc/g++ complier FLAGS. For example, gsl is typically installed
	       	in /opt/local directory in MacOS X. You can add the directory to the search path of gcc/g++ 
		complier by run the following command: 
		
		export CFLAGS="$CFLAGS -I/opt/local/include" CPPFLAGS="$CPPFLAGS -I/opt/local/include" LDFLAGS="$LDFLAGS -L/opt/local/lib"

		2) compile and install 

		./configure --prefix=$HOME
	       	make
	       	make install

EXAMPLES

	In the "examples" subdirectory, the shell script "run_test_examples_SOLiD.sh" gives examples of using
       	art_SOLiD for read simulation.  To test these two examples, just run the script "run_test_examples_SOLiD.sh"

USAGES

	SINGLE-END (F3 READ) SIMULATION
		art_SOLiD [ options ] <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <LEN_READ> <FOLD_COVERAGE>
	
	MATE-PAIR READS (F3-R3 PAIR) SIMULATION
		art_SOLiD [ options ] <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <LEN_READ> <FOLD_COVERAGE> <MEAN_FRAG_LEN> <STD_DEV>
	
	PAIRED-END READS (F3-F5 PAIR) SIMULATION
		art_SOLiD [ options ] <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <LEN_READ_F3> <LEN_READ_F5> <FOLD_COVERAGE> <MEAN_FRAG_LEN> <STD_DEV>
	
	AMPLICON SEQUENCING SIMULATION
		art_SOLiD [ options ] -A s <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <LEN_READ> <READS_PER_AMPLICON>
		art_SOLiD [ options ] -A m <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <LEN_READ> <READ_PAIRS_PER_AMPLICON>
		art_SOLiD [ options ] -A p <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <LEN_READ_F3> <LEN_READ_F5> <READ_PAIRS_PER_AMPLICON>
	
	===== PARAMETERS =====
	
	INPUT_SEQ_FILE            -  filename of DNA/RNA reference sequences in FASTA format
	OUTPUT_FILE_PREFIX        -  prefix or directory for all output data files
	FOLD_COVERAGE             -  fold of read coverage over the reference sequences 
	LEN_READ                  -  length of F3/R3 reads
	LEN_READ_F3               -  length of F3 reads for paired-end read simulation
	LEN_READ_F5               -  length of F5 reads for paired-end read simulation
	READS_PER_AMPLICON        -  number of reads per amplicon
	READ_PAIRS_PER_AMPLICON   -  number of read pairs per amplicon
	MEAN_FRAG_LEN             -  mean DNA/RNA fragment size for matepair/paired-end read simulation
	STD_DEV                   -  standard deviation of the DNA/RNA fragment sizes for matepair/paired-end read simulation
	
	===== OPTIONAL PARAMETERS =====
	
	-A specify the read type for amplicon sequencing simulation (s:single-end, m: matepair, p: paired-end)
	-M indicate to use CIGAR 'M' instead of '=/X' for alignment match/mismatch
	-s indicate to generate a SAM alignment file
	-r specify the random seed for the simulation
	-f specify the scale factor adjusting error rate (e.g., -f 0 for zero-error rate simulation)
	-p specify user's own read profile for simulation
	
	===== EXAMPLES =====
       	1) singl-end 25bp reads simulation at 10X coverage
		art_SOLiD -s seq_reference.fa ./outdir/single_dat 25 10
        2) singl-end 75bp reads simulation at 20X coverage with users' error profile
       		art_SOLiD -s -p ../SOLiD_profiles/pseudo_profile ./seq_reference.fa ./dat_userProfile 75 20
	3) matepair 35bp (F3-R3) reads simulation at 20X coverage with DNA MEAN fragment size 2000bp and STD 50
		art_SOLiD -s seq_reference.fa ./outdir/matepair_dat 35 20 2000 50
        4) matepair reads simulation with a fixed random seed
		art_SOLiD -r 777 -s seq_reference.fa ./outdir/matepair_fs 50 10 1500 50
	5) paired-end reads (75bp F3, 35bp F5) simulation with the MEAN fragment size 250 and STD 10 at 20X coverage
		art_SOLiD -s seq_reference.fa ./outdir/paired_dat 75 35 50 250 10
        6) amplicon sequencing with 25bp single-end reads at 100 reads per amplicon
		art_SOLiD -A s -s amp_reference.fa ./outdir/amp_single 25 100
	7) amplicon sequencing with 50bp matepair reads at 80 read pairs per amplicon
		art_SOLiD -A m -s amp_reference.fa ./outdir/amp_matepair 50 80
        8) amplicon sequencing with paired-end reads (35bp F3, 25bp F5 reads) at 50 pairs per amplicon
		art_SOLiD -A p -s amp_reference.fa ./outdir/amp_pair 35 25 50
	 
OUTPUT DATA FILES

	*.fq   - FASTQ read data files. For matepair/paired-end read simulation, *_R3.fq/*_F5.fq contains data of the first reads, and *_F3.fq for the second reads.
	*.map  - read alignment files. For matepair read simulation, *_R3/*_F5.map has read alignments for the first reads and *_F3.map for the second reads.
	*.sam  - SAM read alignment files 

	FASTQ file format 

		A FASTQ file contains both sequence bases and quality scores of sequencing reads and is in the following format:  

			@read_id 
			sequence_read 
			+ 
			base_quality_scores 
	
		A base quality score is coded by the ASCII code of a single character, where the quality score is equal to ASCII code of the
       		character minus 33.    

		Example: 
		@1_1_1_F3
	       	T10121011322000311302213102311132
	       	+
	       	163=+48./48<347//,=/84-4)77(''-)
	       	@1_1_2_F3
	       	T01213000012110232021321303011332
	       	+
	       	1>7<-653?01:55./.%3+'61-+40(*)#'

	MAP file format 
		A MAP file has a Header and main Body parts. The header part includes the command used to generate this file and reference
	       	sequence id and length. The header @CM tag for command line, and @SQ for reference sequence.  A header always starts with 
		"##ART" and ends with  "##Header End".  ART MAP alignment files can be converted to a UCSC BED file by the map2bed.pl tool.

		HEADER EXAMPLE

		##ART_SOLiD     read_length     32
		@CM     ../../bin/MacOS64/art_SOLiD -s ./testSeq.fa ./single_SOLiD_test 32 10 
		@SQ     seq1    7207
	       	@SQ     seq2    3056
		##Header End

		The body part of a MAP file is a tab-delimited text alignment data where each line is the mapping information record of a read. The format of each line
		is below:
		
		ref_seq_id	read_id		ref_map_pos	ref_seq_strand	num_sequencing_errors	err_pos_1  WrongCorrectColor_1	err_pos_2  WrongCorrectColor_2	...	
	
		ref_map_pos - the alignment start position of reference sequence. ref_map_pos is always relative to the strand of reference
       		sequence. That is, ref_map_pos 10 in the plus (+) strand is different from ref_map_pos 10 in the minus (‐) stand.  

		num_sequencing_errors - number of sequencing call errors

		err_pos_n - the zero-indexed read position of the n_th call error
		WrongCorrectColor_n  - the n_th call error representing by two digits number with the 1st number the wrong called color and the 2nd the correct color

		Example: 
		
		seq1	1_1_1_F3	6028	+	4	18	31	22	03	25	32	27	23

	SAM format

       		SAM is a standard format for next-gen sequencing read alignments. The details of the format and examples are available at the links below:	
		1) http://samtools.sourceforge.net/SAM1.pdf
	       	2) http://genome.sph.umich.edu/wiki/SAM	

		Read sequences in a SAM file are in the regular DNA base space. As one change in SOLiD reads in the color space can result multiple changes in base space,
	        so not all mismatches in base space between reads and its reference sequence are of SOLiD sequencing errors.

	BED format

		See the format at UCSC http://genome.ucsc.edu/FAQ/FAQformat.html#format1

		NOTE: both MAP and BED format files use 0-based coordinate system while SAM format uses 1-based coordinate system.

ERROR PROFILE FILE FORMAT

	An ART SOLiD read error profile has the following five columns:
		1) read position
	       	2) correct base
	       	3) wrong base called 
		4) probability of calling the wrong base in the 1st read
	       	5) probability of calling the wrong base in the 2nd read

	and an error profile end with "*end_of_profile".

	Please see "profile_default" in the folder SOLiD_profiles for a real example 
