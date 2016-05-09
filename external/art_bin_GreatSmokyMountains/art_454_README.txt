ART_454 Version 2.6.0 README, last updated on 03/18/2015, Weichun Huang at whduke@gmail.com                              

DESCRIPTION

	ART_454 is a simulation program to generate sequence read data of Roche 454
      	Pyrosequencing sequencers. ART generates reads according to the empirical read quality profile
       	and the calibrated error profile of uncall/overcall homopolymers from real 454 read data. ART 
	has been using for testing or benchmarking a variety of method or tools for next-generation 
	sequencing data analysis, including read alignment, de novo assembly, detection of SNP, CNV, or
       	other structure variation.

	art_454 can generate both single-end and paired-end of 454 sequencing platform. Besides for regular
	genome DNA and cDNA sequencing simulation, art_454 also supports amplicon sequencing. The reference
	sequences can be either DNA or RNA.  ART outputs include FASTQ read, ALN alignment, STAT read coverage,
       	and optional SAM alignment files.  ALN files can also be converted to UCSC BED files by using 
	aln2bed.pl program included. 

	art_454 comes with art_profiler_454, the tool for generating a art_454 read profile from 454 sequencing
       	data in the fastq format. The tool in the folder ART_profiler_454. Please see the README under in the
       	folder for details 
       	
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

	In the "examples" subdirectory, the shell script "run_test_examples_454.sh" gives two examples of using
       	art_454 for read simulation.  To test these two examples, just run the script "run_test_examples_454.sh"

USAGES

	SINGLE-END SIMULATION
		art_454 [-s] [-a ] [-t] [-r rand_seed] [ -p read_profile ] [ -c num_flow_cycles ] <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <FOLD_COVERAGE>
	
	PAIRED-END SIMULATION
		art_454 [-s] [-a ] [-t] [-r rand_seed] [ -p read_profile ] [ -c num_flow_cycles ] <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <FOLD_COVERAGE> <MEAN_FRAG_LEN> <STD_DEV>
	
	AMPLICON SEQUENCING SIMULATION
		art_454 [-s] [-a ] [-t] [-r rand_seed] [ -p read_profile ] [ -c num_flow_cycles ] <-A|-B> <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <#_READS/#_READ_PAIRS_PER_AMPLICON>
	
	===== PARAMETERS =====
	
	INPUT_SEQ_FILE           -  the filename of DNA/RNA reference sequences in FASTA format
	OUTPUT_FILE_PREFIX       -  the prefix or directory of output read data file (*.fq) and read alignment file (*.aln)
	FOLD_COVERAGE            -  the fold of read coverage over the reference sequences
	MEAN_FRAG_LEN            -  the average DNA fragment size for paired-end read simulation
	STD_DEV                  -  the standard deviation of the DNA fragment size for paired-end read simulation
	#READS_PER_AMPLICON      -  number of reads per amplicon (for 5'end amplicon sequencing)
	#READ_PAIRS_PER_AMPLICON -  number of read pairs per amplicon (for two-end amplicon sequencing)
	
	===== OPTIONAL PARAMETERS =====
	
	-A indicate to perform single-end amplicon sequencing simulation
	-B indicate to perform paired-end amplicon sequencing simulation
	-M indicate to use CIGAR 'M' instead of '=/X' for alignment match/mismatch
	-a indicate to output the ALN alignment file
	-s indicate to output the SAM alignment file
	-d print out warning messages for debugging
	-t indicate to simulate reads from the built-in GS FLX Titanium profile [default: GS FLX profile]
	-r specify a fixed random seed for the simulation (to generate two identical datasets from two different runs)
	-c specify the number of flow cycles by the sequencer [ default: 100 for GS-FLX, and 200 for GS-FLX Titanium ] 
	-p specify user's own read profile for simulation
	   NOTE: the name of a read profile is the directory containing read profile data files.
	         please read the REAME file about the format of 454 read profile data files and.
	         and the default filenames of these data files.
	
	===== EXAMPLES =====
	
	 1) singl-end simulation with 20X coverage
		art_454 -s seq_reference.fa ./outdir/single_dat 20
	 2) paired-end simulation with the mean fragment size 500 and STD 20 using GS FLX Titanium platform
		art_454 -s -t seq_reference.fa ./outdir/paired_dat 10 1500 20
	 3) paired-end simulation with a fixed random seed
		art_454 -s -r 777 seq_reference.fa ./outdir/paired_fxSeed 10 2500 50
	 4) single-end amplicon sequencing with 10 reads per amplicon
		art_454 -A -s amplicon_ref.fa ./outdir/amp_single 10
	 5) paired-end amplicon sequencing with 10 read pairs per amplicon
		art_454 -B -s amplicon_ref.fa ./outdir/amp_paired 10
	 
READ PROFILE FILES 

	CREATE A NEW PROFILE
		Please use the tool art_profiler_454 included in the distribution to generate new art_454 read profiles
	       	from new 454 sequencing data in the fastq format. The tool in the folder ART_profiler_454. Please see 
		the README under in the folder for usage. 

	454 READ PROFILE
	    The name of a read profile is the filename of directory containing read profile data files. a read profile
	    should consist of the three required profile files, and one optional profile file.  

	    The three REQUIRED profiles are: qual_1st_profile, qual_mc_profile, length_dist. 

	    	1) qual_1st_profile   - the empirical homopolymer-length dependent quality distribution of the 1st base of DNA homopolymers 

	    	2) qual_mc_profile    - 1st order Markov Chain-based empirical quality distribution of the remaining bases of DNA homopolymers 

	    	3) length_dist        - the empirical distribution of 454 read length 
	    
	    The OPTIONAL profile is indel_erro_profile. When the profile is not given, the built-in indel_erro_profile will be used. 

	  	  indel_erro_profile - undercall or overcall error profiles of DNA homopolymers 

	FILE FORMATS

       		all read profile files are tab-delimited text files. For detailed formats, please see the real examples the under
	       	the folder "454_profiles/GS_FLX_profile" or "454_profiles/GS_FLX_Titanium_profile" for file formats of these four files. 
	
CREATION OF YOUR OWN 454 READ PROFILE 

	Users can use the program "454_readprofile_art" that comes with the software package to generate 454 read profile from
       	454 sequencing data in the fastq format.

	USAGE:
       		./454_readprofile_art input_fastq_files_dir output_profile_dir [fastq_filename_extension]
      		
	PARAMETERS 

		1) input_fastq_files_dir - directory containing all 454 read data files in the fastq format 

		2) output_profile_dir - directory to put read profile files generated by the program

	EXAMPLES:
       		perl 454_readprofile_art 454_dat_dir 454_profile_dir

	       	perl 454_readprofile_art 454_dat_dir 454_profile_dir fastq
		       
	NOTES:
       		1)the default filename extension for fastq files is fq

	       	2)the program can read gzipped fastq data files with filename extension fq.gz


OUTPUT FILES

	*.fq   - FASTQ read data files. For paired-end read simulation, *1.fq contains data of the first reads, and *2.fq for the second reads.
	*.aln  - read alignment files. For paired-end read simulation, *1.aln has read alignments for the first reads and *2.aln for the second reads.
	*.sam  - SAM alignment file. 
	*.stat - read coverage files with the following three tab-delimited fields per lines:
		 (1)reference position  
		 (2)number of reads starting at the position  
		 (3)number of reads covering the position  

	FASTQ FILE format 
		A FASTQ file contains both sequence bases and quality scores of sequencing reads and is in the following format:  
			@read_id 
			sequence_read 
			+ 
			base_quality_scores 
	
		A base quality score is coded by the ASCII code of a single character, where the quality score is equal to ASCII code of the
       		character minus 33.    

		Example: 
		@refid-4028550-1 
		caacgccactcagcaatgatcggtttattcacgat...
		+ 
		????????????7?????<??>??=&?<<?-<?0?...

	ALN format 
		An ALN file has a Header and main Body parts. The header part includes the command used to generate this file and reference
	       	sequence id and length. The header @CM tag for command line, and @SQ for reference sequence.  A header always starts with 
		"##ART" and ends with  "##Header End".

		HEADER EXAMPLE

		##ART_454
       		@CM     ../../bin/MacOS64/art_454 -s ./testSeq.fa ./single_454_test 10 
		@SQ     seq1    7207
	       	@SQ     seq2    3056
		##Header End

		The body part contains all alignments in the following format 

		>ref_seq_id	read_id	aln_start_pos	ref_seq_strand
	       	ref_seq_aligned
	       	read_seq_aligned 
	
		aln_start_pos is the alignment start position of reference sequence. aln_start_pos is always relative to the strand of reference
       		sequence. That is, aln_start_pos 10 in the plus (+) strand is different from aln_start_pos 10 in the minus (‐) stand.  
	
		ref_seq_aligned is the aligned region of reference sequence, which can be from plus strand or minus strand of the reference sequence.  
		read_seq_aligned is the aligned sequence read, which always in the same orientation of the same read in the corresponding fastq file. 

	SAM format

       		SAM is a standard format for next-gen sequencing read alignments. The details of the format and examples are available at the links below:	
		1) http://samtools.sourceforge.net/SAM1.pdf
	       	2) http://genome.sph.umich.edu/wiki/SAM		


	BED format

		See the format at UCSC http://genome.ucsc.edu/FAQ/FAQformat.html#format1

		NOTE: both ALN and BED format files use 0-based coordinate system while SAM format uses 1-based coordinate system.


