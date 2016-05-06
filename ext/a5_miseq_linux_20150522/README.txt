## Quickstart
    
    a5_pipeline.pl read_1.fastq.gz read_2.fastq.gz mygenome

This will assemble the paired reads in the two files `read_1.fastq.gz` and `read_2.fastq.gz`. These example files are gzip compressed just like standard MiSeq reporter output, but need not be. The final scaffolded assembly will be saved in `mygenome.final.scaffolds.fasta`. 

## Copyright and license

The A5-miseq pipeline is (c) 2011-2014 Andrew Tritt and Aaron Darling. A5-miseq is free, open source software licensed under the GPLv3, please see the file LICENSE for details. The A5-miseq pipeline includes open-source components developed and copyright by 3rd parties: `bwa`, `samtools`, `SGA`, `bowtie`, `IDBA-UD`, `SSPACE`, and `Trimmomatic`. Source code for these components is available from the following locations: 

`bwa`: https://sourceforge.net/projects/bio-bwa/
`samtools`: https://sourceforge.net/projects/samtools/
`SGA`: https://github.com/jts/sga
`bowtie`: https://sourceforge.net/projects/bowtie-bio/
`Trimmomatic`: http://www.usadellab.org/cms/?page=trimmomatic

Please see their license agreements for further details. 
The following two components have been modified from their original versions and the corresponding GPL licensed source code is available in the A5-miseq repository:
`IDBA-UD`: https://sourceforge.net/p/ngopt/code/HEAD/tree/trunk/idba-1.1.1/
`SSPACE`: https://sourceforge.net/p/ngopt/code/HEAD/tree/trunk/SSPACE/

## What is A5-miseq?

_A5-miseq_ is a pipeline for assembling DNA sequence data generated on the Illumina sequencing platform. This README will take you through the steps necessary for running _A5-miseq_. 

## What A5-miseq can't do

There are many situations where A5-miseq is not the right tool for the job. In order to produce accurate results, A5-miseq requires Illumina data with certain characteristics. A5-miseq will likely not work well with Illumina reads shorter than around 80nt, or reads where the base qualities are low in all or most reads before 60nt. A5-miseq assumes it is assembling homozygous haploid genomes. Use a different assembler for metagenomes and heterozygous diploid or polyploid organisms. Use a different assembler if a tool like FastQC reports your data quality is dubious. You have been warned!

Datasets consisting solely of unpaired reads are not currently supported.

## Requirements

A5-miseq requires 64-bit Linux (kernel 2.6.15 or later) or Mac OS X 10.6 or later. A Java Runtime Environment is also required. Mac OS X includes Java. On Linux, check with your distribution provider for details on installing Java. 

## Installation

Once you have downloaded and extracted the pipeline, the `a5_pipeline.pl` script can be run in place. Optionally, the pipeline's bin/ directory can be added to the `$PATH` environment variable so that specifying the full path to `a5_pipeline.pl` is not necessary: 
    
    export PATH=$PATH:/path/to/ngopt/bin

Please change /path/to/ngopt/bin appropriately, and put this command in ~/.bashrc or ~/.profile or another script that gets run at login time. The pipeline does not need to be copied to a specific location in order to be installed. You do not need root or superuser or administrator access to install and use the pipeline with this approach. 

## Usage details
    
    Usage: a5_pipeline.pl [--begin=1-5] [--end=1-5] [--preprocessed] <lib_file> <out_base>
    
    Or:    a5_pipeline.pl <Read 1 FastQ> <Read 2 FastQ> <out_base>
    
    Or:    a5_pipeline.pl <Read 1,2 Interleaved FastQ> <out_base>
    
    <out_base> is the base file name for all output files. When assembling from 
    a single library, the fastq files may be given directly on the command line.
    If using more than one library, a library file must be given as <lib_file>.
    The library file must contain the filenames of all read files.
    
    If --preprocessed is used, <lib_file> is expected to be the library file
    created during step 2 of the pipeline, named <out_base>.preproc.libs. Note 
    that this flag only applies if beginning pipeline after step 2.

#### Example usage 1
    
    a5_pipeline.pl my_libs assembly.out

#### Example usage 2
    
    a5_pipeline.pl my_reads.1.fastq my_reads.2.fastq assembly.out

All output will be contained within the directory &lt;output base&gt;. The above example will produce the following directory and files: 

#### Output files

Running the pipeline using either of the two examples will produce the following output: 
    
    assembly.out.ec.fastq.gz              // Error corrected reads
    assembly.out.contigs.fasta            // Contigs
    assembly.out.crude.scaffolds.fasta    // Crude scaffolds:  Not checked for misassemblies
    assembly.out.broken.scaffolds.fasta   // Broken scaffolds: Checked for misassemblies, but not rescaffolded
    assembly.out.final.scaffolds.fasta    // Final scaffolds:  Checked for misassemblies, and rescaffolded
    assembly.out.final.scaffolds.fastq    // Final scaffolds in FastQ format with base call qualities
    assembly.out.final.scaffolds.qvl      // Quality values for the final scaffolds in QVL format for submission to NCBI
    assembly.out.final.scaffolds.agp      // AGP file describing scaffolding for submission to NCBI 
    assembly.out.assembly_stats.csv       // A tab-separated file with assembly summary statistics


### Compute requirements

A5-miseq can assemble microbial genomes on modern laptop hardware, but will require substantial disk space for temporary files in the process. Most bacterial genomes can be assembled with less than 4GB RAM in a few hours. Several steps in the A5-miseq pipeline are parallelized and the extent of parallelism can be controlled with the `--threads` option to A5-miseq and the OMP_NUM_THREADS environment variable.

### If you are using a library file, please read the following

#### Making a library file

The library file must contain file paths to all libraries being assembled. Separate libraries are delimited in the library file by `[LIB]`. For a paired-end or mate-pair library, reads can be passed as two files, one for each end, or as a single shuffled (aka interleaved) file. Additionally, unpaired reads can be passed as a single file. For a single-end library, reads must be passed as a single file. Reads must be in FASTQ format. In addition to file locations, you can also give an insert size for paired reads. Please note that this is NOT necessary. The given number will only be used when _A5-miseq_ is not able to reliably calculate the insert size on its own. 

Example library file: 
    
    [LIB]
    p1=library1_p1.fastq
    p2=library1_p2.fastq
    up=library1_up.fastq
    ins=350
    [LIB]
    shuf=library2.fastq
    up=library2.unpaired.fq
    ins=6500
    [LIB]
    up=unpaired_library.fastq

### Adapter sequences

A5-miseq screens input sequence reads for all standard Illumina genomic DNA library adapter sequences available as of April 2014. This includes screening for Illumina PE v1.0 adapters, TruSeq adapters, and Nextera (XT) adapters. If your libraries were prepared using alternative adapter sequences and these have not been trimmed from the data prior to assembly, it will be necessary to copy the adapter.fasta file provided with A5-miseq and edit it to include the appropriate adapter set. The path to the new adapter file can be provided with the `--adapter=` command line option when running `a5_pipeline.pl`.

### The assembly stats file

A5-miseq produces summary statistics on the assembly in a tab delimited file with the suffix `.assembly_stats.csv`. This file can be opened in a spreadsheet program such as Excel, Numbers, or LibreOffice Calc. If you have assembled many genomes in a particular directory or directories, the summary statistics from all of them can be combined into a single table with the following command:
`cat *.assembly_stats.csv | sort | uniq > all_assembly_stats.csv`
The resulting table will be in a file called `all_assembly_stats.csv`.

The columns of the assembly stats file are as follows:

    Contigs             // Number of contigs generated
    Scaffolds           // Number of scaffolds generated
    Assembly size       // Sum of lengths of all scaffolds in assembly
    Longest Scaffold    // Length of the longest scaffold
    N50                 // N50 value for the scaffolds (N50 length is defined as the length for which 
                        //    the collection of all scaffolds of that length or longer contains at 
                        //    least half of the total of the lengths of the scaffolds) 
    Raw reads           // Number of raw reads used as input
    EC reads            // Number of reads used for assembly, after quality control and error correction
    % reads passing EC  // Percentage of the raw reads that the EC reads represent  
    Raw nt              // Number of nucleotides used as input
    EC nt               // Number of nucleotides after quality control and error correction
    % nt passing EC     // Percentage of the raw nucleotides that the EC nucleotides represent
    Raw cov             // Average depth of sequence coverage provided by the raw data 
    Median cov          // Median actual depth of coverage in the assembly 
    10th percentile cov // 10th percentile depth of coverage -- 90% of sites have greater coverage
    bases >= Q40        // Number of bases that have a PHRED-scale quality greater or equal to Q40
    Assembler version   // Version used to run the assembly.


### Other notes

Various stages of the A5-miseq pipeline estimate assembly parameters from subsets of the input reads. This parameter estimation can vary from run to run, and can be dependent on the order of the input reads. Therefore it is possible that two runs of A5-miseq on the same data will not generate the exact same assembly. 

