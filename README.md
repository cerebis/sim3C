# sim3C

[![CodeFactor](https://www.codefactor.io/repository/github/cerebis/sim3c/badge)](https://www.codefactor.io/repository/github/cerebis/sim3c)

Read-pair simulation of 3C-based sequencing methodologies (HiC, Meta3C, DNase-HiC)

## Recent Updates

- Python 3 support (requires 3.11)
- Minimal Docker image (`cerebis/sim3c`)
- New optional TOML-format community profile definition
  - finer granularity
  - eliminates parameter redundancy
- Approximately 6x faster read generation 
  - testing system: MacOS Intel i9
  - 50% efficiency, 150bp reads, B.subtilis chrom + plasmid, uniform abundance 
- Read-pair output now uses dnaio (https://github.com/marcelm/dnaio)
  - interleaved or split R1/R2 files
  - fasta or fastq format
  - supports gzip and bzip2 compression
- Cython implementation of performance bottlenecks
- Primary random number generation uses Permuted Congruential Generators (PCG) (https://github.com/imneme/pcg-c)
  - Implementation of a cython wrapper for PCG C-library

## Installation

### Using Docker image

A docker image of sim3C has been added `cerebis/sim3C:latest`.

**Example use**
```
docker run --rm -v $PWD:/data cerebis/sim3c --seed 1234 -e Sau3AI -l 150 -n 10000 --insert-mean 300 --insert-sd 50 --profile /data/profile.tsv /data/ref_genomes.fna /data/output_R1.fq.gz /data/output_R2.fq.gz
```

### Local installs

To install and run sim3C you will require Python >=3.11, C-compiler, Make, and LLVM. We recommend that users employ runtime environments such as virtualenv or conda. In particular, conda will make it easy to satisfy the runtime requirement of LLVM.

The sim3C executable can be installed for an individual user directly from Github using Pip as follows.

```bash
pip install --user git+https://github.com/cerebis/sim3C
```

Python dependencies will automatically be satisfied during installation. 

If you encounter problems please visit and log an issue at the [project site on Github](https://github.com/cerebis/sim3C/issues). 

## Usage

### External files

#### Reference Sequence(s) (mandatory)

At a minimum, Sim3C requires a reference sequence (or sequences) from which to draw reads. This reference must be in FASTA format. For multiple references, all must be contained in the single multi-FASTA file. All sequence identifiers must be unique must be unique in a multi-FASTA file.

#### Community Profile (optional)

A community profile can be supplied, which gives the user more control over the definition. Without this enternal profile file, each individual sequence encountered in the supplied reference will be treated as a separate monochromosomal genome.

A profile is a simple tabular text file with the columns:

1. chromosome
2. cell
3. molecule
4. relative abundance
5. chromosome copy number.

There is a mismatch between the hierarchical nature of a community and the flat nature of this simple format. Despite the repetition that can occur for more complicated profiles, we have chosen to stick with this format for simplicity for the time being.

It is easiest to regard the first column as the primary column, for which each entry must be unique. The second column is inherently redundant when dealing with multi-chromosomal cell definitions. The third column groups sequences (e.g. draft genomes) as a single molecule, which permits simulating interactions between grouped sequences as intra-molecular. The fourth column refers to the abundance of the cell, and so is as equally redundant as column 2. The fifth column allows users to increase the number of copies of a chosen chromosome within a cell. Optional comments are prefixed with a #.

**An example definition with four cells**
- cell: e.coli contains two molecules: chromosome and plasmid.
  -  the molecule "chromosome" is in two pieces. This is an example of the new structure.
- molecule names are up to the user but must be the same for all related sequences.
- relative abundance values need not be normalized to sum to 1.
```
#chrom    cell     molecule      abundance    copy_number
contig1   e.coli   chromosome    0.6           1
contig2   e.coli   chromosome    0.2           1
contig3   e.coli   plasmid       0.1           4
contig4   b.subt   chrom_xyz     0.05          1
contig5   s.aur    foobar        0.05          1
```

##### Column definitions

**1. chromosome:** (string)
 
Each chromosome name must match a sequence id within the reference FASTA file, therefore all constraints on FASTA id fields are imposed on this column. For long IDs, such as those in Refseq, no attempt to parse the namespaces is attempted in Sim3C and so the entire ID must be included. 
  
I.E. For the following reference FASTA fragment:
```
>db|foo|bar|12345
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC
TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA
TATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACC
ATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAG
CCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAA
GTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCC
```

The profile line might be

```
db|foo|bar|1234  mycell  1  1
```

**2. cell:** (string)

Cells act as containers of chromosomes. Users can choose any label they desire, baring whitespace. For multi-chromosome cell/genome definitions, this label will be repeated, as it indicates which container in to which the chromosome is placed.

**3. relative abundance:** (float)

Relative abundances are defined per-cell, therefore this value will be repeated for each chromosome belonging to the cell. The abundances do not need to sum to 1 as the profile is normalised internally.

**4. copy number:** (int)

Copy number is most often set to 1, but gives the user the freedom to increase the abundance of chromosomes independent of the cellular abundance.

### Running sim3C

The simplest runtime scenario would be a strictly mono-chromosomal community, which requires only reference FASTA.

Simulate 500k 150bp read-pairs using traditional HiC, NlaIII as an enzyme and uniformly random abundance across all sequences.
```bash
> sim3C --dist uniform -n 500000 -l 150 -e NlaIII -m hic myref.fasta sim.fq
```

If a community profile has been prepared and we wish to simulate Meta3C.
```bash
> sim3C --profile mycom.txt -n 500000 -l 150 -e HpyCH4IV -m meta3c myref.fasta sim.fq
```

Both a random seed and a output profile name can be specified at runtime. These make reducibility possible. The random seed is used to initialise all number generators within the simulation and, if given, the profile name will allow Sim3C to save the state of the profile when drawn at random from a distribution. Though saving the profile state is not necessary to reproducibly rerun Sim3C, it assists downstream analyses which may wish to know the true state.

### Useful options

#### Ambiguous IUPAC symbols

```--convert```

At present, Art.py is not able to model errors when reference sequenes contain ambiguous symbols other than N (i.e. MRWSYKVHDB). In these cases, if users do not wish to prepare sequences themselves, the `--convert` option will convert all such symbols to N in memory, prior to simulation. Therefore, emitted simulated reads will contain N in these locations.

#### Faster simulation

```--simple-reads```

Users whose work does not require simulated read errors -- or for whom time is very short -- sim3C can be run in a "simple-read" mode. In testing, disabling error modelling results in a 60% increase in simulation speed.

**Please Note:** when error modelling is disasbled, if reference sequences contain ambiguous symbols (i.e. MRWSYKVHDB), then these will be carried through to the simulated reads.

#### Output format

Output reads can be written in either FASTA or FASTQ format, where the format is inferred from the file extension specified at runtime. Eg. `.fq|.fastq` -> FASTQ, `.fa|.fasta` -> FASTA.

#### Compress output

Output reads can be compressed using gzip or bzip2, where the compression type is inferred from the file extension specified at runtime. 
Eg. `.gz` -> gzip compression, `.bz2` -> bzip2 compression.

```--compress```

#### Split or Interleaved output

Output reads can be written as interleaved or split R1/R2 files. At runtime, specifying a single output read file will produce interleaved read-pairs, while specifying two output files will produce split R1/R2 files. 

**Please note:** Only the suffixes of file names are inspected, there is no requirement to adhere to a `_1`/`_2` or `_R1`/`_R2` naming convention with split read output.


**Interleaved**
```
sim3C --dist uniform -n 500000 -l 150 -e NlaIII -m hic myref.fasta sim.fq
```

**Split R1/R2**
```
# conventional syntax
sim3C --dist uniform -n 500000 -l 150 -e NlaIII -m hic myref.fasta sim_R1.fq sim_R2.fq
# odd names
sim3C --dist uniform -n 500000 -l 150 -e NlaIII -m hic myref.fasta foo.fq bar.fq
```

#### Examples
```
# uncompressed, interleaved FASTA output
sim3C --dist uniform -n 500000 -l 150 -e NlaIII -m hic myref.fasta sim.fa

# gzip compressed, interleaved FASTQ output
sim3C --dist uniform -n 500000 -l 150 -e NlaIII -m hic myref.fasta sim.fq.gz

# gzip compressed, split R1/R1 FASTQ output
sim3C --dist uniform -n 500000 -l 150 -e NlaIII -m hic myref.fasta sim_1.fq.gz sim_2.fq.gz
```

#### Specify restriction digest enzyme

```--enzyme [string]``` OR ```-e [string]```

For HiC and Meta3C simulation, an enzyme is required. The default is the 4-cutter NlaIII. The name is case-sensitive and supports most enzymes defined in ReBase[2], as implemented in BioPython Restriction.

## References

1. Huang, Weichun, Leping Li, Jason R. Myers, and Gabor T. Marth. 2012. “ART: A next-Generation Sequencing Read Simulator.” Bioinformatics  28 (4). Oxford University Press: 593–94.

2. Roberts, Richard J., Tamas Vincze, Janos Posfai, and Dana Macelis. 2015. “REBASE--a Database for DNA Restriction and Modification: Enzymes, Genes and Genomes.” Nucleic Acids Research 43 (Database issue): D298–99.
