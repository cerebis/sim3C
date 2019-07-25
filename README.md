# sim3C

Read-pair simulation of 3C-based sequencing methodologies (HiC, Meta3C, DNase-HiC)

## Recent changes

The codebase has been refactored into an installable package format. Users will now find it much easier to get sim3C running on their systems. This work is also in preparation for updating sim3C to Python3.

### Changes

- Version incremented to 0.2
- sim3C is now installed as an executable package and run as `sim3C`
- Logging is now used for simultaneous console and file logs.
- The monolithic sim3C.py has been broken into smaller concerns.


## Installation

To install and run sim3C you will require Python 2.7 and LLVM. We recommend that users employ runtime environments such as virtualenv or conda. In particular, conda will make it easy to satisfy the runtime requirement of LLVM. 

The sim3C can be installed directly from Github using Pip as follows.

```bash
pip install git+https://github.com/cerebis/sim3C
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
3. relative abundance
4. chromosome copy number.

There is a mismatch between the hierarchical nature of a community and the flat nature of this simple format. Despite the repetition that can occur for more complicated profiles, we have chosen to stick with this format for simplicity for the time being.

It is easiest to regard column 1 as the primary column, for which each entry must be unique. Column 2 is inherently redundant when dealing with multi-chromosomal cell definitions. The third column refers to the abundance of the cell, and so is as equally redundant as column 2. The forth column allows users to increase the number of copies of a chosen chromosome within a cell. Optional comments are prefixed with a #.

Eg. A simple 2 cell definition, where the first has 2 replicons.
```
#chrom  cell  abundance  copy_number
seq1    bac1  0.4        1
seq2    bac1  0.4        4
seq3    bac2  0.6        1
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
> sim3C --dist uniform -n 500000 -l 150 -e NlaIII -m hic myref.fasta sim.fastq
```

If a community profile has been prepared and we wish to simulate Meta3C.
```bash
> sim3C --profile mycom.txt -n 500000 -l 150 -e HpyCH4IV -m meta3c myref.fasta sim.fastq
```

Both a random seed and a output profile name can be specified at runtime. These make reducibility possible. The random seed is used to initialise all number generators within the simulation and, if given, the profile name will allow Sim3C to save the state of the profile when drawn at random from a distribution. Though saving the profile state is not necessary to reproducibly rerun Sim3C, it assists downstream analyses which may wish to know the true state.

### Useful options

#### Ambiguous IUPAC symbols

```--convert```

At present, Art.py is not able to model errors when reference sequenes contain ambiguous symbols other than N (i.e. MRWSYKVHDB). In these cases, if users do not wish to prepare sequences themselves, the `--convert` option will convert all such symbols to N in memory, prior to simulation. Therefore, emitted simulated reads will contain N in these locations.

#### Faster simulation

```--simple-reads```

Although Sim3C can simulate read-errors, by use of art_illumina[1] machine profiles, there is currently a significant performance hit. If users are interested in faster simulations, possibly to explore a wider space more quickly before a more thorough validation, simple reads without error are possible.

Reference sequences can contain ambiguous symbols (i.e. MRWSYKVHDB) when using the simple read mode.

#### Compress output

```--compress {gzip, bzip2}``` OR ```-C {gzip, bzip2}```

Write the output FASTQ in either gzip or bzip2 compressed format.

#### Specify restriction digest enzyme

```--enzyme [string]``` OR ```-e [string]```

For HiC and Meta3C simulation, an enzyme is required. The default is the 4-cutter NlaIII. The name is case-sensitive and supports most enzymes defined in ReBase[2], as implemented in BioPython Restriction.

## References

1. Huang, Weichun, Leping Li, Jason R. Myers, and Gabor T. Marth. 2012. “ART: A next-Generation Sequencing Read Simulator.” Bioinformatics  28 (4). Oxford University Press: 593–94.

2. Roberts, Richard J., Tamas Vincze, Janos Posfai, and Dana Macelis. 2015. “REBASE--a Database for DNA Restriction and Modification: Enzymes, Genes and Genomes.” Nucleic Acids Research 43 (Database issue): D298–99.
