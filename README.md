# Sim3C

Read-pair simulation of 3C-based sequencing methodologies (HiC, Meta3C, DNase-HiC)

## Dependencies

*Python 2.7*

*Python modules*

- biopython
- intervaltree
- numpy
- scipy
- tqdm
- PyYAML

Dependencies can be satisfied through pip, using the supplied requirements file.

```bash
pip install -U -r requirements.txt
````

## Intro

Analogous to the well established process of simulating whole-genome shotgun reads, sim3C.py simulates read-pairs as if generated from a sequencing library prepared using a HiC/3C methodology.

Experimental noise is considered both in terms of Illumina-based sequencing error and that inherent in the production of HiC/3C library generation. Sequencing error is afforded by a reimplementation of art_illumina (Huang et al, 2011) in Python (Art.py). HiC/3C noise, in the form of spurious ligation products, is modeled as the uniformly random association of any two sites across the entire extent of the source genome(s), the rate of which is user controlled.

To support community sampling (metagenomes), an abundance profile is supplied along with the set of reference sequences at runtime. The profile can either take the form of an explicit table or be drawn at random from a user chosen probability distribution (equal, uniform random or log-normal).

The tool conceptualises the process of HiC/3C read-pair generation, as a series of conditional probabilities. For intra-chromosomal read-pairs, genomic separation is constrained to follow an empirically determined long-tailed composition of the geometric and uniform distributions. For inter-chromosomal (but same cell) pairs, genomic position is constrained only by proximity to a restriction site. Lastly, inter-cellular read-pairs are not considered.

Fine scale structurally related features that have been observed via contact maps in real experiments (Tung et al, 2013) are also reproduced in our simulation. Namely, contacts between the two arms of the chromosome and chromosomal interacting domains (CID). Within bacterial HiC contact maps, inter-arm contacts are responsible for the fainter anti-diagonal observed (y=-x rather than y=x), while the tightly folded CID domains act to modulate contact frequencies over their local extent, resulting in blocks of density. In our simulation, CIDs are randomly generated at runtime using the supplied seed, and the strength of their effect adjusted or disabled by the user.

Pseudocode
```
If a spurious event:
  unconstrained pairing of any two sites across entire (meta)genome
Else if an inter-chromosomal event: 
  unconstrained positions on Chr_A, Chr_B from Genome_N
Else is an intra-chromosomal event:
  contrained positions x1, x2 ~ g(x) on Chr_A
Where g(x) reflects both the empirical constraint on distance and fine-scale features mentioned above.
```
