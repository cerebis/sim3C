meta-sweeper
============

By expressing microbial community composition and the details of metagenomic
experiments in parametric form, meta-sweeper aims to permit the assessment of
analytical methods under variation in both community and experiment.


=== prerequisites

The following software must be installed prior to using the workflow:

[biopython](http://biopython.org/)
[pandas](http://pandas.pydata.org/)

=== usage

Standard (local) execution
```bash
nextflow -C test.config run sweep.nf
```

SGE execution
```bash
nextflow -C test.config run sweep.nf --profile sge
```

PBSPro execution
```bash
nextflow -C test.config run sweep.nf --profile pbspro
```

[Darling Lab](http://darlinglab.org/)

[i3 institute, UTS
Sydney](http://www.uts.edu.au/research-and-teaching/our-research/ithree-institute)
