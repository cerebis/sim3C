meta-sweeper
============

By expressing microbial community composition and the details of metagenomic
experiments in parametric form, meta-sweeper aims to permit the assessment of
analytical methods under variation in both community and experiment.


### prerequisites
This workflow requires a system running linux on an x86-64 machine.

The workflow depends on the following software, which must be installed prior to using it:

* [biopython](http://biopython.org/)

* [pandas](http://pandas.pydata.org/)

* [PyYaml](http://pyyaml.org/)

* [networkx](https://networkx.github.io/)

* [pysam](https://github.com/pysam-developers/pysam)

On ubuntu these can be installed by running:
```bash
sudo apt-get install python-biopython python-pandas python-yaml python-networkx python-pysam
```

On other systems pip may be the preferred method for installation.

### usage

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
