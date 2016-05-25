meta-sweeper
============

By expressing microbial community composition and the details of metagenomic experiments in parametric form, meta-sweeper aims to permit the assessment of analytical methods under variation in both community and experiment.


### Prerequisites
This workflow requires a system running linux on an x86-64 machine due to dependent binaries, while Nextflow itself requires Java 7 or higher and should be installed and run using the same variant.

### Dependencies

The [Nextflow](http://www.nextflow.io/) framework, which can be installed easily using either of the following:

```wget -qO- get.nextflow.io | bash``` or ```curl -fsSL get.nextflow.io | bash```

In addition, meta-sweeper expects that the main executable ``nextflow``` is accessible on the path. Users can move this file to a location already on the path or add its parent directory to the path.

The workflow depends on the following Python modules, which must be installed prior to using it:

* [biopython](http://biopython.org/)

* [pandas](http://pandas.pydata.org/)

* [PyYaml](http://pyyaml.org/)

* [networkx](https://networkx.github.io/)

* [pysam](https://github.com/pysam-developers/pysam)

On ubuntu these can be installed by running:
```bash
sudo apt-get install python-biopython python-pandas python-yaml python-networkx python-pysam
```

On other systems [pip](https://pip.pypa.io/en/stable/) may be the preferred method for installation.
```bash
pip install --upgrade biopython pandas PyYAML networkx pysam
```

### usage

We recommend that users start sweeps by using the ```meta-sweeper.sh``` launch script. It is the easiest way to start meta-sweeper, providing a tiny bit of boiler plate to the subordinate workflows. Namely, it obtains the installed path of meta-sweeper and checks that nextflow exists on the path. 

Using regular local processes.
```bash
meta-sweeper.sh -c sweep.config run hic-sweep.nf
```
With Nextflow it is easy to submit the work to a grid architecture. For Meta-sweeper, these details are organised as execution profiles in the file named ```execution.config```. We have defined a few examples which may work out of the box on your system.

SGE execution, where the target queue is ```all.q```.
```bash
meta-sweeper.sh -c sweep.config run hic-sweep.nf --profile sge
```

PBS execution, where the target queue is ```batch```.
```bash
meta-sweeper.sh -c sweep.config run hic-sweep.nf --profile pbs
```

[Darling Lab](http://darlinglab.org/)

[i3 institute, UTS
Sydney](http://www.uts.edu.au/research-and-teaching/our-research/ithree-institute)
