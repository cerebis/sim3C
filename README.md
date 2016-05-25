meta-sweeper
============

By expressing microbial community composition and the details of metagenomic experiments in parametric form, meta-sweeper aims to permit the assessment of analytical methods under variation in both community and experiment.


### Prerequisites
- Nextflow requires Java 7+.
- Workflows require a *Linux x86-64* runtime environment due to dependence on external statically linked tools.

### Dependencies

#### Nextflow

The [Nextflow](http://www.nextflow.io/) framework, which can be installed easily using either of the following:

```wget -qO- get.nextflow.io | bash``` or ```curl -fsSL get.nextflow.io | bash```

assuming you have met the prerequisite of installing Java 7+.

In addition, meta-sweeper expects that the main executable ```nextflow``` is accessible on the path. Users can move this file to a location already on the path or add its parent directory to the path.

#### Python modules

The workflows depend on the following Python modules, which must be installed prior to using it:

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

### Usage

We recommend that users start sweeps by using the ```meta-sweeper.sh``` launch script. It is the easiest way to start meta-sweeper, providing a tiny bit of boiler plate to the subordinate workflows. Namely, it obtains the installed path of meta-sweeper and checks that nextflow exists on the path. 

How parameters are varied over the sweep are defined in a Nextflow configuration file. An example ```sweep.config``` along with supporting files in ```test``` has been provided.

#### Local execution

Using regular local processes.
```bash
meta-sweeper.sh -c sweep.config run hic-sweep.nf
```

#### Grid execution

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
