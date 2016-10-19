meta-sweeper
============

By expressing microbial community composition and the details of metagenomic experiments in parametric form, meta-sweeper aims to permit the assessment of analytical methods under variation in both community and experiment.

Prerequisites
-------------

### Nextflow 
- Presently requires **Oracle** Java 7+. 
- Workflows require a *Linux x86-64* runtime environment due to dependence on external statically linked tools.

### Biopython from Pip
- GCC or equivalent
- Python libraries and header files needed for Python development (Debian: python-dev, Redhat: python-devel)

### Pysam from Pip
- zlib library and headers (Debian: zlib1g-dev, Redhat: zlib-devel)

Installation
------------

### Nextflow

The [Nextflow](http://www.nextflow.io/) framework, which can be installed easily using either of the following:

```wget -qO- get.nextflow.io | bash``` or ```curl -fsSL get.nextflow.io | bash```

assuming you have met the prerequisite of installing Java 7+. Note, you will need to have ```wget``` or ```curl``` installed depending on your choice above. Whether either is installed by default is dependent on which distribution of Linux you are using.

In addition, meta-sweeper expects that the main executable ```nextflow``` is accessible on the path. Users can move this file to a location already on the path or add its parent directory to the path.

#### Python modules

The workflows depend on the following Python modules, which must be installed prior to using it:

* biopython
* dendropy
* intervaltree
* networkx
* numpy
* pandas
* python-louvain
* pysam
* PyYaml
* scipy

Some of these dependencies can be installed using distributional package manager, but not all will be found.

On ubuntu
```bash
sudo apt-get install python-biopython python-pandas python-yaml python-networkx python-pysam
```
On Redhat
```bash
sudo yum install python-biopython python-pandas python-yaml python-networkx python-pysam
```
A more general and complete solution is to use [pip](https://pip.pypa.io/en/stable/).

Using '--upgrade' will ensure that the current version of each module is installed. Though we have not encountered this problem, a note of caution. Upgrading Python modules through Pip, and thus outside of your system's package manager, does incur the potential risk of version conflicts if your installed system packages depend on obsolete or otherwise deprecated functionality with a given Python module. If in doubt, you could try the same command below but omit the '--upgrade' option.

```bash
pip install --upgrade -r requirements.txt
```

### Usage

We recommend users start sweeps using the ```meta-sweeper.sh``` launch script. It is the easiest way to start meta-sweeper, providing a tiny bit of boiler plate to the subordinate workflows. Namely, it obtains the installed path of meta-sweeper and checks that nextflow exists on the path. 

How parameters are varied over the sweep are defined in a Nextflow configuration file. An example ```sweep.config``` along with supporting files in ```test``` has been provided.

#### Local execution

Using regular local processes.
```bash
meta-sweeper.sh -c sweep.config run hic-sweep.nf
```

**Note:** the nature of mixing concurrency and potentially resource hungry processes (such as genome assembly) can mean that a basic local execution strategy may result in resource stravation and subsequently premature program termination. It is recommended that, in the long run, it is worthwhile for users to configure a [Nextflow supported distributed resource manager (DRM)](https://www.nextflow.io/docs/latest/executor.html) such as SGE, SLURM, etc. to take responsibility for managing available local resources.

#### Distributed execution

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
