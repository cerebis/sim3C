Meta-Sweeper
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

### Basic Usage

**Environmental Configuration**

Before running a meta-sweeper workflow, you must initialise the shell environment. 

1. Set meta-sweepers home to your installed location.
  - If meta-sweeper was installed in /home/user/meta-sweeper.
  - ```export METASWEEPER_HOME=/home/user/meta-sweeper```
2. Append the meta-sweeper home patht to nextflow's classpath.
  - ```export NXF_CLASSPATH=$NXF_CLASSPATH:$METASWEEPER_HOME```
3. Check that nextflow is installed and has been added to the path.
  - Either of the following should return the canonical path to the nextflow executable.
  - ```command -v nextflow``` or ```which nextflow```
4. A further dependency on [beagle-lib](https://github.com/beagle-dev/beagle-lib) exists for the [timeseries](#2-time-series-deconvolution) workflow.
    
    To use this workflow, users must set the environmental variable BEAGLE_LIB to point to the directory containing ```libhmsbeagle-jni.so``` shared library file. 
    
    E.g. ```BEAGLE_LIB=/usr/lib```
    
    For some Linux distributions, beagle-lib can be satisfied through the system package manager. In other cases, users will need to download and [install beagle-lib](https://github.com/beagle-dev/beagle-lib/wiki/LinuxInstallInstructions). Please be certain that all prerequisites described therein are met prior to attempting compilation.

We have provided a Bash script which attempts to automate and verify this process. Sourcing this script is done as follows. __Note__: as we are trying to initialise environmental variables in your shell, you must source this script (```. bash_configure``` OR ```source bash_configure```) rather than execute it (```./bash_configure```). Execution will prevent environmental variables from being set in your shell.

```bash
. bash_configure
```

Users should pay attention to the output from the script. Problems will be highlighted in red.


#### Invocation

After initialisation, users may start sweeps using standard Nextflow invocation syntax.

[See below](#sweep-definition) for an explanation of how a sweep is defined.

__Default execution__

Any workflow can be started either explicitly via nextflow's syntax:

```bash
nextflow run hic-sweep.nf
```

or by treating any of the workflows as executables:


```bash
./hic-sweep.nf
```

**Note:** the nature of mixing concurrency and potentially resource hungry processes (such as genome assembly) can mean that a basic local execution strategy may result in resource stravation and subsequently premature program termination. It is recommended that, in the long run, it is worthwhile for users to configure a [Nextflow supported distributed resource manager (DRM)](https://www.nextflow.io/docs/latest/executor.html) such as SGE, SLURM, etc. to take responsibility for managing available local resources.

__Command-line options__

Nextflow's base command structure can be seen by invoking ```nextflow``` with no options. While help with each internal command can be found by ```nextflow help [command]```.

When invoking the ```nextflow run```, two important options are ```-profile``` and ```-resume```. The former permits a shorthand for including additional configuration details at runtime (perhaps when using different execution environments), while the latter informs Nextflow attempt and resume an interrupted workflow.

__Execution Target__

For even a relatively shallow parametric sweep, the sum total of computational resources required can be significant. Scheduling systems (SGE, PBS, SLURM, etc) are ideally suited to managing this task and Nextflow makes directing execution to them easy. 

Specific scheduler submission details vary, as can the resources required for individual tasks. This information (executor, queue name, cpu, memory, disk) can be encapsulated either broadly in a [Config profile](https://www.nextflow.io/docs/latest/config.html#config-profiles) or per-Process using [Process directives](https://www.nextflow.io/docs/latest/process.html#directives).

We have provided a few simple examples of such profiles within the default nextflow configuration file ```nextflow.config```, which unless overridden, is automatically sourced by Nextflow at invocation.

__Submission Examples__

Submit to SGE queue manager.

```bash
./hic-sweep.nf -profile sge
```

Submit to a PBS queue manager.

```bash
./hic-sweep.nf -profile pbs
```

### Implemented Workflows

In an effort to highlight meta-sweeper's utility, three workflows have been implemented covering different analysis topics, applicable to studied by parametric sweep. Each can be adjusted through their YAML configuration, and therefore can be tailored to user requirements.

We encourage users to modify these examples for their own purposes.

#### Sweep Definition

For each workflow that follows, a configuration file exists which permits users to modify how the parameters involved vary. This text file is easy to understand and follows YAML syntax. There is a separate configuration file for each workflow, since they do not all share the same parameter set.
  
Users are free to extend or reduce the number of values taken on by any parameter but must define at least one value. Defining a single value essentially fixes a parameter in the sweep. Users should be aware that excessively fine sampling of even a few parameters can lead to an exponential explosion in the full parameter space, potentially outstripping the computational resources at hand.

Replicates can be performed by defining more than one seed value.

####1. Metagenomic HiC

This topic was our original motivation for creating meta-sweeper. The work culminated in the publication [**Deconvoluting simulated metagenomes: the performance of hard- and soft- clustering algorithms applied to metagenomic chromosome conformation capture (3C)**](https://doi.org/10.7717/peerj.2676) and meta-sweeper is intended to allow for straightforward reproduction of that work.

__Configuration and sweep definition__ 

The sweep and how parameters are varied are defined in the configuration file. 

+ Configuration file: *hic.yaml*

The composition of a community is defined under the ```community``` label, where each ```clade``` represents a single explicit phylum, with independent ancestral (VGT), and donar (HGT) sequences. The random seed is used to generate each clade's phylogenetic tree (birth-death) and abundance profile (log-norm). Users can set birth and death rates for trees, as well as log-normal parameters *µ* and *σ* for abundance profiles.  

```yaml
variables:
  # level 0
  # A complete sweep is run for each random seed
  seed: [1, 2, 3]
  # level 1
  # communities are generated for each seed.
  community: !com
    name: trial
    # This community has two clades. 
    # Clades can have any ancestral/donar sequence, but here they use the same.
    # Altogether there will be 8 taxa in the community.
    clades:
      - !clade
        prefix: clade1
        ancestor: test/ancestor.fa
        donor: test/donor.fa
        ntaxa: 5
        tree: {birth: 0.9, death: 0.5}
        profile: {mu: 0.1, sigma: 1}
      - !clade
        prefix: clade2
        ancestor: test/ancestor.fa
        donor: test/donor.fa
        ntaxa: 3
        tree: {birth: 0.7, death: 0.3}
        profile: {mu: 0.2, sigma: 1.5}
  # level 2
  # phylogenetic scale factor [0..1]. 
  # the smaller this value, the more closely related each
  # taxon in a clade become (shorter branches)
  alpha: [1, 0.5]
  # level 3
  # WGS sequencing depth, measured in times coverage
  xfold: [1, 2]
  # level 4
  # The number of HiC read-pairs to generate 
  n3c: [5000, 10000]
options:
  # how many samples. Non-timeseries, this is fixed to one
  num_samples: 1
  # probabiltiy factors for SGEvoler stage.
  evo:
    indel_freq: 1e-4
    small_ht_freq: 1e-4
    large_ht_freq: 1e-4
    inversion_freq: 1e-4
  # WGS sequencing factors.
  wgs:
    read_len: 150
    ins_len: 450
    ins_std: 100
  # 3C/HiC sequencing factors
  n3c:
    inter_prob: 0.9
    read_len: 150
  # the outut directory
  output: out
```

The complete workflow is actually broken into three smaller stages:

1. __Data Generation__ 
    
+ Script: *hic-sweep.nf*

    Creation of communities, WGS and HiC read simulation, Metagenome assembly and read mapping. How each each parameter should vary within the sweep can be adjusted in the configuration file.

2. __Clustering__

+ Script: *hic-cluster.nf*

   After data generation, for each sample point within the sweep, 3C-contig clustering is preformed by Louvain-hard, Louvain-soft and OClustR algorithms. Afterwards, performance and quality metrics are applied. The BCubed external metric is used to assess the performance of each algorithm relative to the ground truth, while simple assembly (N50, L50) and graph (size, order) statistics are compiled alongside an entropic measure of graphical complexity (H<sub>L</sub>) 

3. __Aggregation__

+ Script: *hic-aggregate.nf*

   This is the simplest stage. Here the results from potentially many permuted outcomes are collected and aggregated into a single results file *all_stats.yaml*. The resulting text file in YAML syntax is structured as an array of associative collections, one per sweep point. Validation results are grouped by the queried target: whether that be the assembly, graph, clustering.
   
   The file can be easily deserialized to an object within any language where YAML support exists, which is widely avaiable: [Python](pyyaml.org), [Java & Groovy](www.snakeyaml.org), [Ruby](https://ruby-doc.org/stdlib-1.9.3/libdoc/yaml/rdoc/YAML.html), [C++](https://github.com/jbeder/yaml-cpp), etc. 

   __Example of all_stats.yaml__
    
    A single line from the resulting aggreation file, contains the following entries:
    
    - *params* - the set of input sweep parameters
    - *asmstat* - common assembly statistics
    - *bc* - the weighted BCubed external index (see manuscript)
    - *geigh* - a entropic measure of graphical complexity (see manuscript)
    - *gstat* - common graph statistics

```yaml
params: {seed: '2', alpha: '1', xfold: '1', n3c: '5000', algo: louvsoft}
asmstat: {L50: 12, N50: 455}
bc: {completeness: 0.964, f: 0.141, pre: 1.0, rec: 0.0757}
geigh: {method: eigh, value: 0.1}
gstat: {density: 0.074, inclusives: 28, isolates: 0, mean_deg: 2.0, median_deg: 2, modularity: 0.964, order: 28, size: 28}
```


####2. Time-series Deconvolution

####3. Euk 

* * *

[Darling Lab](http://darlinglab.org/)

[i3 institute, UTS
Sydney](http://www.uts.edu.au/research-and-teaching/our-research/ithree-institute)
