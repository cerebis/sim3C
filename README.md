Meta-Sweeper
============

By expressing microbial community composition and the details of metagenomic experiments in parametric form, meta-sweeper aims to permit the assessment of analytical methods under variation in both community and experiment.

Prerequisites
-------------

### Nextflow 

- Requires **Java SE 8** or later.
  - due to Groovy's use of java.util.function.BiFunction
- Workflows require a *Linux x86-64* runtime environment 
  - only external binary tools for this architecture are supplied

### Biopython from Pip

- GCC
- Python libraries and header files needed for Python development
  - Package for Debian: python-dev
  - Package for Redhat: python-devel

### PySAM from Pip

- zlib library and headers
  - Package for Debian: zlib1g-dev
  - Package for Redhat: zlib-devel
  
### Beagle Library

- Timeseries workflow requires beagle-lib
  - Requires autoconf tool-chain.
  - Well documented at: [beagle-lib](https://github.com/beagle-dev/beagle-lib/wiki/LinuxInstallInstructions)

Installation
------------

### Nextflow

The [Nextflow](http://www.nextflow.io/) framework, which can be installed easily using either of the following:

```wget -qO- get.nextflow.io | bash``` or ```curl -fsSL get.nextflow.io | bash```

assuming you have met the prerequisite of installing Java 7+. Note, you will need to have ```wget``` or ```curl``` installed depending on your choice above. Whether either is installed by default is dependent on which distribution of Linux you are using.

In addition, meta-sweeper expects that the main executable ```nextflow``` is accessible on the path. Users can move this file to a location already on the path or add its parent directory to the path.

### Python modules

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

#### On ubuntu
```bash
sudo apt-get install python-biopython python-pandas python-yaml python-networkx python-pysam
```
#### On Redhat
```bash
sudo yum install python-biopython python-pandas python-yaml python-networkx python-pysam
```
A more general and complete solution is to use [pip](https://pip.pypa.io/en/stable/).

Using '--upgrade' will ensure that the current version of each module is installed. Though we have not encountered this problem, a note of caution. Upgrading Python modules through Pip, and thus outside of your system's package manager, does incur the potential risk of version conflicts if your installed system packages depend on obsolete or otherwise deprecated functionality with a given Python module. If in doubt, you could try the same command below but omit the '--upgrade' option.

#### Installing requirements with pip
```bash
pip install --upgrade -r requirements.txt
```

Setup
-----

### Configuration Steps

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
    
    For some Linux distributions, beagle-lib can be satisfied through the system package manager. In other cases, users will need to download and [install beagle-lib](https://github.com/beagle-dev/beagle-lib/wiki/LinuxInstallInstructions) manaully. Please be certain to their documentation and make sure that all listed prerequisites described therein are met prior to attempting compilation.

### Automated Setup

We have provided a Bash script which attempts to automate and verify the setup process. __Please Note__: as we are initialising environmental variables for your current shell, you must source this script (```. bash_configure``` OR ```source bash_configure```). Do not instead make it executable, as running it (```./bash_configure```) will not result in changes to your shell's environment.

#### Sourcing setup script
```bash
. bash_configure
```

Users should pay attention to the output from the script. We attempt to highlight problems in red.


Workflow Invocation
-------------------

After initialisation, users may start sweeps using standard Nextflow invocation syntax or simply execute the chosen workflow scrip itself. Users will most likely wish to modify sweep to suit their own purposes, [see below](#sweep-definition) for an explanation of how a sweep is defined.

### Running a workflow

Any workflow can be started either using the following syntax:
```bash
nextflow run [options] <workflow>
```

For instance, running Metagenomic-Hic stage 1 can be like so:
```bash
nextflow run hic-sweep.nf
```

or even more simply, by treating any workflows as an executable:
```bash
./hic-sweep.nf
```

It is worth mentioning here that not all of our workflow scripts below are completely independent. For instance, Metagenomic-HiC is broken into three stages, with the second and third stages depending on the results of the previous stage. Therefore users must begin with stage 1 in this case.  

**Note:** the nature of mixing concurrency and potentially resource hungry processes (such as genome assembly) can mean that a simple local execution strategy may result in resource starvation and subsequently premature program termination. It is recommended, in the long run, that users make use of a [Nextflow supported distributed resource manager (DRM)](https://www.nextflow.io/docs/latest/executor.html) such as SGE, SLURM, etc. 

When a DRM is available as a [target of execution](#execution-targets), workflows can be easily submitted with only minimal configurational change to execute on potentially many physical nodes. Most frequently, this entails changes to the queue name for the target system. Even with only a single physical machine however, a well configured DRM can take better responsibility for managing available compute resources in completing a workflow's set of tasks.

### Common options

Four options which users may frequently employ are (note: single hyphen) as follows.

#### ```-profile [name]``` 

Shorthand for including additional configuration details at runtime (perhaps when using different execution environments) and predefined in a nextflow configuration file under a given label **[name]**. For simplicity, these can be placed in the default nextflow.config and we have provided a few examples for a few different execution targets [described below](#execution-target).
    
#### ```-resume``` 

For many reasons, long running workflows may be interrupted either inadvertently or intentionally. For convenience and potentially a large time-saver, it is possible to restart or resume such jobs.
 
#### ```-queue-size [Int]```
 
Limit concurrency to an integer number **[Int]** of simultaneous running tasks.
 
#### ```-work-dir```
 
Specify a working directory other than the default ```work```. This can he useful when separately running multiple workflows and you wish to inspect the contents of working directories without having to also determine which sub-folders pertain to which workflow.

###  Further command-line options

Nextflow will describe its' base interface by invoking it alone on the command-line:
```bash
nextflow
``` 

Detailed help with each command can be accessed as follows: 
```bash
nextflow help [command]
```

### Execution Targets

For even a relatively shallow parametric sweep, the sum total of computational resources required can be significant. Scheduling systems (SGE, PBS, SLURM, etc) are ideally suited to managing this task and Nextflow makes directing execution to them easy. 

Specific scheduler submission details vary, as can the resources required for individual tasks. This information (executor, queue name, cpu, memory, disk) can be encapsulated either broadly in a [config profile](https://www.nextflow.io/docs/latest/config.html#config-profiles) or per-process using [process directives](https://www.nextflow.io/docs/latest/process.html#directives).

We have provided a few simple examples of such profiles within the default nextflow configuration file [nextflow.config](#predefined-profiles), which is automatically sourced by Nextflow at invocation.

#### Predefined profiles

We have attempted to use the simplest DRM-appropriate defaults, but queue names in particular frequently vary between deployments. Additionally, it is possible to specify [many more criteria](https://www.nextflow.io/docs/latest/config.html#scope-executor) for your particular queuing system.

```
profiles {
        
        standard {
                process.executor = 'local'
        }

        sge {
                process.executor = 'sge'
                queue = 'all.q'
                clusterOptions = '-S /bin/bash'
        }

        pbs {
                process.executor = 'pbs'
                queue = 'batch'
                clusterOptions = '-S /bin/bash'
        }

        pbspro {
                process.executor = 'pbspro'
                queue = 'workq'
                clusterOptions = '-S /bin/bash'
        }
}
```

#### Submission Examples

__Submit to SGE queue manager__
```bash
./hic-sweep.nf -profile sge
```

__Submit to a PBS queue manager__
```bash
./hic-sweep.nf -profile pbs
```

### Sweep Definition

Primarily, a parametric sweep is defined by a set of variables which are sampled over a predefined range and granularity. In addition, there can be many other values held constant throughout the sweep, which are otherwise necessary and relevant to declare.

Effectively, the sweep is a sample-space where it is up to the user to decide on how finely the space is explored, whether the change between points is linear or otherwise. Once decided, the user declares the set of values for each variable in a simple configuration file. 

Within the sweep configuration file, all varied parameters are grouped under the '*variables*' label, while the remaining fixed values are collected under the '*options*' label. 

For the workflows implemented below, we have opted for a separate file for each worflow. Although not be strictly necessary, it permits keeping the size and complexity of the configuration file to a minimum, while more importantly avoiding unintended side-effects when changes with regard to one workflow are inadvertently taken up in another.
  
To tailor a sweep your preferences, users simply to extend or reduce the number of values taken on by any parameter. Note that at least one value must be defined for any variable used in a workflow. By defining a single value, a parameter is effectively fixed in the sweep. Ambition can quickly get the best of you and users should keep in mind that fine sampling of even a few parameters can lead to an exponential explosion in the full parameter space, potentially outstripping the computational resources at hand.

For any workflow, replicates are easily generated by defining multiple seeds.

#### Sweep variables 

How parameters vary in a sweep are defined in the [configuration file](#configuration-file-example). For the implemented workflows, there are a few more common parameters which we will now mention.

+ ```variables:``` This top most label is where swept parameters are defined. **Please note**: unfortunately, it is not possible to create new variables and have them dynamically added to a workflow. Only the values taken on in each case  may be changed.
    
    + ```seed:``` Where possible, all executables which rely on random seeds receive this value. Defining multiple seeds will result in replicates of the entire sweep.
    
    + ```community:``` The critical entry defining the community to be analyzed. Here, a community is compose of an arbitrary number of phylogenetic clades. Each clade may use independent ancestral and donar (LGT) sequences. The per-clade phylogenetic trees are generated by birth-death simulation and pre-clade abundance profiles modeled as log-normal distributions. The size of each clade (number of taxa), the birth and death rates and log-normal parameters *µ* and *σ* are all user accessible. When multiple clades are defined, a final step merges them into a single normalized community.
     
    + ```alpha:``` The branch length scale factor *α<sub>BL</sub>* controls the degree of evolutionary divergence between related taxa, with an recommended range of [0..1]. For a semi-qualitative perspective, evolutionary divergence measured in *ANI<sub>b</sub>* (average nucleotide identity by BLAST alignment) for an isotropic star topology will see *ANI<sub>b</sub>* < 85% (species divergence) when *α<sub>BL</sub>* ≅ 1 and *ANI<sub>b</sub>* > 95% (strain divergence) when *α<sub>BL</sub>* < 0.2 (Figure 1).
    
    + ```xfold:``` Coverage, as measured in "times-genome-size" and used for workflows which make reference to conventional whole-genome shotgun sequencing of a community. For a community of 50 Mbp, xfold=10 would mean 500 Mbp of total sequencing or ~ 1.6 million 150 bp paired-end reads.
    
    + ```n3c:``` For those workflows which deal with HiC sequencing, this variable defines the number of HiC read-pairs generated from a given community. In the case of 3C-contig graphs, 100,000 to 1 million pairs may prove sufficient for a simple community of four bacterial species. For larger or more complex communities, resulting in more fragmented WGS assemblies, extending as sweep beyond 1 million pairs is recommended.

![Alpha vs ANI][fig1-ani]

<p>
<img src="/doc/ani.png" alt="Figure 1: ANIb vs Alpha" height="364" width="450">
<strong>Figure 1</strong>. Average nucleotide identity as a function of scale factor <strong>α<sub>BL</sub></strong>
</p>

#### Configuration file example

Taken from Metagenomic-HiC workflow, the following is an example of a sweep configuration file in YAML syntax.

```yaml
variables:
  # Level 0
  # A complete sweep is run for each random seed defined here. In this case
  # we will have 3 replicates of the entire sweep.
  seed: [1, 2, 3]
  # Level 1
  # Communities are generated for each seed.
  community: !com
    name: trial
    # A community of two clades. 
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
  # Level 2
  # Phylogenetic scale factor [0..1]. 
  # The smaller this value, the more closely related
  # taxa in each clade become.
  alpha: [1, 0.5, 0.1]
  # Level 3
  # WGS sequencing depth, measured in times coverage
  xfold: [1, 5, 10]
  # Level 4
  # The number of HiC read-pairs to generate 
  n3c: [50000, 100000, 1000000]

options:
  # How many samples. For most workflows, this is set to 1.
  num_samples: 1
  # Probabilties supplied to sgEvolver when generating community sequences.
  evo:
    indel_freq: 1e-4
    small_ht_freq: 1e-4
    large_ht_freq: 1e-4
    inversion_freq: 1e-4
  # WGS sequencing.
  wgs:
    read_len: 150
    ins_len: 450
    ins_std: 100
  # 3C/HiC sequencing 
  n3c:
    inter_prob: 0.9
    read_len: 150
  # The output directory
  output: out
```

Implemented Workflows
---------------------

In an effort to highlight meta-sweeper's utility, three workflows have been implemented covering different analysis topics, applicable to studied by parametric sweep. Each can be adjusted through their YAML configuration, and therefore can be tailored to user requirements.

We encourage users to modify these examples for their own purposes.

### 1. Metagenomic HiC

This topic was our original motivation for creating meta-sweeper. The work culminated in the publication [**Deconvoluting simulated metagenomes: the performance of hard- and soft- clustering algorithms applied to metagenomic chromosome conformation capture (3C)**](https://doi.org/10.7717/peerj.2676) and meta-sweeper is intended to allow for straightforward reproduction of that work.


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

#### Explanation of results
    
The resulting aggreation file ```all_stats.yaml``` lists both the input parameters  and their resulting outcomes one sweep point per line. The structure of each line is a map of maps (Java/Groovy) or dict of dicts (Python), which can be easily deserialized. We chose not to create a flat table by default as there is a loss of context in doing so. 

Users can create a flat CSV table using the supplied script ```tabular_results.py```, which can then be inspected in your chosen software (Excel, R, etc). The results, where each row refers to the input parameters and results for a single point in the sweep.

The serialized YAML object contains the following entries:
    
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

### 2. Time-series Deconvolution


* * *

[Darling Lab](http://darlinglab.org/)

[i3 institute, UTS
Sydney](http://www.uts.edu.au/research-and-teaching/our-research/ithree-institute)
