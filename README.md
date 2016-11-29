Meta-Sweeper
============

### Table of Contents

1. [__Introduction__](#introduction)
2. [__Prerequisites__](#prerequisites)
3. [__Installation__](#installation)
    - [Nextflow](#nextflow)
    - [Python modules](#python-modules)
4. [__Setup__](#setup)
    - [Configuration steps](#configuration-steps)
    - [Automated setup](#automated-setup)
5. [__Workflow Invocation__](#workflow-invocation)
    - [Running a workflow](#running-a-workflow)
    - [Common options](#common-options)
    - [Further options](#further-command-line-options)
    - [Execution targets](#execution-targets)
    - [Predefined profiles](#predefined-profiles)
    - [Submission examples](#submission-examples)
    - [Trouble shooting](#trouble-shooting)
6. [__Sweep Definition__](#sweep-definition)
    - [Sweep variables](#sweep-variables)
    - [Configuration file](#configuration-file-example)
    - [Output files](#output-file-naming)
7. [__Implemented Workflows__](#implemented-workflows)
    - [Metagenomic HiC](#1-metagenomic-hic)
    - [Time-series deconvolution](#2-time-series-deconvolution)
8. [__Included Tools__](#included-tools)
    - [Read simulation](#read-simulation)
    - [Clustering algorithms](#clustering-algorithms)

Introduction
------------

Meta-Sweeper is a tool for conducting parametric (combinatorial) sweeps of selected experimental variables relevant to metagenomics research. In particular, we have implemented workflows which enable the generation of simulated microbial communties and their simulated sequecing (WGS and HiC/3C). Downstream of simulation, specific topics of analysis in HiC/3C clustering and Non-negative matrix factorization have been provided. Built upon the [Nextflow framework](http://nextflow.io), Meta-Sweeper is remains general enough that there is no reason that alternative workflows or extensions to existing workflows could not be implemented.

Meta-Sweeper is the culmination and refinement of a central process employed in our ongoing work exploring how to obtain the maximum value when integrating HiC/3C sequencing data with conventional whole-genome shotgun DNA sequencing. To accomplish this, we conducted a parametric sweep of select experimental and model parameters such as sequencing depth, community composition and evolutionary divergence. At each sampled point in the sweep, we pitted a set of algorithms against each-other, where their aim was (in a sense) the deconvolution of the simulated community. This work culminated in the publication: [**Deconvoluting simulated metagenomes: the performance of hard- and soft- clustering algorithms applied to metagenomic chromosome conformation capture (3C)**](https://doi.org/10.7717/peerj.2676).

Parametric sweeps, and the sampling of parameter spaces in general, have reoccurring applicability within bioinformatics beyond our own topic of research. We believe that this systematic approach is potentially useful to a broad range of research topics, offering a means for more thorough quality assurance and performance testing. Aside from explicit testing, the process of sampling a parameter space can also be used simply as an exploratory technique. In all cases, a reproducible approach is preferable to one-off ad-hoc approaches that can prove difficult to reproduce even for the original authors. 

Some suggested applications are such things as: algorithm development, the exploration of experimental requirements for new techniques or simply the comparative assessment of existing tools. 

Although numerous workflow management systems exist, we feel there is a lack of facility for parametric sweeps. Further, ease of deployment to various high-throughput computing environments is a crucial feature, when bioinformatics tasks are often computationally expensive in terms of CPU and/or memory resources.

To that end, we invested the time to resolve our processes into a more coherent and easily deployed system. We hope this will prove useful to others.

Prerequisites
-------------

- __Python 2.7__
- __Java SE 8+__
  - due to Groovy's use of java.util.function.BiFunction
- __Linux x86-64 runtime environment__ 
  - supplied external binary tools
- __GCC__
  - if building Biopython (via pip)
- __Python development library and headers__
  - Debian package: _python-dev_
  - Redhat package: _python-devel_
- __zlib library and headers__
  - Package for Debian: _zlib1g-dev_
  - Package for Redhat: _zlib-devel_
- __Beagle Library__
  - if using time-series workflow
  - requires autoconf tool-chain.
  - process well documented at: [beagle-lib](https://github.com/beagle-dev/beagle-lib/wiki/LinuxInstallInstructions)

Installation
------------

### Nextflow

The [Nextflow](http://www.nextflow.io/) framework, which can be installed easily using either of the following:

```wget -qO- get.nextflow.io | bash``` __or__ ```curl -fsSL get.nextflow.io | bash```

This assumes you have met the prerequisite of installing Java 8+. Note, you will need to have ```wget``` or ```curl``` installed depending on your choice above. Whether either is installed by default is dependent on which distribution of Linux you are using.

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

### Predefined profiles

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

### Submission Examples

__Submit to SGE queue manager__
```bash
./hic-sweep.nf -profile sge
```

__Submit to a PBS queue manager__
```bash
./hic-sweep.nf -profile pbs
```

### Trouble Shooting

#### ```--debug```

All of our [implemented workflows below](#implemented-workflows) accept the double-hyphen runtime option ```--debug```. 

Including this option at invocation time will test your environmental setup and the workflow processing logic but will not execute any actual bioinformatics tools. Workflows will execute very quickly in this mode and, if all is well, will complete without emitting an error. The mode creates only mock output files to meet inter-process dependencies. This mode is helpful when first configuring your system, particularly when configuring submission to a queue manager.

#### Nextflow log

Nextflow writes to a hidden log file (.nextflow.log), which it manages in a rotation akin to Linux logrotate. If an error occurs during a workflow, it is very likely that Nextflow will expose it as a message to stdout, prompting the user to inspect the log file. In most cases, Nextflow will provide details on which task failed and its working directory. This directory should be the first port of call when troubleshooting a new problem. We would recommend inspecting both command.err (.command.err) and command.out (.command.out) to see if what went wrong can be ascertained. If not, inspect command.sh (.command.sh) and verify that the issued command is properly formed. 

It is our experience that runtime errors are often caused by unset environmental variables. This can be as simple as inadvertantly switching shell terminals, where the second shell has not been configured. Therefore, even when you think you have done so already, please double check that you have initialized the environment as [outlined above](#configuration-steps).

#### Some Error Sources
1. Unset environmental variables
2. Incorrect runtime architecture on execution host (requires x86_64)
3. Incorrect Java JVM. Please try using Java 8 or later.
4. Out of storage space. Sweeps can occupy significant space.
5. Exceeding CPU or other resource limits on execution host.
6. Missing Python dependencies for invoked interpreter.
   
   If using a Python interpreter _other than_ the system default, make sure that the alternate's path _(/home/foobar/bin)_ preceeds that of the system's default _(/usr/bin)_.

    e.g. PATH=/home/foobar/bin:/usr/bin
    
7. No internet access.

    If necessary, Meta-Sweeper uses [Groovy Grape](http://docs.groovy-lang.org/latest/html/documentation/grape.html) to automatically satisfy a few dependencies from internet repositories. Without access, compliation of our workflows will fail (see below).

#### Script Compilation Errors.

In more severe cases, the Nextflow script might itself have failed to compile _(groovy scripts are built at invocation time)_. In this situation, errors can be harder to decrypt and we would encourage users to contact us.

Sweep Definition
----------------

Primarily, a parametric sweep is defined by a set of variables which are sampled over a predefined range and granularity. In addition, there can be many other values held constant throughout the sweep, which are otherwise necessary and relevant to declare.

Effectively, the sweep is a sample-space where it is up to the user to decide on how finely the space is explored, whether the change between points is linear or otherwise. Once decided, the user declares the set of values for each variable in a simple configuration file. 

Within the sweep configuration file, all varied parameters are grouped under the '*variables*' label, while the remaining fixed values are collected under the '*options*' label. 

For the workflows implemented below, we have opted for a separate file for each worflow. Although not be strictly necessary, it permits keeping the size and complexity of the configuration file to a minimum, while more importantly avoiding unintended side-effects when changes with regard to one workflow are inadvertently taken up in another.
  
To tailor a sweep your preferences, users simply to extend or reduce the number of values taken on by any parameter. Note that at least one value must be defined for any variable used in a workflow. By defining a single value, a parameter is effectively fixed in the sweep. Ambition can quickly get the best of you and users should keep in mind that fine sampling of even a few parameters can lead to an exponential explosion in the full parameter space, potentially outstripping the computational resources at hand.

For any workflow, replicates are easily generated by defining multiple seeds.

### Sweep variables 

How parameters vary in a sweep are defined in the [configuration file](#configuration-file-example). For the implemented workflows, there are a few more common parameters which we will now mention.

+ ```variables:``` This top most label is where swept parameters are defined. **Please note**: unfortunately, it is not possible to create new variables and have them dynamically added to a workflow. Only the values taken on in each case  may be changed.
    
    + ```seed:``` Where possible, all executables which rely on random seeds receive this value. Defining multiple seeds will result in replicates of the entire sweep.
    
    + ```community:``` The critical entry defining the community to be analyzed. Here, a community is compose of an arbitrary number of phylogenetic clades. Each clade may use independent ancestral and donar (LGT) sequences. The per-clade phylogenetic trees are generated by birth-death simulation and pre-clade abundance profiles modeled as log-normal distributions. The size of each clade (number of taxa), the birth and death rates and log-normal parameters *µ* and *σ* are all user accessible. When multiple clades are defined, a final step merges them into a single normalized community.
     
    + ```alpha:``` The branch length scale factor *α<sub>BL</sub>* controls the degree of evolutionary divergence between related taxa, with an recommended range of [0..1]. For a semi-qualitative perspective, evolutionary divergence measured in *ANI<sub>b</sub>* (average nucleotide identity by BLAST alignment) for an isotropic star topology will see *ANI<sub>b</sub>* < 85% (species divergence) when *α<sub>BL</sub>* ≅ 1 and *ANI<sub>b</sub>* > 95% (strain divergence) when *α<sub>BL</sub>* < 0.2 (Figure 1).
    
    + ```xfold:``` Coverage, as measured in "times-genome-size" and used for workflows which make reference to conventional whole-genome shotgun sequencing of a community. For a community of 50 Mbp, xfold=10 would mean 500 Mbp of total sequencing or ~ 1.6 million 150 bp paired-end reads.
    
    + ```n3c:``` For those workflows which deal with HiC sequencing, this variable defines the number of HiC read-pairs generated from a given community. In the case of 3C-contig graphs, 100,000 to 1 million pairs may prove sufficient for a simple community of four bacterial species. For larger or more complex communities, resulting in more fragmented WGS assemblies, extending as sweep beyond 1 million pairs is recommended.

<p>
<br>
<img src="/docs/ani.png" alt="Figure 1: ANIb vs Alpha" height="364" width="450">
<br><strong>Figure 1</strong>. Average nucleotide identity as a function of scale factor <strong>α<sub>BL</sub></strong>
</p>

### Configuration file example

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

### Output file naming

At each point in the sweep, files which are considered outputs at various steps within a workflow are copied to the output folder [default: out]. For each such output file, the file name is used to pass-on information about the the parameters used at that given point in the sweep. This information is encoded as a structured string using two delimiters. These delimiters were chosen to avoid conflicts with system constraints, avoid unicode and be at least *somewhat* human readable. 

This string is then prepended to the respective file name, which afterwards follows the regular convention of encoding information about the type of file.

__Name/Value encoding:__

Per-parameter syntax is delimited by a ```#``` character as follows: 
    
    [parameter_name]#[parameter_value]

__Chaining parameters:__

Multiple parameters are joined by the three character string ```-+-```:

    [parameter_name1]#[parameter_value1]-+-[parameter_name2]#[parameter_value2]

__Example:__

As an example, the following would be the file names for WGS reads R1/R2 pertaining to the sweep point seed=1, *α<sub>BL</sub>*=0.5, xfold=10:
    
     seed#1-+-alpha#0.5-+-xfold#10.wgs.r1.fq.gz
     seed#1-+-alpha#0.5-+-xfold#10.wgs.r2.fq.gz


Implemented Workflows
---------------------

In an effort to highlight the utility of Meta-Sweeper, we have implemented three workflows covering different metagenomic analysis topics we feel applicable to study by parametric sweep. In each case, the workflow definition can be adjusted through the respective configuration file, tailoring it to user requirements.

We encourage users to modify these examples for their own purposes.

### 1. Metagenomic-HiC

This topic was our original motivation for creating Meta-Sweeper. The work culminated in the publication [**Deconvoluting simulated metagenomes: the performance of hard- and soft- clustering algorithms applied to metagenomic chromosome conformation capture (3C)**](https://doi.org/10.7717/peerj.2676), where Meta-Sweeper represents a refinement of the methods employed in that work, allowing for its straightforward reproduction.

The sweep is varied over five levels:

1. Random seed
2. Community structure
3. Community divergence (*α<sub>BL</sub>*)
4. WGS coverage (xfold)
5. HiC depth (n3c)

__Sweep file:__ *hic.yaml*

The complete workflow has been broken into three stages, which we feel are often of separate interest. 

#### Stage 1: Data Generation 
    
+ __Script:__ *hic-sweep.nf*

    This first stage is responsible for the creation of communities, WGS and HiC read simulation, metagenome assembly and read mapping. The results from this stage are copied to the output folder [default: out]. How each each parameter should vary within the sweep can be adjusted in the configuration file. We would recommend users consider storage and CPU requirements prior to expanding the size of the sweep.
    
    __Note__: the definition as it stands in only toy-like in size.
    
#### Stage 2: Clustering

+ __Script:__ *hic-cluster.nf*

   After stage 1 has completed, stage 2 out clustering of the 3C-contig graphs resulting from stage 1. At present, clustering is performed with the algorithms: Louvain-hard, Louvain-soft and OClustR. Afterwards, performance and quality metrics are applied. The BCubed external metric is used to assess the performance of each algorithm relative to the ground truth, while simple statistics for assembly (N50, L50) and graphs (size, order) are compiled alongside an entropic measure of graphical complexity (H<sub>L</sub>) 

#### Stage 3: Aggregation

+ __Script:__ *hic-aggregate.nf*

   Acting effectively as a *reduce* function, this is the simplest stage. Here the results from potentially many permuted outcomes are collected and aggregated into a single results file *all_stats.yaml*.
   
   This results file (e.g. ```out/all_stats.yaml```) is organized as a list, where each element encompasses both the input parameters and resulting outcomes from one sweep point. The elements are stored as associative collections, one per sweep point, where validation results from stage 2 are grouped by the statistical test performed. 
   
   The file can be deserialized to an object within any language where support for YAML exists (i.e. [Python](pyyaml.org), [Java & Groovy](www.snakeyaml.org), [Ruby](https://ruby-doc.org/stdlib-1.9.3/libdoc/yaml/rdoc/YAML.html), [C++](https://github.com/jbeder/yaml-cpp), etc.) If users do not wish to deal with a serialized object, the result file can be converted to a flat table using the script ```tabular_results.py```. When converted to tabular form, rows represent the results for each point in the sweep.

#### Results file definition

The serialized YAML object contains the following entries:
    
- *params* - the set of input sweep parameters
- *asmstat* - common assembly statistics
- *bc* - the weighted BCubed external index (see manuscript)
- *geigh* - a entropic measure of graphical complexity (see manuscript)
- *gstat* - common graph statistics

Example: two sweep points from all_stats.yaml
```yaml
- params: {seed: '2', alpha: '1', xfold: '1', n3c: '5000', algo: louvsoft}
  asmstat: {L50: 12, N50: 455}
  bc: {completeness: 0.964, f: 0.141, pre: 1.0, rec: 0.0757}
  geigh: {method: eigh, value: null}
  gstat: {density: 0.0741, inclusives: 28, isolates: 0, mean_deg: 2.0, median_deg: 2, modularity: 0.964, order: 28, size: 28}
- params: {seed: '3', alpha: '1', xfold: '2', n3c: '5000', algo: louvsoft}
  asmstat: {L50: 130, N50: 568}
  bc: {completeness: 0.953, f: 0.0152, pre: 1.0, rec: 0.008}
  geigh: {method: eigh, value: -0.0}
  gstat: {density: 0.006, inclusives: 342, isolates: 0, mean_deg: 2.006, median_deg: 2, modularity: 0.997, order: 342, size: 343}
```

which when converted to tabular format would look like:

```csv
,params-xfold,params-alpha,params-n3c,params-seed,params-algo,asmstat-L50,asmstat-N50,gstat-mean_deg,gstat-density,gstat-modularity,gstat-median_deg,gstat-inclusives,gstat-isolates,gstat-order,gstat-size,geigh-method,geigh-value,bc-pre,bc-rec,bc-completeness,bc-f
0,1,1,5000,2,louvsoft,12,455,2.0,0.0741,0.964,2,28,0,28,28,eigh,,1.0,0.0757,0.964,0.141
1,2,1,5000,3,louvsoft,130,568,2.006,0.006,0.997,2,342,0,342,343,eigh,-0.0,1.0,0.008,0.953,0.0152
```

### 2. Time-series Deconvolution

A test case for metagenome deconvolution by non-negative matrix factorisation of a time-series data-set.

This workflow is implemented as a single script. Internally, the stages of execution involve:
 
1. Community creation
2. WGS read simulation
3. WGS read mapping
4. NNMF based deconvolution
 
The sweep is varied over four levels:

1. Random seed
2. Community structure
3. Community divergence (*α<sub>BL</sub>*)
4. WGS coverage (xfold)

This workflow does not involve HiC sequencing data.

__Script:__ *timeseries-deconvolute.nf*

__Sweep file:__ *timeseries.yaml*

The final outcome of the deconvolution process for each sweep point can be found in the files following the naming strategy: ${key}.truth.report.txt. Here, the variable ${key} represents the structured string encoding information about the parameters used at a given point in the sweep. [This is explained above](#output-file-naming)
 
 E.g. For seed=1, *α<sub>BL</sub>*=0.5 and xfold=10, the report file would be named:
    
     seed#1-+-alpha#0.5-+-xfold#10.truth.report.txt


Included Tools
--------------

### Read Simulation

- __metaART –__ Simulation of whole-genome shotgun paired-end reads from communities. 
  
  metaART simply wraps the art_illumina binary from ART ([Huang et al, 2012](https://doi.org/10.1093/bioinformatics/btr708)), which among many machine types and sequencing modes, simulates Illumina paired-end reads. Our extension allows users to associate relative abundances with the supplied reference sequences. Profiles can be either supplied in the form of a pre-cast table or determined at runtime from a chosen distribution. The association of abundance profiles with the reference sequences, permits the simulated sampling of metagenomic communities or multi-chromosomal clonal sequencing projects with variable copy-number.

- __sim3C –__ Simulation of HiC/3C read-pairs.
Art module -- full read error model. noise. duplication of sites, etc.

### Clustering Algorithms
- __louvain_cluster –__ Louvain community detection based graph clustering. Both traditional hard-clustering solutions of the best partition and a naive soft-clustering option ([DeMaere and Darling, 2016](https://doi.org/10.7717/peerj.2676)).

- __oclustr –__ An implementation of the OClustR graph hard-clustering algorithm ([Pérez-Suárez et al, 2013](http://dx.doi.org/10.1016/j.neucom.2013.04.025)).




* * *

[Darling Lab](http://darlinglab.org/)

[i3 institute, UTS
Sydney](http://www.uts.edu.au/research-and-teaching/our-research/ithree-institute)
