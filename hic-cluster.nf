#!/usr/bin/env nextflow
/**
 * Clustering of HiC derived contig graphs (3C-contig graph)
 *
 * Input to this workflow are sourced from the output directory of the data generation
 * workflow -- hic-sweep.nf.
 *
 * As the workflow is much simpler than the primary workflow,
 * channels here are generated in a simpler fashion than with the sweep.
 *
 * A complication is the joining of truth tables and contig graphs into a
 * single channel. This simplifies downstream processes which depend on
 * both results; such as external scoring metrics./
 *
 * Usage: hic-cluster.nf [--debug]
 */
/*
 * meta-sweeper - for performing parametric sweeps of simulated
 * metagenomic sequencing experiments.
 * Copyright (C) 2016 "Matthew Z DeMaere"
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
import MetaSweeper

MetaSweeper ms = MetaSweeper.fromFile(new File('hic.yaml'))

sweep = MetaSweeper.createSweep()

// fetch truth tables from sweep output directory
truths = ms.keyedFrom(sweep, file("${ms.options.output}/*.truth"))

// fetch contig graphs from output dir and pair with respective truth table
graphs = ms.keyedFrom(sweep, file("${ms.options.output}/*.graphml"))
graphs = sweep.joinChannels(graphs, truths, 3)

//
// Perform Clustering
//

// prepare input channels, one per clustering algorithm
(cl_in, graphs) = graphs.into(2)

cl_in = sweep.withVariable('algo', ['louvsoft', 'louvhard', 'oclustr'])
            .describe('Clustering Algorithms')
            .extendChannel(cl_in, 'algo')

chanMap = sweep.forkOnVariable(cl_in, 'algo')

//
// Cluster using Louvain-soft
//

process LouvSoft {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('g.graphml'), truth, algo from chanMap['louvsoft']

    output:
    set key, file("${key}.cl"), truth into ls_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.cl
        """
    }
    else {
        """
        louvain_cluster.py --otype soft --ofmt mcl g.graphml ${key}.cl
        """
    }
}

//
// Cluster using Louvain-hard
//

process LouvHard {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('g.graphml'), truth, algo from chanMap['louvhard']

    output:
    set key, file("${key}.cl"), truth into lh_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.cl
        """
    }
    else {
        """
        louvain_cluster.py --otype hard --ofmt mcl g.graphml ${key}.cl
        """
    }
}


//
// Cluster using OClustR
//

process Oclustr {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('g.graphml'), truth, algo from chanMap['oclustr']

    output:
    set key, file("${key}.cl"), truth into oc_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.cl
        """
    }
    else {
        """
        oclustr.py -f mcl g.graphml ${key}.cl
        """
    }
}

//
// Compute graph statistics
//

// copy channel
(ginfo_in, graphs) = graphs.into(2)

process GraphInfo {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('g.graphml') from ginfo_in

    output:
    set file("${key}.gstat"), file("${key}.geigh") into ginfo_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.gstat
        echo $key > ${key}.geigh
        """
    }
    else {
        """
        graph_stats.py --ofmt yaml g.graphml ${key}.gstat
        graph_complexity.py --ofmt yaml --method eigh g.graphml ${key}.geigh
        """
    }
}


//
// Compute BCubed for reach clustering against its ground truth
//

// Collect all clustering results into a single channel
bc_in = ls_out.mix(lh_out, oc_out)

process Bcubed  {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('clust'), file('truth') from bc_in

    output:
    set key, file("${key}.bc") into bc_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.bc
        """
    }
    else {
        """
        bcubed.py --weighted --tfmt yaml --ofmt yaml truth clust ${key}.bc
        """
    }
}


//
// Compute assembly statistics
//

// fetch contig fastas from output dir
contigs = ms.keyedFrom(sweep, file("${ms.options.output}/*.contigs.fasta"), 2)

process AssemblyStats {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('contigs.fa') from contigs

    output:
    set key, file("${key}.asmstat") into astat_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.asmstat
        """
    }
    else {
        """
        calc_N50_L50.py --ofmt yaml contigs.fa ${key}.asmstat
        """
    }
}
