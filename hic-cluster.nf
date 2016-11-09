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

/**
 * Workflow to cluster the results of the primary hic sweep (hic-sweep.nf)
 *
 * Input to this workflow are sourced from the output directory of the primary
 * sweep workflow. As the workflow is much simpler than the primary workflow,
 * channels here are generated in a simpler fashion than with the sweep.
 *
 * A complication is the joining of truth tables and contig graphs into a
 * single channel. This simplifies downstream processes which depend on
 * both results; such as external scoring metrics.
 */
import MetaSweeper

MetaSweeper ms = MetaSweeper.fromFile(new File('hic.yaml'))

// fetch truth tables from sweep output directory
truths = ms.keyedFrom(file("${ms.options.output}/*.truth"))

// fetch contig graphs from output dir and pair with respective truth table
graphs = ms.keyedFrom(file("${ms.options.output}/*.graphml"))
graphs = ms.joinChannels(graphs, truths, 3)

/**
 * Cluster using louvain-soft
 */
// prepare input channel
(ls_in, graphs) = graphs.into(2)

        // extend the key to include cluster algorithm.
ls_in = ls_in.map { [ ms.extendKey(it.getKey(), 'algo', 'louvsoft'), *it.dropKey() ] }

process LouvSoft {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('g.graphml'), truth from ls_in

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

/**
 * Cluster using louvain-hard
 */
// prepare input channel
(lh_in, graphs) = graphs.into(2)

        // extend the key to include cluster algorithm.
lh_in = lh_in.map { [ ms.extendKey(it.getKey(), 'algo', 'louvhard'), *it.dropKey() ] }

process LouvHard {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('g.graphml'), truth from lh_in

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

/**
 * Cluster using OClustR
 */
// prepare input channel
(oc_in, graphs) = graphs.into(2)

        // extend the key to include cluster algorithm.
oc_in = oc_in.map { [ ms.extendKey(it.getKey(), 'algo', 'ocluster'), *it.dropKey() ] }

process Oclustr {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('g.graphml'), truth from oc_in

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

/**
 * Compute graph statistics
 */
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

/**
 * Compute BCubed for reach clustering against its ground truth
 */

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
        bcubed.py --tfmt yaml --ofmt yaml truth clust ${key}.bc
        """
    }
}

/**
 * Compute assembly statistics
 */

        // fetch contig fastas from output dir
contigs = ms.keyedFrom(file("${ms.options.output}/*.contigs.fasta"), 2)

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
