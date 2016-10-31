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
 * channels here are generated from conventional Nextflow methods, rather than
 * reliance on the MetaSweeper classes.
 *
 * The only complication is the initial joining of truth tables and contig graphs
 * into a single channel. This simplifies downstream processes which depend on
 * both results; such as external scoring metrics.
 */
import MetaSweeper

MetaSweeper ms = MetaSweeper.fromFile(new File('sweep.yaml'))

// fetch truth tables from sweep output directory
truths = Channel.from(file("${ms.options.output}/*.truth"))
                // add a string key for each iteration step
                .map { f -> [MetaSweeper.dropSuffix(f.name), f] }

// fetch contig graphs from output dir and pair with respective truth table
graphs = Channel.from(file("${ms.options.output}/*.graphml"))
                // add a string key for each iteration step
                .map { f -> [MetaSweeper.dropSuffix(f.name), f] }
                // prepare a key relevant to truth table from clustering (one less iterative depth)
                .map { t -> [MetaSweeper.removeLevels(t[0], 1), *t ] }
                // join channels, using new key as index
                .phase(truths)
                // unwrap lists and remove redundant elements
                .map{ t -> [*t[0][1..-1], t[1][1]] }


// copy channel
(cl_in, graphs) = graphs.into(2)

process LouvSoft {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('g.graphml'), truth from cl_in

    output:
    set key, file(cl_file), truth into ls_out

    script:
    cl_file = MetaSweeper.appendTo(key, "louvsoft.cl")

    if (params.debug) {
        """
        echo $key > $cl_file
        """
    }
    else {
        """
        louvain_cluster.py --otype soft --ofmt mcl g.graphml $cl_file
        """
    }
}


// copy channel
(cl_in, graphs) = graphs.into(2)

process LouvHard {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('g.graphml'), truth from cl_in

    output:
    set key, file(cl_file), truth into lh_out

    script:
    cl_file = MetaSweeper.appendTo(key, "louvhard.cl")

    if (params.debug) {
        """
        echo $key > $cl_file
        """
    }
    else {
        """
        louvain_cluster.py --otype hard --ofmt mcl g.graphml $cl_file
        """
    }
}


// copy channel
(cl_in, graphs) = graphs.into(2)

process Oclustr {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('g.graphml'), truth from cl_in

    output:
    set key, file(cl_file), truth into oc_out

    script:
    cl_file = MetaSweeper.appendTo(key, "oclustr.cl")

    if (params.debug) {
        """
        echo $key > $cl_file
        """
    }
    else {
        """
        oclustr.py -f mcl g.graphml $cl_file
        """
    }
}


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



// Collect all clustering results and score
bc_in = ls_out.mix(lh_out)

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

// fetch contig fastas from output dir
contigs = Channel.from(file('out/*.contigs.fasta'))
                 // add a string key for each iteration step
                 .map { f -> [MetaSweeper.dropSuffix(f.name), f] }

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
