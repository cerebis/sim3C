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

MetaSweeper ms = MetaSweeper.fromFile(new File('sweep.yaml'))

// fetch truth tables from sweep output directory
truths = MetaSweeper.keyedFrom(file("${ms.options.output}/*.truth"))

// fetch contig graphs from output dir and pair with respective truth table
graphs = MetaSweeper.keyedFrom(file("${ms.options.output}/*.graphml"))
                // reduce the key depth to that of the truth tables
                .map { t -> [t[0].popLevels(1), *t] }
                // join channels, using new reduced key as the index
                .phase(truths)
                // unwrap lists and remove redundant elements
                .map { t -> [*t[0][1..-1], t[1][1]] }


// prepare input channel
(cl_in, graphs) = graphs.into(2)
cl_in = cl_in.map { t -> t.extendKey('algo', 'louvsoft') }

process LouvSoft {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('g.graphml'), truth from cl_in

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

ls_out = ls_out.map { it.nameify(1, 'ls_clustering') }

// prepare input channel
(cl_in, graphs) = graphs.into(2)
cl_in = cl_in.map { t -> t.extendKey('algo', 'louvhard') }

process LouvHard {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('g.graphml'), truth from cl_in

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

lh_out = lh_out.map { it.nameify(1, 'lh_clustering') }

// prepare input channel
(cl_in, graphs) = graphs.into(2)
cl_in = cl_in.map { t -> t.extendKey('algo', 'oclustr') }

process Oclustr {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('g.graphml'), truth from cl_in

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

oc_out = oc_out.map { it.nameify(1, 'oc_clustering') }

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

ginfo_out = ginfo_out.map { it.nameify(1, 'gstat'); it.nameify(2, 'geigh') }


// Collect all clustering results and score
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

bc_out = bc_out.map { it.nameify(1, 'bc_scores') }

// fetch contig fastas from output dir
contigs = MetaSweeper.keyedFrom(file("${ms.options.output}/*.contigs.fasta"), 2)

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

astat_out = astat_out.map { it.nameify(1, 'asmstat') }
