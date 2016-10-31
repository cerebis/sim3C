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

MetaSweeper ms = MetaSweeper.fromFile(new File('sweep.yaml'))

// truth tables
truths = Channel.from(file("${ms.options.output}/*.truth"))
                .map { f -> [MetaSweeper.dropSuffix(f.name), f] }

// contig graphs
graphs = Channel.from(file("${ms.options.output}/*.graphml"))
                .map { f -> [MetaSweeper.dropSuffix(f.name), f] }

// copy channel
(cl_in, graphs) = graphs.into(2)

process LouvSoft {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('g.graphml') from cl_in

    output:
    set key, file(cl_file) into ls_out

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
    set key, file('g.graphml') from cl_in

    output:
    set key, file(cl_file) into lh_out

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
    set key, file('g.graphml') from cl_in

    output:
    set key, file(cl_file) into oc_out

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

//
//bc_sweep = truths.onCopy().cross(louvsoft_cl.map { f-> [helper.dropSuffix(f.name), f]}
//    .mix(louvhard_cl.map{ f-> [helper.dropSuffix(f.name), f]})
//    .mix(oclustr_cl.map{ f-> [helper.dropSuffix(f.name), f]})
//    .map { t -> [ helper.removeLevels(t[0],1), *t ] })
//    .map { t -> [ t[1][1], t[1][2], t[0][1]] }
//
//
//process Bcubed  {
//    cache 'deep'
//    publishDir params.output, mode: 'symlink', overwrite: 'true'
//
//    input:
//    set oname, file('clust'), file('truth') from bc_sweep
//
//    output:
//    set file("${oname}.bc") into cluster_bc
//
//    """
//    bcubed.py --tfmt json --ofmt json truth clust "${oname}.bc"
//    """
//}
//

contigs = Channel.from(file('out/*.contigs.fasta'))
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
        calc_N50_L50.py --ofmt yaml contigs.fa "${key}.asmstat"
        """
    }
}
