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
import static Globals.*
class Globals {
    static final String SEPARATOR = '%'
}
helper = this.class.classLoader.parseClass(new File('Helper.groovy')).newInstance()
duplicator = this.class.classLoader.parseClass(new File('ChannelDuplicator.groovy')).newInstance()

// truth tables
truths = Channel.from(file('out/*.truth'))
        .map { f -> [f.name[0..-7], f] }
truths = duplicator.createFrom(truths)

// contig graphs
graphs = Channel.from(file('out/*.graphml'))
        .map { f -> [helper.dropSuffix(f.name, '%'), f] }
graphs = duplicator.createFrom(graphs)

gr_sweep = graphs.onCopy()

process LouvSoft {
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'true'

    input:
    set oname, file('g.graphml') from gr_sweep

    output:
    set file("${oname}${SEPARATOR}louv-soft.cl") into louvsoft_cl

    """
    louvain_cluster.py --otype soft --ofmt mcl g.graphml "${oname}.louv-soft.cl"
    """
}

gr_sweep = graphs.onCopy()

process LouvHard {
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'true'

    input:
    set oname, file('g.graphml') from gr_sweep

    output:
    set file("${oname}${SEPARATOR}louv-hard.cl") into louvhard_cl

    """
    louvain_cluster.py --otype hard --ofmt mcl g.graphml "${oname}.louv-hard.cl"
    """
}

gr_sweep = graphs.onCopy()

process Oclustr {
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'true'

    input:
    set oname, file('g.graphml') from gr_sweep

    output:
    set file("${oname}${SEPARATOR}oclustr.cl") into oclustr_cl

    """
    oclustr.py -f mcl g.graphml "${oname}.oclustr.cl"
    """
}

gr_sweep = graphs.onCopy()

process GraphStats {
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'true'

    input:
    set oname, file('g.graphml') from gr_sweep

    output:
    set file("${oname}.gstat"), file("${oname}.geigh") into graph_info

    """
    graph_stats.py --ofmt json g.graphml "${oname}.gstat"
    graph_complexity.py --ofmt json --method eigh g.graphml "${oname}.geigh"
    """
}

louvsoft_cl.map { f-> [helper.dropSuffix(f.name), f]}
    .mix(louvhard_cl.map{ f-> [helper.dropSuffix(f.name), f]})
    .mix(oclustr_cl.map{ f-> [helper.dropSuffix(f.name), f]})
    //.subscribe{println it}

bc_sweep = truths.onCopy().cross(louvsoft_cl.map { f-> [helper.dropSuffix(f.name), f]}
    .mix(louvhard_cl.map{ f-> [helper.dropSuffix(f.name), f]})
    .mix(oclustr_cl.map{ f-> [helper.dropSuffix(f.name), f]})
    .map { t -> [ helper.removeLevels(t[0],1), *t ] })
    .map { t -> [ t[1][1], t[1][2], t[0][1]] }

process Bcubed  {
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'true'

    input:
    set oname, file('clust'), file('truth') from bc_sweep

    output:
    set file("${oname}.bc") into cluster_bc

    """
    bcubed.py --tfmt json --ofmt json truth clust "${oname}.bc"
    """
}

contigs = duplicator.createFrom(Channel.from(file('out/*.contigs.fasta'))
    .map { f -> [f.name[0..-15], f] })

contig_sweep = contigs.onCopy()

process AssemblyStats {
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'true'

    input:
    set oname, file('contigs.fa') from contig_sweep

    output:
    set file("${oname}.asmstat") into asm_stats

    """
    calc_N50_L50.py --ofmt json contigs.fa "${oname}.asmstat"
    """
}

/*
cluster_bc.map { f -> [helper.dropSuffix(f.name), f] }
    .subscribe{println it}
*/