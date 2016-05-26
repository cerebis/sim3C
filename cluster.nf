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

class Globals {
    static String separator = '%'
}
helper = this.class.classLoader.parseClass(new File('Helper.groovy')).newInstance()


graphs = Channel.from(file('out/*graphml'))
        .map { f -> [f, helper.dropSuffix(f.name)] }
        .subscribe { println it }

process LouvSoft {
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'true'

    input:
    set file('g.graphml'), oname from graphs

    output:
    set file("${oname}.louv-soft.cl") into clusterings

    """
    louvain_cluster.py --otype soft --ofmt mcl g.graphml "${oname}.louv-soft.cl"
    """
}

process LouvHard {
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'true'

    input:
    set file('g.graphml'), oname from graphs

    output:
    set file("${oname}.louv-hard.cl") into clusterings

    """
    louvain_cluster.py --otype hard --ofmt mcl g.graphml "${oname}.louv-hard.cl"
    """
}

process Oclustr {
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'true'

    input:
    set file('g.graphml'), oname from graphs

    output:
    set file("${oname}.louv-hard.cl") into clusterings

    """
    oclustr.py -f mcl g.graphml "${oname}.oclustr.cl"
    """
}
