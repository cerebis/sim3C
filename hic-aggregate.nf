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

//MetaSweeper ms = MetaSweeper.fromFile(new File('sweep.yaml'))


import groovyx.gpars.dataflow.DataflowQueue
import java.util.regex.Pattern
import groovy.json.JsonSlurper
import groovy.json.JsonBuilder
import static Globals.*

class Globals {
    static final String SEPARATOR = '%'
}

class DataflowUtils {

    /**
     * Cross this channel with another, but flatten the resulting rows
     * and keep only one key element in the first position.
     */
    DataflowQueue flatCross(DataflowQueue b) {
        return this.cross(b).map{ t -> [*t[0], t[1][1..-1]].flatten() }
    }

    /**
     * Remove specified number of sweep levels from key.
     */
    DataflowQueue removeLevels(int levels) {
        return this.map { t -> [ DataflowUtils.removeLevels(t[0], levels), t[1..-1] ].flatten() }
    }

    /**
     * Remove specified number of sweep levels from key, but retain them in the row.
     */
    DataflowQueue popLevels(int levels) {
        return this.map { t -> [ DataflowUtils.popLevels(t[0], levels), t[1..-1] ].flatten() }
    }

    static String[] popLevels(Path path, int n) {
        def name = path.toAbsolutePath().toString()
        return popLevels(name, n)
    }

    static String[] popLevels(String name, int n) {
        def elems = name.split(Pattern.quote(MetaSweeper.SEPARATOR))
        def a = elems[0..-(n+1)].join(MetaSweeper.SEPARATOR)
        def b = elems[-n..-1].join(MetaSweeper.SEPARATOR)
        return [a, b]
    }

    static String removeLevels(Path path, int n) {
        def name = path.toAbsolutePath().toString()
        return removeLevels(name, n)
    }

    static String removeLevels(String name, int n) {
        return name.split(Pattern.quote(MetaSweeper.SEPARATOR))[0..-(n+1)].join(MetaSweeper.SEPARATOR)
    }

}

// Runtime mixin
DataflowQueue.mixin DataflowUtils


// collect grpahs and their statistics
g_stats = MetaSweeper.keyedFrom(file('out/*.graphml'))
    .flatCross(MetaSweeper.keyedFrom(file('out/*.gstat')))
    .flatCross(MetaSweeper.keyedFrom(file('out/*.geigh')))


// collect clusterings and their statistics
cl_stats = MetaSweeper.keyedFrom(file('out/*.cl'))
    .flatCross(MetaSweeper.keyedFrom(file('out/*.bc')))
    // split off the algorithm level so that cl_stats key is equivalent to g_stats key
    .popLevels(1)

// join these two channels on their shared key
stat_sweep = g_stats.flatCross(cl_stats)

process Aggregate {

    input:
    set key, file('graph'), file('gstat'), file('geigh'), algo, file('cl'), file('bc') from stat_sweep
    
    output:
    set val(algo), val(key), stdout into all_stats
    
    """
    aggregate.py --fmt yaml gstat geigh bc
    """
}

/*
def slurper = new JsonSlurper()

def sweep_labels = ['ancestor','donor','alpha_BL','tree','profile','xfold','n3c']

all_stats = all_stats
    .map { t -> [ t[0], [sweep_labels, t[1].split(SEPARATOR)].transpose().collectEntries{ k, v -> [k, v] }, t[2] ] }
    .map { t -> t[1]['algo']=t[0]; t[1].putAll(slurper.parseText(t[2])); t[1]}
    .reduce([]) {a, b -> a.push(b); return a}


fout = file('all_stats.json')
fout << new JsonBuilder(all_stats.val).toString()
fout << '\n'
*/