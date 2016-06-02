import groovyx.gpars.dataflow.DataflowQueue
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
        def elems = name.split(SEPARATOR)
        def a = elems[0..-(n+1)].join(SEPARATOR)
        def b = elems[-n..-1].join(SEPARATOR)
        return [a, b]
    }

    static String removeLevels(Path path, int n) {
        def name = path.toAbsolutePath().toString()
        return removeLevels(name, n)
    }

    static String removeLevels(String name, int n) {
        return name.split(SEPARATOR)[0..-(n+1)].join(SEPARATOR)
    }

}

// Runtime mixin
DataflowQueue.mixin DataflowUtils

helper = this.class.classLoader.parseClass(new File('Helper.groovy')).newInstance()
duplicator = this.class.classLoader.parseClass(new File('ChannelDuplicator.groovy')).newInstance()

def keyedFrom(path) {
    Channel.from(file(path)).map { f -> [helper.dropSuffix(f.name, SEPARATOR), f]}
}

g_stats = keyedFrom('out/*.graphml')
    .flatCross(keyedFrom('out/*.gstat'))
    .flatCross(keyedFrom('out/*.geigh'))

cl_stats = keyedFrom('out/*.cl')
    .flatCross(keyedFrom('out/*.bc'))
    .popLevels(1)

stat_sweep = g_stats.flatCross(cl_stats)

process Aggregate {

    input:
    set oname, file('graph'), file('gstat'), file('geigh'), method, file('cl'), file('bc') from stat_sweep
    
    output:
    set val(method), val(oname), stdout into all_stats
    
    """
    aggregate.py --fmt gstat geigh bc
    """
}

def slurper = new JsonSlurper()

def sweep_labels = ['ancestor','donor','alpha_BL','tree','profile','xfold','n3c']

all_stats = all_stats
    .map { t -> [ t[0], [sweep_labels, t[1].split(SEPARATOR)].transpose().collectEntries{ k, v -> [k, v] }, t[2] ] }
    .map { t -> t[1]['algo']=t[0]; t[1].putAll(slurper.parseText(t[2])); t[1]}
    .reduce([]) {a, b -> a.push(b); return a}


fout = file('all_stats.json')
fout << new JsonBuilder(all_stats.val).toString()
fout << '\n'
