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
import groovy.transform.AutoClone
import groovy.transform.Synchronized
import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Nextflow
import nextflow.Channel
import java.nio.file.Path

/**
 * Utility methods used in the construction of meta-sweep workflows.
 */
class Helper {
    static String SEPARATOR = '-+-'
    static def delimiters = /[ ,\t]/

    static {
        List.metaClass.pick = { ... picks ->
            if (picks instanceof Object[]) {
                picks = picks.collect()
            }
            else if (! (picks instanceof Collection)) {
                picks = [picks]
            }
            else {
                picks = picks.unique()
            }
            picks.inject([] as Set){acc, pi -> pi
                owner.delegate.each {
                    if (!(it instanceof NamedValue) || it.name == pi) {
                        acc << it
                    }
                }
                acc
            }
        }

        List.metaClass.addKey = {
            Key key = new Key()
            delegate.each { key.put(it) }
            [key, *delegate]
        }

        List.metaClass.nameify = { ix, name ->
            def val = delegate[ix]
            delegate[ix] = new NamedValue(name, val)
            return delegate
        }
    }


    @AutoClone
    static class NamedValue {
        private static final int HASH_PRIME = 7
        String name
        Object value

        NamedValue(name, value) {
            this.name = name
            this.value = value
        }

        String toString() {
            if (value instanceof Collection) {
                value.toSorted().join(' ')
            }
            else {
                value.toString()
            }
        }

        public static NamedValue create(name, value) {
            new NamedValue(name,value)
        }

        boolean equals(Object obj) {
            if (!(obj instanceof NamedValue)) {
                return false
            }
            if (this.is(obj)) {
                return true
            }
            NamedValue other = obj as NamedValue
            return other.name == name && other.value == value
        }

        int hashCode() {
            return HASH_PRIME * name.hashCode() * value.hashCode()
        }
    }

    static class Key {
        private static final int HASH_PRIME = 37
        LinkedHashMap varMap = [:]

        @Synchronized
        Key copy() {
            Key newKey = new Key()
            varMap.each{ k, v -> newKey.put(v.clone())}
            return newKey
        }

        @Synchronized
        Map putAll(Collection v) {
            def ret = [:]
            v.each{ ret[it.name] = put(it) }
            return ret
        }

        @Synchronized
        Object put(NamedValue v) {
            def old = varMap.putIfAbsent(v.name, v)
            assert old == null: "attempted duplicate key insertion \"$v.name\" with existing keys: ${getKeys()}"
            return old
        }

        @Synchronized
        Object put(String k, Object v) {
            assert !(v instanceof NamedValue): 'NamedValue objects should be passed alone'
            return put(new NamedValue(k, v))
        }

        @Synchronized
        Object getAt(String k) {
            return varMap[k]
        }

        String toString() {
            return id()
        }

        String id(Collection vars) {
            return vars.collect{ simpleSafeString(it.value) }.join(SEPARATOR)
        }

        String id(int maxDepth) {
            def values = varMap.values() as List
            assert maxDepth > 0 : 'Depth must be a positive integer'
            assert maxDepth <= values.size() : "Requested depth [$maxDepth] exceeds defined variable count [${values.size()}]"
            return id(values[0..maxDepth-1])
        }

        String id() {
            return id(varMap.values())
        }

        boolean equals(Object other) {
            if (!(other instanceof Key)) {
                return false
            }
            if (this.is(other)) {
                return true
            }
            Key _o = other as Key
            return _o.varMap.keySet().equals(varMap.keySet())
        }

        int hashcode() {
            return HASH_PRIME * varMap.keySet().hashcode()
        }

        public def getKeys() {
            varMap.keySet() as List
        }

        List orderedJoin(Key other) {
            List a = getKeys() as List
            List b = other.getKeys() as List
            List ret = []
            for (int i = 0; i < a.size() && i < b.size(); i++) {
                if (a[i] != b[i]) {
                    break
                }
                ret[i] = a[i]
            }
            return ret
        }

        List join(Key other) {
            getKeys().intersect(other.getKeys()) as List
        }

    }


    static class Sweep {
        @Delegate
        Map varRegistry = [:]

        public Object put(Object key, Object value) {
            varRegistry[key] = value.collect { new NamedValue(key, it) }
        }

        public String description() {
            int width = varRegistry.keySet().inject(0) {
                acc, it -> acc = it.length() > acc ? it.length() : acc
            }

            StringBuffer desc = new StringBuffer()
            def spc = '\nName '.padRight(width+1)
            desc.append("Variations defined for this sweep\n${spc} Values\n")
            varRegistry.collect{ k, v ->
                def label = "${k}:"
                desc.append("${label.padRight(width+1)} $v\n")
            }
            return desc.toString()
        }

        public Collection permuteAll() {
            permute(values())
        }

        public Collection permuteNames(Object... varNames) {
            permute(subMap(varNames.collect()).values())
        }

        public Collection permuteLevels(int begin = 0, int end = -1) {
            def names = keySet() as List
            permute(subMap(names[begin..end]).values())
        }

        protected Collection permute(Collection collections) {
            if (collections == null || collections.isEmpty()) {
                return []
            } else {
                List res = []
                permuteImpl(collections, res, 0, [])
                return res
            }
        }

        /** Recursive implementation for {@link #permutations(List, Collection)} */
        private static void permuteImpl(ori, res, d, current) {
            // if depth equals number of original collections, final reached, add and return
            if (d == ori.size()) {
                res.add(current)
                return
            }

            // iterate from current collection and copy 'current' element N times, one for each element
            def currentCollection = ori[d];
            for (element in currentCollection) {
                List copy = current.collect()
                copy.add(element);
                permuteImpl(ori, res, d + 1, copy)
            }
        }

        public DataflowQueue permutedChannel(Object... varNames) {
            Channel.from(permuteNames(varNames))
                    .map{ it.addKey() }
        }

        public DataflowQueue extendChannel(DataflowQueue df, Object... varNames) {
            assert varNames.size() > 0 : "Error, no variables named in extendChannel"
            def p = permuteNames(varNames)
            if (p.size() == 0) {
                return df
            }

            df = df.spread(p)
                    .map { row ->
                        def key = row[0].copy();
                        def vars = row[-varNames.size()..-1]
                        vars.each { key.put(it) }

                        proc_vars = []
                        int num_proc_vars = row.size() - varNames.size() - 1
                        if (num_proc_vars > 0) {
                            row[1..num_proc_vars].each { proc_vars << it }
                        }

                        [key, *proc_vars]
                    }

            return df
        }

        public DataflowQueue joinChannels(DataflowQueue df1, DataflowQueue df2, Integer... level) {
            assert level != null &&
                    level.size() >= 1 && level.size() <= 2 : 'Level can be either 1 or 2 integer values'

            if (level.size() == 1) {
                level = [level[0], level[0]]
            }

            df1 = df1.map { it ->
                // remake key
                [ it[0].id(level[0]), *it[1..-1] ] }

            df2 = df2.map { it ->
                // remake key
                [ it[0].id(level[1]), *it[1..-1] ] }

            return df1.cross(df2).map { it.flatten() }
        }

    }

    static def selectByName(keys,list) {
        def m = list.inject([:]) {acc, it -> acc[it.name] = it; acc}
        keys.inject([]) { acc, k -> assert m[k] != null, "element [$k] not found in map"; acc << m[k]; acc }
    }

    /**
     * Add a name to every object in a list.
     */
    static def addName(name, list) {
        list.collect{new NamedValue(name, it)}
    }

    /**
     * Return absolute path from string.
     *
     * @param path - path as string, wildcards allowed
     * @return java.nio.file.Path or list of
     */
    static def absPath(path) {
        Nextflow.file(path)*.toAbsolutePath()
    }

    /**
     * Create string representation of value. If this is an instance of {@see java.nio.file.Path},
     * then use only the file name. In the context of filesystem restrictions, mask problematic
     * characters if safe=true.
     *
     * @param val - object to represent as a string
     * @param safe - a flag that when true masks problematic characters
     * @return String
     */
    static def simpleSafeString(val, safe=true) {
        if (val instanceof java.nio.file.Path) {
            val = val.name - ~/[\.][^\.]+$/
        }
        else {
            val = val
        }
        return safe ? safeString(val) : val.toString()
    }

    /**
     * When performing combinatorial sweeps, channels represent tables of varied parameters,
     * where each row corresponds to an individual combination of variable values.
     * <p>
     * String representations of a value-set (one value for each variable involved in sweep stage) are
     * joined together to act as an proxy identifier or key. These keys are used both in naming output
     * files of a process and as a means of joining channels together along relevant combinations.
     * <p>
     * All columns in a row are returned, prepended by the new key.
     *
     * @param row - channel row, where columns between firstCol and lastCol contain the relevant values from
     * which to generate a key.
     * @param firstCol - the first column in the row to use in the key.
     * @param lastCol - the last column in the row to use in the key.
     * @return [key, *row]
     */
    static def addKey(row, firstCol=0, lastCol=-1) {
        [row[firstCol..lastCol].collect { simpleSafeString(it.value) }.join(SEPARATOR), *row]
    }

    /**
     * Select a subset of elements from a collection.
     * <p>
     * This is mostly for readability and is used to reduce the outcome of channel joins
     * to only the desired columns. Where often, joins contain redundant or superfluous
     * columns.
     *
     * @param row - the collection to slice
     * @param elems - the element indices to return
     * @return the sliced collection
     */
    static def selectById(row, elems) {
        row[elems]
    }

    /**
     * Return a safer representation of a string. Where "safe" is in the context
     * of filesystem restrictions. Therefore, masking of backslash characters and
     * the use only of the base file name when provided with an instance of
     * {@see java.nio.file.Path}
     *
     * @param val - object to return as a stirng
     * @return String representation
     */
    static String safeString(val) {
        def s
        if (val instanceof Path) {
            s = val.name
        }
        else {
            s = val.toString()
        }
        return s.replaceAll(/[\\\/]/, "_")
    }

    static int[] stringToInts(String str) {
        return (str.split(delimiters) - '').collect { elem -> elem as int }
    }

    static float[] stringToFloats(String str) {
        return (str.split(delimiters) - '').collect { elem -> elem as float }
    }

    /**
     * Convert a whitespace or comma delimited String to a List.
     * @param str - the string to split
     * @return List
     */
    static String[] stringToList(String str) {
        return str.split(delimiters) - ''
    }

    static String dropSuffix(str, sep) {
        def p1 = str.lastIndexOf(sep)
        if (p1 == -1) {
            return str
        }
        return str.indexOf('.', p1).with { it != -1 ? str[0..<it] : str }
    }

    static String dropSuffix(str) {
        return str.lastIndexOf('.').with {it != -1 ? str[0..<it] : str}
    }

    static String[] splitSampleName(Path path) {
        def m = (path.name =~ /^(.*)_?([rR][0-9]+).*$/)
        return m[0][1..2]
    }

    static String removeLevels(Path path, int n) {
        def name = path.toAbsolutePath().toString()
        return name.split(SEPARATOR)[0..-(n+1)].join(SEPARATOR)
    }

    static String removeLevels(String name, int n) {
        return name.split(SEPARATOR)[0..-(n+1)].join(SEPARATOR)
    }

    static Object[] product(A, B) {
        return A.collectMany{a->B.collect{b->[a, b]}}
    }
}

