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
package mzd

@Grab('org.yaml:snakeyaml:1.17')
@Grab('com.google.guava:guava:19.0')
import groovy.transform.AutoClone
import groovy.transform.Synchronized
import groovyx.gpars.dataflow.DataflowChannel
import groovyx.gpars.dataflow.DataflowQueue
import java.nio.file.Path
import java.io.File
import java.io.Writer
import java.nio.file.Files
import nextflow.Channel
import nextflow.Nextflow
import nextflow.util.KryoHelper
import org.yaml.snakeyaml.Yaml
import org.yaml.snakeyaml.TypeDescription
import org.yaml.snakeyaml.constructor.Constructor
import com.google.common.hash.Hashing
import com.google.common.hash.Hasher
import com.google.common.hash.HashFunction
import com.google.common.hash.Funnels

/**
 * Utility methods used in the construction of meta-sweep workflows.
 */
class MetaSweeper {
    public Map<String, Object> variables
    public Map<String, Object> options
    public sweep

    static String SEPARATOR = '-+-'
    static def delimiters = /[ ,\t]/

    static {

        // For Nextflow caching to function, we must register our custom classes
        // with the chosen serialization library (Kryo). We do not provide any
        // custom Serializer implementation.
        KryoHelper.register(Key);
        KryoHelper.register(NamedValue);

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

    static MetaSweeper fromFile(File config) {
        ConfigLoader.read(config)
    }

    public String describeSweep() {
        ConfigLoader.write(this)
    }

    public MetaSweeper createSweep() {
        sweep = new Sweep()
        return this
    }

    public DataflowChannel extend(DataflowQueue df, Object... varNames) {
        sweep.extendChannel(df, varNames)
    }

    public MetaSweeper withVariable(String name, Boolean stepInto=false) {
        if (stepInto) {
            sweep[name] = variables[name].collect()
        }
        else {
            sweep[name] = variables[name]
        }
        return this
    }

    public DataflowQueue permute() {
        sweep.permutedChannel()
    }


    @AutoClone
    static class NamedValue {
        String name
        Object value

        /**
         * A default constructor is required for Kryo deserialization.
         */
        NamedValue() {/*...*/}

        /**
         * Basic constructor.
         * @param name
         * @param value
         */
        NamedValue(name, value) {
            this.name = name
            this.value = value
        }

        /**
         * String representations are of the underlying values only, collections are
         * concatenated by whitespace.
         * @return String
         */
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

        /**
         * NamedValue identity is determined by both name and value.
         * @param other
         * @return
         */
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

        /**
         * NamedValue identity is determined by both name and value.
         * @return
         */
        int hashCode() {
            return Objects.hash(name, value)
        }
    }

    static class Key {
        LinkedHashMap varMap = [:]

        /**
         * A default constructor is required for deserialization.
         *
         * NOTE. Though not required when no other constructors have been implemented, it has been added
         * to avoid the possibility of introducing an unintended side-effect -- if one was implemented
         * later. This would then break Kryo deserialization, which expects to be able toinstantiate
         * objects without parameters.
         */
        Key() {/*...*/}

        /**
         * Analogous to String concatenation, the + operator is not reflexive. (a+b != b+a)
         *
         * @param other
         * @return Key
         */
        @Synchronized
        Key plus(Key other) {
            Key cpy = new Key()
            cpy.varMap.putAll(this.varMap)
            cpy.varMap.putAll(other.varMap)
            return cpy
        }

        @Synchronized
        List values() {
            varMap.collect{ k, v -> v}
        }

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
            List values = varMap.values() as List
            assert maxDepth > 0 : 'Depth must be a positive integer'
            assert maxDepth <= values.size() : "Requested depth [$maxDepth] exceeds defined variable count [${values.size()}]"
            return id(values[0..maxDepth-1])
        }

        String id() {
            return id(varMap.values())
        }

        /**
         * Key identity considers only the keys in the contained Map (varMap)
         * @param other
         * @return
         */
        boolean equals(Object other) {
            if (!(other instanceof Key)) {
                return false
            }
            if (this.is(other)) {
                return true
            }
            Key _o = other as Key

            return _o.toString().equals(toString())
        }

        /**
         * Key identity considers only the keys in the contained Map (varMap)
         * @return
         */
        int hashCode() {
            return Objects.hash(toString())
        }

        public def getKeys() {
            varMap.keySet() as List
        }

        List join(Key other) {
            getKeys().intersect(other.getKeys()) as List
        }

        /**
         * Return a subkey, where the depth is less than the maximum depth of the sweep.
         * @param depth - the depth to retain
         * @return reduced key
         */
        public Key subKey(int depth) {
            List values = varMap.values() as List
            assert depth > 0 : 'Depth must be a positive integer'
            assert depth <= values.size() : "Requested depth [$depth] exceeds defined variable count [${values.size()}]"
            Key k = new Key()
            k.putAll(values[0..(depth-1)])
            return k
        }

        /**
         * At a given sweep level, split a key into two parts.
         *
         * Neither part can be empty, that is you cannot split at the lowest or highest level.
         *
         * @param level - the sweep level at which to split (0 > level < max_level)
         * @return Map containing two keys: 'hi' and 'lo'
         */
        public Map splitKey(int level) {
            def values = varMap.values() as List
            assert level > 0 : 'Depth must be > 0'
            assert level < values.size() : "Requested cut point [$level] would produce an empty split for max [${values.size()}] levels"
            def keys = ['lo':new Key(), 'hi':new Key()]
            keys.lo.putAll(values[0..level-1])
            keys.hi.putAll(values[level..-1])
            return keys
        }
    }

    /**
     * Community represents a collection of clades/groups of related genomes.
     * Together it can be regarded as a simulated community.
     */
    static class Community implements Iterable {
        // just a human centric name for a community
        public String name
        // the list of clades/groups to generate
        public List<Clade> clades

        public Iterator iterator() {
            clades.iterator()
        }

        /**
         * Hash file content from file name
         * @param hasher
         * @param fileName
         * @return
         */
        static private Hasher hashFileContent(Hasher hasher, String fileName) {
            hashFileContent(hasher, new File(fileName))
        }

        /**
         * Hash file content. Appropriated from nextflow.util.CacheHelper as it
         * has been defined private. (removed further nextflow specific helper
         * classes)
         *
         * @param hasher
         * @param file
         * @return
         */
        static private Hasher hashFileContent(Hasher hasher, File file) {
            OutputStream output = Funnels.asOutputStream(hasher);
            try {
                Files.copy(file.toPath(), output);
            }
            catch( IOException e ) {
                throw new IllegalStateException("Unable to hash content: ${file}", e);
            }
            finally {
                if (output) {
                    output.close()
                }
            }
            return hasher;
        }

        /**
         * Create a hash for a given community definition. This is used
         * as a shorthand identifier when naming files. Everything involved
         * in the definition of a community is used except the name -- which would
         * play no role in structure.
         *
         * Both ancestor and donor sequences are deeply hashed.
         *
         * Only mumur3_32 is used since we do not desire an enormously long
         * string and we doubt that there would be clashes for any realistic
         * range of communities in a given experiment.
         *
         * @return String representation of community hash
         */
        public String hash() {

            HashFunction hf = Hashing.murmur3_32()
            Hasher hshr = hf.newHasher()

            // hash clades definitions
            clades.each { cl ->
                // deeply hash ancestor and donor
                hashFileContent(hshr, cl.ancestor)
                hashFileContent(hshr, cl.donor)

                hshr.putInt(cl.ntaxa)

                // include the definition of tree and profile
                cl.tree.collect { e ->
                    hshr.putUnencodedChars(e.getKey())
                            .putFloat(e.getValue())
                }
                cl.profile.collect { e ->
                    hshr.putUnencodedChars(e.getKey())
                            .putFloat(e.getValue())
                }
            }
            hshr.hash().toString()
        }

        /**
         * Represented by the name and hash
         * @return String
         */
        @Override
        String toString() {
            "${name}-${hash()}".toString()
        }

    }

    /**
     * A Clade represents set of related genomes evolved from
     * a common ancestor. How these sequences are different
     * is dictated by a phylogenetic tree and how they would
     * manifest in a ecosystem by an abundance profile.
     *
     * Both trees and profiles are expected to be generated
     * through simulation -- which itself is defined by a set
     * of parameters. These simulation parameters are stored
     * as a map for either property and intended simulation
     * tool/algorithm.
     *
     * How many taxa in exist in the clade is defined by ntaxa.
     */
    static class Clade {
        // common ancestor
        private Path ancestor
        // donor used for htg
        private Path donor
        // number of leaves/tips.
        public Integer ntaxa
        // parameters involved in generating the tree
        public Map<String, Float> tree
        // parameters involved in generating the profile
        public Map<String, Float> profile

        public void setAncestor(String ancestor) {
            this.ancestor = Nextflow.file(ancestor)
        }
        public String getAncestor() {
            return ancestor.toString()
        }
        public Path getAncestorPath() {
            return ancestor
        }

        public void setDonor(String donor) {
            this.donor = Nextflow.file(donor)
        }
        public String getDonor() {
            return donor.toString()
        }
        public Path getDonorPath() {
            return donor
        }

        /**
         * String representation of the Clade as a concatenated set of its parameters
         * @return String
         */
        @Override
        public String toString() {
            def l = [simpleSafeString(ancestor), simpleSafeString(donor), ntaxa] + mapString(tree) + mapString(profile)
            l.flatten().join(SEPARATOR)
        }

        protected static List mapString(Map<?,?> m) {
            def sorted_keys = m.sort()*.key
            sorted_keys.collect{ k -> "${m[k]}" }
        }
    }

    /**
     * Read sweep configurations in YAML format
     *
     * Should be threadsafe.
     */
    static class ConfigLoader {

        private static ThreadLocal<Yaml> threadLocal

        static {
            /**
             * ThreadLocal instance of Yaml read/write object
             */
            threadLocal = new ThreadLocal<Yaml>() {
                @Override
                protected Yaml initialValue() {
                    // constructor root is MetaSweeper
                    Constructor cnstr = new Constructor(MetaSweeper.class)
                    // two tags are used to instantiate custom classes from
                    // within the map hierarchy
                    cnstr.addTypeDescription(new TypeDescription(Community.class, '!com'))
                    cnstr.addTypeDescription(new TypeDescription(Clade.class, '!clade'))
                    cnstr.addTypeDescription(new TypeDescription(File.class))
                    new Yaml(cnstr)
                }
            }
        }

        /**
         * Read a configuration from a file.
         * @param file -- YAML file to read
         * @return MetaSweeper
         */
        public static MetaSweeper read(File file) {
            Yaml yaml = threadLocal.get()
            yaml.load(new FileReader(file))
        }

        /**
         * Read a configuration from a string
         * @param doc -- String to read
         * @return MetaSweeper
         */
        public static Map read(String doc) {
            Yaml yaml = threadLocal.get()
            yaml.load(doc)
        }

        /**
         * Convert MetaSweeper object to a string in YAML syntax
         * @param ms -- Object to convert
         * @return String (YAML)
         */
        public static String write(MetaSweeper ms) {
            Yaml yaml = threadLocal.get()
            yaml.dump(ms)
        }

        /**
         * Write a YAML representation of a MetaSweeper object to a stream
         * @param ms - Object to write
         * @param writer - output writer
         */
        public static void write(MetaSweeper ms, Writer writer) {
            Yaml yaml = threadLocal.get()
            yaml.dump(ms, writer)
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
            def ln = varNames.collect()
            ln = subMap(ln)
            permute(ln.values())
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

        /** Recursive implementation for permutations(List, Collection) */
        private static void permuteImpl(ori, res, d, current) {
            // if depth equals number of original collections, final reached, add and return
            if (d == ori.size()) {
                res.add(current)
                return
            }

            // iterate from current collection and copy 'current' element N times, one for each element
            def currentCollection = ori[d];
            for (element in currentCollection) {
                List cpy = current.collect()
                cpy.add(element);
                permuteImpl(ori, res, d + 1, cpy)
            }
        }

        public  DataflowQueue permutedChannel() {
            List allVars = varRegistry.keySet() as List
            permutedChannel(*allVars)
        }

        public DataflowQueue permutedChannel(Object... varNames) {
            def pn = permuteNames(varNames)
            Channel.from(pn).map{ it ->
                try {
                    it.addKey()
                }
                catch (Throwable t) {
                    println "ITER $it"
                    println "EXCEPTION: $t"
                    throw ex
                }
            }
        }

        public DataflowQueue extendChannel(DataflowQueue df, Object... varNames) {
            assert varNames.size() > 0 : "Error, no variables named in extendChannel"
            def p = permuteNames(varNames)
            if (p.size() == 0) {
                return df
            }

            df = df.spread(p)
                    .map { row ->
                        // create a copy of existing key
                        def key = row[0].copy();
                        // get the set of new variables added to each row
                        def new_vars = row[-varNames.size()..-1]
                        // add each new variable to the key
                        new_vars.each { key.put(it) }

                        def all_vars = row[1..-1].collect() // { acc, it -> acc << it }
                        //println "ALL $all_vars"
                        [key, *all_vars]
                    }

            return df
        }

        /**
         * Group the elements of a channel by a common key at a specific level in the
         * sweep hierarchy.
         *
         * It is expected that the specified level would be such that multiple rows in
         * the table will share a common key value. i.e 1 < level < max_level
         *
         * The returned map is indexed by the common key. Its elements are lists which
         * hold the "overhang" (remaining key levels dropped for this grouping) and the
         * row values for each grouped set of values.
         *
         * @param df - channel to group
         * @param level - sweep level to group at
         * @return map of grouped rows.
         */
        static Map groupAtLevel(DataflowQueue df, int level) {
            List list = df.toList().get()
            return list.inject([:]) { acc, row ->
                def k = row[0].splitKey(level)
                if (!acc[k.lo]) {
                    acc[k.lo] = []
                }
                acc[k.lo].add( ['overhang': k.hi, 'values':row[1..-1]] )
                acc
            }
        }

        /**
         * Join two channels together at a given level in the sweep hierarchy. This is used
         * when merging the results of two branches in the sweep.
         *
         * Such cases arise when inputs from a single point are used by two different tasks (branch point)
         * and then their output channels merged together for further processing (merge point).
         *
         * The merged channel elements are lists (table rows) composed of a leading key (which reflects the join point)
         * and followed by all the per-task output NamedValueS which were emitted by both channels. As such, it can be
         * a good idea to reduce rows to only the output values which will be required in subsequent tasks. This can be
         * accomplished using helper meta-method "List.pick(varname1, varname2...).
         *
         * Eg. If two output columns were 'reads' and 'contigs', using the map operator in Nextflow
         *
         *      chan.map{ it.pick('reads', 'contigs') }
         *
         * @param df1 - channel 1
         * @param df2 - channel 2
         * @param level - the level of the sweep hierarchy to merge at.
         * @return channel of merged values
         */
        public DataflowQueue joinChannels(DataflowQueue df1, DataflowQueue df2, int level) {
            Map g1 = groupAtLevel(df1, level)
            Map g2 = groupAtLevel(df2, level)

            // combine values from each input channel for common keys
            List ret = g1.keySet().inject([]) { acc, k ->
                if (g2[k]) {
                    def x = [g1[k], g2[k]].combinations()
                    x = [[k], x].combinations{ a, b ->
                        [a, b['overhang'].flatten(), b['values'].flatten()]
                    }
                    acc.addAll(x)
                }
                acc
            }

            // join keys in sweep order and return rows as [key, output_values]
            ret = ret.collect { row ->
                Key joined = orderKey(row[0] + row[1][0] + row[1][1])
                [joined, *row[2]]
            }

            return Channel.from(ret)
        }

        /**
         * Impose the sweep variable order on a key
         * @param key - the key to reorder
         * @return the reordered key
         */
        public Key orderKey(Key key) {
            return varRegistry.inject(new Key()) { ordKey, name, values ->
                if (key[name]) {
                    ordKey.put(key[name])
                }
                ordKey
            }
        }

    }

    static def selectByName(keys, list) {
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
        else if (val instanceof java.io.File) {
            val = val.toPath()
            val = val.name - ~/[\.][^\.]+$/
        }
        else { val = val }
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
