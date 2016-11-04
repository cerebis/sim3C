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
import java.util.regex.Pattern
import nextflow.Channel
import nextflow.Nextflow
import nextflow.util.KryoHelper
import org.yaml.snakeyaml.Yaml
import org.yaml.snakeyaml.DumperOptions
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
    public Map<String, Object> variables = [:]
    public Map<String, Object> options = [:]
    public sweep

    // NOTE: File system dependent! Both separator patterns must be composed of
    // legal, non-escaped characters of the underlying filesystem. Escaped characters
    // appear to be handled badly when wildcards are used within Nextflow file()
    static String KEYVALUE_SEP = '#'
    static String PARAM_SEP = '-+-'

    static def charsToReplace = /[\\\/ \t;!?*"']/
    static def delimiters = /[ ,\t]/
    static int MAX_LINE_LENGTH = 1024

    static {

        // Mixin some convenience methods.
        DataflowQueue.mixin DataflowUtils

        // For Nextflow caching to function, we must register our custom classes
        // with the chosen serialization library (Kryo). We do not provide any
        // custom Serializer implementation.
        KryoHelper.register(Key)
        KryoHelper.register(NamedValue)
        KryoHelper.register(Community)
        KryoHelper.register(Clade)

        List.metaClass.getKey = { ->
            if (delegate[0] instanceof Key) {
                return delegate[0]
            }
            throw GroovyRuntimeException('The leading element was not an instance of Key')
        }

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

        List.metaClass.pickWithoutKeys = { ... picks ->
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
                acc << delegate[0]
                owner.delegate[1..-1].each {
                    if (it instanceof NamedValue && it.name == pi) {
                        acc << it
                    }
                }
                acc
            }

        }

        /**
         * Add a key to a list which does not currently have one.
         */
        List.metaClass.addKey = {
            assert ! delegate[0] instanceof Key : 'This list already has a leading key/index element [${delegate[0]}]'
            Key key = new Key()
            delegate.each { key.put(it) }
            [key, *delegate]
        }

        /**
         * Update an existing key
         */
        List.metaClass.addOrUpdateKey = {
            Key key = new Key()
            int start = delegate[0] instanceof Key ? 1 : 0
            delegate[start..-1].each{ key.put(it) }
            [key, *(delegate[start..-1])]
        }

        /**
         * Add a name to a given list element, by wrapping it inside a NamedValue object
         */
        List.metaClass.nameify = { ix, name ->
            def val = delegate[ix]
            delegate[ix] = new NamedValue(name, val)
            return delegate
        }

        List.metaClass.unwrap = {
            delegate.collect { (it instanceof NamedValue ? it.value : it) }
        }
    }

    static class DataflowUtils {

        /**
         * Cross this channel with another, but flatten the resulting rows
         * and keep only one key element in the first position.
         */
        DataflowQueue flatCross(DataflowQueue b) {
            return this.cross(b).map{ t -> [*t[0], t[1][1..-1]].flatten() }
        }

    }

    /**
     * NamedValue wrap objects with a name.
     */
    @AutoClone
    static class NamedValue implements Comparable<NamedValue> {
        String name
        Object value

        /**
         * A default constructor is required for Kryo deserialization.
         */
        NamedValue() {/*...*/}

        /**
         * Basic constructor.
         * @param name - variable name
         * @param value - value to wrap
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

        int compareTo(NamedValue other) {
            return value.compareTo(other.value)
        }
    }

    @AutoClone
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
         * Construct a key from a given file. It is assume this file has been named using the delimited
         * convention used in sweeps and that the same separators are in use.
         * @param keyedFile -- a file whose name is a set of delimtied key/value pairs.
         * @param numSuffix -- the number of '.' appended suffixes to find and remove
         */
        Key(Path keyedFile, int numSuffix=1) {
            // make sure string is treated literally
            String safe_kv_sep = Pattern.quote(KEYVALUE_SEP)
            // drop any suffix from file's name
            String keyStr = keyedFile.name
            for (int i=0; i<numSuffix; i++) {
                keyStr = dropSuffix(keyStr)
            }
            // doubly split the remainder, first for values then key/value pairs.
            // using they pairs, initialise the instance's map.
            keyStr.split(Pattern.quote(PARAM_SEP))
                    .inject(this.varMap) { map, it ->
                        def (k, v) = it.split(safe_kv_sep)
                        map << [(k): new NamedValue(k, v)]
                    }
        }

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
            assert old == null: "attempted duplicate insertion of key [$v.name]"
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

        String id() {
            varMap.collect { k, v -> "$k$KEYVALUE_SEP${simpleSafeString(v)}" }.join(PARAM_SEP)
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
         * Remove levels of the key from the right side (lowest level first).
         * @param numLevels -- the number of levels to remove
         * @return the Key after removing @{numLevels} levels.
         */
        public Key popLevels(int numLevels) {
            List values = varMap.values() as List
            assert values.size() > numLevels : "Attempted to pop of [$numLevels] greater than total depth [${values.size()}"
            Key k = new Key()
            k.putAll(values[0..-(numLevels+1)])
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
            assert level <= values.size() : "Requested cut point [$level] beyond end. Max [${values.size()}] levels"
            def keys = ['lo':new Key(), 'hi':new Key()]
            keys.lo.putAll(values[0..level-1])
            if (level < values.size()) {
                keys.hi.putAll(values[level..-1])
            }
            return keys
        }

        public Key selectedKey(String... names) {
            assert names.size() > 0 : 'Must have at least one element selected'
            Map m = varMap.subMap(names)
            assert m.size() > 0 : "No elements matched [$names]"
            assert m.size() == names.size() : "Not all specified key names [$names] existed in [${varMap.keySet()}]"
            Key k = new Key()
            k.putAll(m.values())
            return k
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

                hshr.putUnencodedChars(cl.prefix)
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

        @Override
        public int hashCode() {
            Objects.hash(hash())
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
        //
        public String prefix
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
        public String describe() {
            def l = [prefix, simpleSafeString(ancestor), simpleSafeString(donor), ntaxa] +
                    mapString(tree) + mapString(profile)
            l.flatten().join(PARAM_SEP)
        }

        @Override
        public String toString() {
            prefix
        }

        @Override
        public int hashCode() {
            return Objects.hash(this.prefix)
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
            def colnames = 'Name '.padRight(width+1) + 'Values '
            def hline = ''.padRight(colnames.size(),'-')
            desc.append(hline + '\n')
            desc.append(colnames + '\n')
            desc.append(hline + '\n')
            int num = varRegistry.inject(1){ acc, k, v ->
                def label = "${k}:"
                desc.append("${label.padRight(width+1)} $v\n")
                acc *= v.size()
                acc
            }
            desc.append(hline + '\n')
            desc.append("Total combinations: $num\n")
            return desc.toString()
        }

        public Collection permuteAll() {
            permute(values())
        }

        /**
         * Return the permutation of the variables, as referenced by their names.
         * @param varNames the variables to permut
         * @return all permutations as a list
         */
        public Collection permuteNames(Object... varNames) {
            def ln = varNames.collect()
            // delegated subMap
            ln = subMap(ln)
            permute(ln.values())
        }

        public Collection permuteLevels(int begin = 0, int end = -1) {
            def names = keySet() as List
            permute(subMap(names[begin..end]).values())
        }

        /**
         * Permute the list of lists, return all permutations
         * @param collections -- a list of lists to permute
         * @return all permutations
         */
        protected Collection permute(Collection collections) {
            if (collections == null || collections.isEmpty()) {
                return []
            } else {
                List res = []
                permuteImpl(collections, res, 0, [])
                return res
            }
        }

        /**
         * Recursive implementation for permutations(List, Collection)
         */
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

        /**
         * Permute the entire sweep.
         * @return the permutations as a channel
         */
        public  DataflowQueue permutedChannel() {
            List allVars = varRegistry.keySet() as List
            permutedChannel(*allVars)
        }

        /**
         * Create channel embodying the permutations for the specified
         * variables of a sweep.
         * @param varNames the specified variables to include
         * @return the permutations as a channel
         */
        public DataflowQueue permutedChannel(Object... varNames) {
            def pn = permuteNames(varNames)
            Channel.from(pn).map{ it ->
                try {
                    it.addOrUpdateKey()
                }
                catch (Throwable t) {
                    println "ITER $it"
                    println "EXCEPTION: $t"
                    throw ex
                }
            }
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

        /**
         * Given an existing channel, extend it by those sweep variables specified
         * by their names. These variables must already have been defined within
         * the sweep instance.
         *
         * @param df -- the channel to extend
         * @param varNames -- the variables, referenced by their sweep names
         * @return an extended channel
         */
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

    }

    /**
     * Extend an existing list key (index) by adding the supplied key/value pair.
     *
     * This also updates the sweep definition.
     *
     * @param key -- target key to extend
     * @param name -- variable name to add
     * @param value -- value corresponding to this variable
     * @return a new instance of Key
     */
    public Key extendKey(Key key, String name, Object value) {
        Class clazz
        if (variables[name]) {
            clazz = variables[name][0].getClass()
        }
        else {
            variables[name] = []
        }
        variables[name].add(clazz ? clazz.cast(value) : value)
        withVariable(name)
        Key newKey = key.clone()
        newKey.put(name, value)
        return newKey
    }

    /**
     * Read a configuration from file and return an instance of MetaSweeper
     * @param config -- YAML format configuration file.
     * @return new MetaSweeper instance
     */
    static MetaSweeper fromFile(File config) {
        ConfigLoader.read(config)
    }

    /**
     * Print a table describing the current sweep.
     * @param msg -- an additional message to include
     */
    public void describeSweep(String msg=null) {
        if (msg) {
            println msg
        }
        println this.sweep.description()
    }

    /**
     * Initialize the sweep.
     * @return This MetaSweeper instance
     */
    public MetaSweeper createSweep() {
        sweep = new Sweep()
        return this
    }

    /**
     * Extend a channel by adding the named configuration variable to the sweep and
     * then permute the passed Channel.
     *
     * @param df -- the existing channel to extend
     * @param varName -- the named variable with which to extend the channel
     * @return the extended Channel
     */
    public DataflowChannel extend(DataflowQueue df, String... varName) {
        sweep.extendChannel(df, varName)
    }

    /**
     * Add a configuration variable to sweep of this MetaSweeper instance.
     * @param name - the configuration variable to add
     * @param stepInto - for non-collections which are iterable, step into the object and retrieve
     * the elements instead
     * @return MetaSweeper instance
     */
    public MetaSweeper withVariable(String name, Boolean stepInto=false) {
        if (stepInto) {
            // step into a variable which itself iterable but not a plain collection
            sweep[name] = variables[name].collect()
        }
        else {
            sweep[name] = variables[name]
        }
        return this
    }

    /**
     * For a given sweep defintion, create a channel which represents all the permutations.
     *
     * @return channel of perumuted sweep parameters.
     */
    public DataflowQueue permute() {
        sweep.permutedChannel()
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
    public static Map groupAtLevel(DataflowQueue df, int level) {
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

        // rejoin keys and return rows as [key, output_values]
        ret = ret.collect { row ->
            Key joined = sweep.orderKey(row[0] + row[1][0] + row[1][1])
            [joined, *row[2]]
        }

        return Channel.from(ret)
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
        return s.replaceAll(charsToReplace, "_")
    }

    /**
     * Append a suffix to a given key string. This will use
     * the defined separator to join the key string and new suffix.
     * @param keyString -- the key string to append to
     * @param suffix -- the suffix to append to key
     * @return joined keyString and suffix
     */
    static String appendKey(String keyString, String suffix) {
        "${keyString}${PARAM_SEP}${suffix}"
    }

    static DataflowQueue appendKey(DataflowChannel channel, String suffix) {
        channel.map{ f -> [appendKey(f[0], suffix), *f[1..-1]]}
    }

    /**
     * Convert a whitespace or comma delimited String to a List.
     * @param str - the string to split
     * @return List
     */
    static String[] stringToList(String str) {
        return str.split(delimiters) - ''
    }

    /**
     * Create a channel with a leading key derived from the file name. This assumes
     * the sweep variable delimited filename syntax is being used by the file.
     * @param path sweep file or list of files
     * @return Channel with a leading key for use in joins
     */
    public static DataflowChannel simpleKeyedFrom(path, int numSuffix=1) {
        assert path instanceof Path || path instanceof List<Path> : 'Error: supplied path must be an instance of ' +
                'java.nio.Path or Collection<Path>'

        Channel.from(path).map { f -> [new Key(f, numSuffix), f] }
    }

    /**
     * Attempts to cast a string first to an integer then to a double.
     * If both cause a NumberFormatException, the string is returned unchanged.
     * @param s -- the string to cast
     * @return Integer, Double or String representation of the input.
     */
    static def implicitValueOf(s) {
        if (! s instanceof String) {
            // we only want to try and cast strings
            return s
        }
        try {
            return Integer.valueOf(s)
        } catch (NumberFormatException ex) {/*...*/}
        try {
            return Double.valueOf(s)
        } catch (NumberFormatException ex) {/*...*/}
        return s
    }

    /**
     * Use implicit casting on an entire collection. It is assume that a collection should
     * be only one type. For cases where more than one type has been returned by implicit
     * conversion (such as 1 returning as an Integer, when other numbers in the list were
     * Doubles), the list will be cast to the most specific type.
     *
     * An exception is when Strings are mixed with non-strings. In this case, the original
     * collection is returned.
     *
     * @param c -- the collection to cast
     * @return the cast collection
     */
    static def implicitValueOf(Collection stringValues) {
        def castValues = stringValues.collect{ implicitValueOf(it) }
        def castTypes = castValues.inject([] as Set){ acc, it -> acc.add(it.getClass()); acc }
        if (castTypes.size() == 1) {
            // were done, everything is one type
            return castValues
        }
        else {
            if (Double.class in castTypes) {
                // when doubles exist, everything is converted to doubles.
                return castValues.collect { (Double) it }
            }
            else {
                // any other mixed-type case is presently abandoned.
                return stringValues as List
            }
        }
    }

    /**
     * Create a channel using the keys as found on a path or list of paths.
     * @param path -- one or a list of @{link java.nio.Path}
     * @param numSuffix -- the number of dot '.' suffixes to remove, leaving only the key data.
     * @param permute -- if true, create a permutation of all parameters. Otherwise just return those
     * iterations which were observed.
     * @return a channel of iterations, possibly part or all of a sweep depending on the paths parsed
     */
    public DataflowChannel keyedFrom(path, int numSuffix=1, boolean permute=false) {
        assert path instanceof Path || path instanceof List<Path> : 'Error: supplied path must be an instance of ' +
                'java.nio.Path or Collection<Path>'

        if (! sweep) {
            createSweep()
        }

        List files = Channel.from(path).toList().get()
        def toAdd = [:].withDefault{[] as Set}

        // Get the key from each file name, collect a map of variables and values
        def keyedFiles = files.collect { f ->
            Key key = new Key(f, numSuffix)
            key.varMap.collect{ k, v ->
                toAdd[k].add(v.value)
            }
            [key, f]
        }

        // convert string representations to more specific types if possible. Override any
        // existing definition of that variable or create a new one. Finally, add the variable
        // to the sweep.
        toAdd.each { k, v ->
            variables[k] = implicitValueOf(v).sort(false)
            withVariable(k)
        }

        // create a channel which is a permutation of all the discovered variables and their values.
        if (permute) {
            return permute()
        }
        else {
            return Channel.from(keyedFiles)
        }
    }

    /**
     * Remove the right-most dot '.' suffix from a String
     * @param str -- the string to parse
     * @return string without a suffix if found
     */
    static String dropSuffix(str) {
        return str.lastIndexOf('.').with {it != -1 ? str[0..<it] : str}
    }

    /**
     * Create the cartesian product of two collections
     * @param A -- collection A
     * @param B -- collection B
     * @return cartesian product as an array of objects
     */
    static Object[] product(Collection A, Collection B) {
        return A.collectMany{a->B.collect{b->[a, b]}}
    }

    static Yaml getYamlParser() {
        DumperOptions options = new DumperOptions();
        options.setWidth(MAX_LINE_LENGTH)
        return new Yaml(options)
    }

}
