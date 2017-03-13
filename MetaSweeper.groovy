package mzd
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

@Grab('com.google.guava:guava:19.0')
@Grab('org.yaml:snakeyaml:1.17')
import groovy.transform.AutoClone
import groovy.transform.Synchronized
import groovyx.gpars.dataflow.DataflowChannel
import groovyx.gpars.dataflow.DataflowQueue
import java.nio.file.Path
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
        KryoHelper.register(Community)
        KryoHelper.register(Clade)

        /**
         * Get just the key from the list (table row)
         */
        List.metaClass.getKey = { ->
            assert delegate[0] instanceof Key : 'Leading element was not an instance of Key'
            delegate[0]
        }

        List.metaClass.dropKey = { ->
            assert delegate[0] instanceof Key : 'Leading element was not an instance of Key'
            delegate[1..-1]
        }

        /**
         * Convenience method for picking out list elements by index, where lists are conceptualised
         * as rows in a table, each containing a leading @{link Key} element.
         *
         * The method always returns the first element and tests that it is an instance of Key.
         * <br>
         * Both integer indexes and lists of such integers can be supplied together. Lists imply that the
         * returned elements should be wrapped in a list.
         *
         * Supplying no indices still returns the 1 element list containing the key.
         *
         * Ranges can be supplied, but will be interpreted as a list of integers and so be returned
         * wrapped in a list.
         */
        List.metaClass.pick = { Object... indices ->
            def key = delegate[0]
            assert key instanceof Key : 'First element of row is not a Key. Perhaps unwrap nested elements first'
            def out = [key]
            indices.each { rr ->
                if (rr instanceof List) {
                    assert rr.any{ ri -> ri instanceof Integer } : 'List contains non-integers'
                    out << delegate[rr]
                }
                else if (rr instanceof Integer) {
                    out << delegate[rr]
                }
                else {
                    throw new GroovyRuntimeException('arguments may ony be a combination of List<Integer> or Integer')
                }
            }
            out
        }

        List.metaClass.pickWithoutKey = { Object... indices ->
            indices.collect { rr ->
                if (rr instanceof List) {
                    assert rr.any{ ri -> ri instanceof Integer } : 'List contains non-integers'
                    delegate[rr]
                }
                else if (rr instanceof Integer) {
                    delegate[rr]
                }
                else {
                    throw new GroovyRuntimeException('arguments may ony be a combination of List<Integer> or Integer')
                }
            }
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
     * A {@link Key} acts as an identifier for a single iteration within a parameter sweep.
     *
     * Within {@link MetaSweeper} channels are best conceptualised as tables. Each element in
     * a channel element corresponds to a row and each element within a row corresponds to a
     * column. MetaSweeper attempts to maintan a key as the first element of each row, providing
     * a means of recording variable state for each result in the workflow and permit the joining
     * of output channels which split at an earlier stage.
     *
     * As the depth of the sweep varies, the number of elements stored in the Key varies in accordance.
     * When represented as a string, keys are doubly delimited, first as variable-name/variable-value
     * pairs for each active variable, and then as the set of active variables.
     *
     * The presently used delimiters have been chosen so as to avoid conflicts with common filesystems,
     * but provide a basic degree of human-readability.
     *
     * Eg. "seed#123-+-xfold#10-+-n3c#100000"
     */
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
                        map.put(k, v)
                        map
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
        Map map() {
            this.varMap
        }

        @Synchronized
        List values() {
            return this.varMap.values()
        }

        @Synchronized
        List keyList() {
            return this.varMap.keySet() as List
        }

        @Synchronized
        Key copy() {
            Key newKey = new Key()
            this.varMap.each{ k, v -> newKey.put(k, v)}
            return newKey
        }

        @Synchronized
        Map putAll(Collection kv) {
            def ret = [:]
            v.each{ ret[it.name] = put(it) }
            return ret
        }

        @Synchronized
        Object put(String k, Object v) {
            assert this.varMap.putIfAbsent(k, v) == null : "attempted duplicate insertion of key [$k]"
        }

        @Synchronized
        Object getAt(String k) {
            return this.varMap[k]
        }

        String toString() {
            return id()
        }

        String id() {
            return this.varMap.collect { k, v ->
                "$k$KEYVALUE_SEP${simpleSafeString(v)}" }.join(PARAM_SEP)
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

        List join(Key other) {
            return keyList().intersect(other.keyList()) as List
        }

        /**
         * Return a subkey, where the depth is less than the maximum depth of the sweep.
         * @param depth - the depth to retain
         * @return reduced key
         */
        Key subKey(int depth) {
            assert depth > 0 : 'Depth must be a positive integer'
            List keys = keyList()
            assert depth <= keys.size() : "Requested depth [$depth] exceeds defined variable count [${keys.size()}]"
            Key k = new Key()
            k.varMap << this.varMap.subMap(keys[0..(depth-1)])
            return k
        }

        /**
         * Pop a number of levels from the key on the right side (lowest level first).
         * @param numLevels -- the number of levels to remove
         * @return the Key after removing @{numLevels} levels.
         */
        Key popLevels(int numLevels) {
            assert numLevels >= 0 : 'Number of levels must be >= 0'
            assert numLevels < this.varMap.size() : 'Action would remove all levels from key'
            return subKey(this.varMap.size() - numLevels)
        }

        /**
         * At a given sweep level, split a key into two parts.
         *
         * Neither part can be empty, that is you cannot split at the lowest or highest level.
         *
         * @param level - the sweep level at which to split (0 > level < max_level)
         * @return Map containing two keys: 'hi' and 'lo'
         */
        Map splitKey(int level) {
            List keys = keyList()
            assert level > 0 : 'Depth must be > 0'
            assert level <= keys.size() : "Requested cut point [$level] beyond end. Max [${keys.size()}] levels"
            def splitKey = ['lo': new Key(), 'hi': new Key()]
            splitKey.lo.varMap << this.varMap.subMap(keys[0..<level])
            splitKey.hi.varMap << this.varMap.subMap(keys[level-1..-1])
            return splitKey
        }

        Key selectedKey(String... names) {
            assert names.size() > 0 : 'Must have at least one element selected'
            Map m = varMap.subMap(names)
            assert m.size() > 0 : "No elements matched [$names]"
            assert m.size() == names.size() : "Not all specified key names [$names] existed in [${varMap.keySet()}]"
            Key key = new Key()
            m.each { name, value ->
                key.put(name, value)
            }
            return key
        }

    }

    /**
     * Primary class representing a parameter sweep.
     */
    static class Sweep {
        @Delegate
        Map varRegistry = [:]

        /**
         * Introduce a new variable and its values to the sweep.
         *
         * @param name - the variable name
         * @param values - the values associated with this variable
         * @return this {@link Sweep}
         */
        Sweep withVariable(String name, Object values) {
            put(name, values)
            return this
        }

        /**
         * Drop (or remove) a variable from this instance of {@link Sweep}.
         * @param name - the variable name
         * @return this {@link Sweep} without the variable
         */
        Sweep dropVariable(String name) {
            varRegistry.remove(name)
            return this
        }

        /**
         * Append a value to a variable in sweep. If the variable doesn't exist,
         * then a new entry is created.
         * @param name - the variable name in the sweep
         * @param value - the value to add to this variable
         * @return this instance of {@link Sweep}
         */
        Sweep appendValue(String name, Object value) {
            if (name in varRegistry) {
                varRegistry[name].add(value)
            }
            else {
                put(name, value)
            }
            return this
        }

        /**
         * Add a new entry to the variable registry. If the supplied values are not
         * a collection, then first wrap it with a list.
         *
         * @param key - the variable name for this entry
         * @param value - the value(s) for this variable
         * @return the added values
         */
        Object put(Object key, Object value) {
            if (!(value instanceof Collection)) {
                varRegistry[key] = [value]
            }
            else {
                varRegistry[key] = value.collect()
            }
        }

        /**
         * Print a tabular description of the present state of the sweep.
         * @param title - an optional title
         * @return this {@link Sweep}
         */
        Sweep describe(String title=null) {

            if (title) {
                println title
            }
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
            println desc.toString()

            return this
        }

        /**
         * Generate a permutation of all variables in the sweep.
         * @return the permutation as a collection
         */
        Collection permuteAll() {
            permute(values())
        }

        /**
         * Return the permutation of the variables, as referenced by their names.
         * @param varNames the variables to permut
         * @return all permutations as a collection
         */
        Collection permuteNames(Object... varNames) {
            def ln = varNames.collect()
            // delegated subMap
            ln = subMap(ln)
            // make name/value pairs for each variable, then permute the pairs.
            permute(ln.collect{ k, v -> GroovyCollections.combinations([k], v) })
        }

        /**
         * Carry-out the permutation between a range of levels within the sweep.
         * @param begin - the first level at which to begin
         * @param end - the last level, at which to stop
         * @return a permutation as a collection
         */
        Collection permuteLevels(int begin = 0, int end = -1) {
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
        DataflowQueue permutedChannel() {
            List allVars = varRegistry.keySet() as List
            permutedChannel(*allVars)
        }

        /**
         * Create channel embodying the permutations for the specified
         * variables of a sweep.
         * @param varNames the specified variables to include
         * @return the permutations as a channel
         */
        DataflowQueue permutedChannel(Object... varNames) {
            def pn = permuteNames(varNames)

            Channel.from(pn).map{ row ->
                try {
                    row.inject([new Key()]) { row_out, it ->
                        def (name, val) = it
                        row_out[0].put(name, val)
                        row_out << val
                        row_out
                    }
                }
                catch (Throwable t) {
                    println "ITER $row"
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
        Key orderKey(Key key) {
            return varRegistry.inject(new Key()) { ordKey, name, values ->
                if (key[name]) {
                    ordKey.put(name, key[name])
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
        DataflowQueue extendChannel(DataflowQueue df, String... varNames) {
            assert varNames.size() > 0 : "Error, no variables named in extendChannel"

            def p = permuteNames(varNames)
            if (p.size() == 0) {
                println 'Warning: the specified variables may not be referenced. Permutation length was 0.'
                return df
            }

            df = df.spread(p).map { row ->
                // the key
                def key = row[0].copy();
                // the output row begins with the key adn preexisting variables
                def row_out = [key] + row[1..<-(2*varNames.size())].collect()
                // the new variables, we'll pair up the name/values too
                def new_vars = pairListElements( row[-(2*varNames.size())..-1] )
                // add each new variable to the key, add append it to the row
                new_vars.each { it ->
                    def (name, val) = it
                    key.put(name, val)
                    row_out << val
                }
                row_out
            }

            return df
        }

        /**
         * Assuming an ordered list of successive pairs, [A1,A2,B1,B2], combine each pair
         * as a 2-element list.
         *
         * f([A1,A2,B1,B2]) -> [[A1,A2], [B1,B2]]
         *
         * @param list - flat interleaved list of pairs
         * @return nested list of pairs.
         */
        static List pairListElements(List list) {
            (0..<list.size()/2).collect{ list[it*2..<(it+1)*2] }
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
         * and followed by all the per-task output values which were emitted by both channels. As such, it can be
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
        DataflowQueue joinChannels(DataflowQueue df1, DataflowQueue df2, int level) {
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
                Key joined = orderKey(row[0] + row[1][0] + row[1][1])
                [joined, *row[2]]
            }

            return Channel.from(ret)
        }


        /**
         * Fork the contents of an existing channel into multiple output channels,
         * where there will be one output channel per-value of the specified
         * sweep variable. Matches are complete (^..$) on the values.
         *
         * @param df -- DataflowChannel to fork
         * @param name -- Sweep variable name to fork on
         * @return the split channel as a list of channels.
         */
        Map<String, DataflowChannel> forkOnVariable(DataflowChannel df, String name) {
            assert name in varRegistry : "Error: ${name} is not a registered variable name"
            def values = this.varRegistry[name]
            // create the empty destination channels
            def forkedChans = values.collect { Channel.create() }
            // the result is returned a map, named by the values for the specified variable
            def chanMap = values.inject([:]) { acc, it -> acc[it] = []; acc}
            // fork the channel
            df.choice( forkedChans ) {
                def key = it.getKey()
                def idx = values.findIndexOf { vi -> key[name] =~ /$vi/ }
                assert idx != -1 : "Error: failed to find ${vi} in key ${key}"
                idx
            }
            // put the split channels into map. Assumes order is preserved
            values.eachWithIndex { it, n -> chanMap[it] = forkedChans[n] }
            chanMap
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
        Key extendKey(Key key, String name, Object value) {
            appendValue(name, value)
            Key newKey = key.clone()
            newKey.put(name, value)
            return newKey
        }

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
     * Community represents a collection of clades/groups of related genomes.
     * Together it can be regarded as a simulated community.
     */
    static class Community implements Iterable {
        // just a human centric name for a community
        public String name
        // the list of clades/groups to generate
        public List<Clade> clades
        // parameters involved in generating the profile
        public Map<String, Float> profile

        Iterator iterator() {
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
            OutputStream output = Funnels.asOutputStream(hasher)
            try {
                Files.copy(file.toPath(), output)
            }
            catch( IOException e ) {
                throw new IllegalStateException("Unable to hash content: ${file}", e)
            }
            finally {
                if (output) {
                    output.close()
                }
            }
            return hasher
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
        String hash() {

            HashFunction hf = Hashing.murmur3_32()
            Hasher hshr = hf.newHasher()

            // include the definition of profile
            profile.collect { e ->
                hshr.putUnencodedChars(e.getKey())
                        .putFloat(e.getValue())
            }

            // hash clades definitions
            clades.each { cl ->
                // deeply hash clade's ancestor and donor
                hashFileContent(hshr, cl.ancestor)
                hashFileContent(hshr, cl.donor)

                hshr.putUnencodedChars(cl.prefix)
                hshr.putInt(cl.ntaxa)

                // include the definition of clade's tree
                cl.tree.collect { e ->
                    hshr.putUnencodedChars(e.getKey())
                    if (e.value instanceof String)
                        hshr.putUnencodedChars(e.getValue())
                    else if (e.value instanceof Float)
                        hshr.putFloat(e.getValue())
                    else if (e.value instanceof Double)
                        hshr.putDouble(e.getValue())
                    else if (e.value instanceof Integer)
                        hshr.putInt(e.getValue())
                    else if (e.value instanceof Long)
                        hshr.putLong(e.getValue())
                    else
                        throw new RuntimeException('Unsupported type')
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
            "${name}_${hash()}".toString()
        }

        @Override
        int hashCode() {
            Objects.hash(hash())
        }

    }

    /**
     * A Clade represents set of related genomes evolved from
     * a common ancestor. How these sequences are different
     * is dictated by a phylogenetic tree.
     *
     * Both trees are expected to be generated either through simulation -- which
     * itself is defined by a set of parameters -- or defined explicitly either
     * as a string directly in the YAML file or as a reference to an external
     * URL (file, http, etc).
     *
     * How many taxa in exist in the clade is defined by ntaxa.
     */
    static class Clade {
        private static final String REF_TAG = 'ref'
        private static final String ALGO_TAG = 'algo'

        // a short id or name
        public String prefix
        // common ancestor
        private Path ancestor
        // donor used for htg
        private Path donor
        // number of leaves/tips.
        public Integer ntaxa
        // parameters involved in generating the tree
        public Map<String, Object> tree

        void setAncestor(String ancestor) {
            this.ancestor = Nextflow.file(ancestor)
        }
        String getAncestor() {
            return ancestor.toString()
        }
        Path getAncestorPath() {
            return ancestor
        }

        void setDonor(String donor) {
            this.donor = Nextflow.file(donor)
        }
        String getDonor() {
            return donor.toString()
        }
        Path getDonorPath() {
            return donor
        }

        /**
         * Check that the algorithm is supported by the workflow.
         * Note: At present we only support birth_death
         */
        boolean isSupportedAlgorithm() {
            ALGO_TAG in tree && tree[ALGO_TAG] == 'birth_death'
        }

        /**
         * Test whether this clade uses a user defined tree.
         *
         * Checks for the existence of {@link Clade#REF_TAG} in the tree map.
         */
        boolean isDefined() {
            REF_TAG in tree
        }

        /**
         * Test whether this clade uses a procedurally generated tree.
         *
         * Checks for the existence of {@link Clade#ALGO_TAG} in the tree map.
         */
        boolean isProcedural() {
            ALGO_TAG in tree
        }

        /**
         * Return the user defined NEWICK format tree as a string.
         *
         * Tree definitions may be provided either as a plain string in the YAML config.
         * or as a URL.
         *
         * When using a URL, please note that a http-get is performed each time this method
         * is called. Therefore, attention should be paid to not make redundant calls to
         * this method in cases where a remote server might interpret repeated calls as
         * something to block.
         *
         * Examples.
         *
         * ref: "(A:0.5, B:0.1);"
         * ref: "file:a/relative/path/tree.nwk"
         * ref: "file:/an/absolute/path/tree.nwk"
         * ref: "http://remote.server.io/a/tree.nwk"
         *
         * @return a string representaing a NEWICK format tree
         */
        String getDefined() {
            try {
                tree[REF_TAG].toURL().getText()
            }
            catch (MalformedURLException ex) {
                tree[REF_TAG]
            }
        }

        /**
         * String representation of the Clade as a concatenated set of its parameters
         * @return String
         */
        String describe() {
            def l = [prefix, simpleSafeString(ancestor), simpleSafeString(donor), ntaxa] + mapString(tree)
            l.flatten().join(PARAM_SEP)
        }

        @Override
        String toString() {
            prefix
        }

        @Override
        int hashCode() {
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
        static MetaSweeper read(File file) {
            Yaml yaml = threadLocal.get()
            yaml.load(new FileReader(file))
        }

        /**
         * Read a configuration from a string
         * @param doc -- String to read
         * @return MetaSweeper
         */
        static Map read(String doc) {
            Yaml yaml = threadLocal.get()
            yaml.load(doc)
        }

        /**
         * Convert MetaSweeper object to a string in YAML syntax
         * @param ms -- Object to convert
         * @return String (YAML)
         */
        static String write(MetaSweeper ms) {
            Yaml yaml = threadLocal.get()
            yaml.dump(ms)
        }

        /**
         * Write a YAML representation of a MetaSweeper object to a stream
         * @param ms - Object to write
         * @param writer - output writer
         */
        static void write(MetaSweeper ms, Writer writer) {
            Yaml yaml = threadLocal.get()
            yaml.dump(ms, writer)
        }
    }

    /**
     * Initialize a sweep.
     * @return This MetaSweeper instance
     */
    static Sweep createSweep() {
        new Sweep()
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
     * Create a sweep, initialized using the keys as found on a path or list of paths. This
     * also overrides variables with the same entry keys in this MetaSweeper instance.
     *
     * @param sweep -- sweep instance to initialize from path parsing
     * @param path -- one or a list of @{link java.nio.Path}
     * @param numSuffix -- the number of dot '.' suffixes to remove, leaving only the key data.
     * @param permute -- create a permutation channel from all variables
     * @return a channel of iterations, possibly part or all of a sweep depending on the paths parsed
     */
    DataflowChannel keyedFrom(sweep, path, int numSuffix=1, boolean permute=false) {
        assert path instanceof Path || path instanceof List<Path> : 'Error: supplied path must be an instance of ' +
                'java.nio.Path or Collection<Path>'

        // why make a list when we shold be able to iterate on channel
        List files = Channel.from(path).toList().get()
        def toAdd = [:].withDefault{[] as Set}

        // Get the key from each file name, collect a map of variables and values
        def keyedFiles = files.collect { f ->
            Key key = new Key(f, numSuffix)
            key.varMap.collect{ k, v ->
                toAdd[k].add(v)
            }
            [key, f]
        }

        // convert string representations to more specific types if possible. Override any
        // existing definition of that variable or create a new one. Finally, add the variable
        // to the sweep.
        toAdd.each { k, v ->
            // add to set of defined variables for this instance
            variables[k] = implicitValueOf(v).sort(false)
            // also add immediately to sweep
            sweep.withVariable(k, variables[k])
        }

        // create a channel which is a permutation of all the discovered variables and their values.
        if (permute) {
            return sweep.permutedChannel()
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
