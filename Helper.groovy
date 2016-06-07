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
import nextflow.Nextflow
import java.nio.file.Path

/**
 * Utility methods used in the construction of meta-sweep workflows.
 */
class Helper {
    static String SEPARATOR = '-+-'
    static def delimiters = /[ ,\t]/

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
     * Create string representation of value. If this is an instance of {@link java.nio.file.Path},
     * then use only the file name. In the context of filesystem restrictions, mask problematic
     * characters if safe=true.
     *
     * @param val - object to represent as a string
     * @param safe - a flag that when true masks problematic characters
     * @return String
     */
    static def str(val, safe=true) {
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
        [row[firstCol..lastCol].collect { str(it) }.join(SEPARATOR), *row]
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
    static def select(row, elems) {
        row[elems]
    }

    /**
     * Return a safer representation of a string. Where "safe" is in the context
     * of filesystem restrictions. Therefore, masking of backslash characters and
     * the use only of the base file name when provided with an instance of
     * {@link java.nio.file.Path}
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
