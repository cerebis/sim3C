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

class Helper {
    static String SEPARATOR = '-+-'
    static def delimiters = /[ ,\t]/

    static def absPath(path) {
        Nextflow.file(path)*.toAbsolutePath()
    }

    static def str(val, safe=true) {
        if (val instanceof java.nio.file.Path) {
            val = val.name - ~/[\.][^\.]+$/
        }
        else {
            val = val
        }
        return safe ? safeString(val) : val
    }

    static def addKey(row, firstCol=0, lastCol=-1) {
        [row[firstCol..lastCol].collect { str(it) }.join(SEPARATOR), *row]
    }

    static def select(row, elems) {
        row[elems]
    }

    static int[] stringToInts(String str) {
        return (str.split(delimiters) - '').collect { elem -> elem as int }
    }

    static float[] stringToFloats(String str) {
        return (str.split(delimiters) - '').collect { elem -> elem as float }
    }

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
