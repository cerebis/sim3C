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
import java.nio.file.Path

class Helper {
    static def separators = /[ ,\t]/

    static int[] stringToInts(String str) {
        return (str.split(separators) - '').collect { elem -> elem as int }
    }

    static float[] stringToFloats(String str) {
        return (str.split(separators) - '').collect { elem -> elem as float }
    }

    static String[] stringToList(String str) {
        return str.split(Helper.separators) - ''
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
        return name.split(Globals.separator)[0..-(n+1)].join(Globals.separator)
    }

    static String removeLevels(String name, int n) {
        return name.split(Globals.separator)[0..-(n+1)].join(Globals.separator)
    }

    static Object[] product(A, B) {
        return A.collectMany{a->B.collect{b->[a, b]}}
    }
}
