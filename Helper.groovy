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
