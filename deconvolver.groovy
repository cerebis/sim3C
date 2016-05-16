class Globals {
    static String separator = '%'
}

trees = file(params.trees).collectEntries{ [it.name, it] }
tables = file(params.tables).collectEntries{ [it.name, it] }
ancestor = file(params.ancestor)
donor = file(params.donor)

alpha_BL = Helper.stringToList(params.alpha)
xfold = Helper.stringToList(params.xfold)
nhic = Helper.stringToList(params.hic_pairs)

/**
 * Helper methods
 */

import groovyx.gpars.dataflow.DataflowQueue

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

class ChannelDuplicator {
    DataflowQueue orig

    ChannelDuplicator(DataflowQueue orig) {
        this.orig = orig
    }

    DataflowQueue onCopy() {
	    def copied, keep
        (copied, keep) = this.orig.into(2)
        this.orig = keep
        return copied
    }

    static ChannelDuplicator createFrom(Object[] o) {
        return new ChannelDuplicator(Channel.from(o))
    }

    static ChannelDuplicator createFrom(DataflowQueue q) {
        return new ChannelDuplicator(q)
    }
}


/**
 * Generate simulated communities
 */

evo_sweep = Channel
        .from([ancestor])
        .spread([donor])
        .spread(alpha_BL)
        .spread(trees.values())
        .map { it += it.collect { Helper.safeString(it) }.join(Globals.separator) }

process Evolve {
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'false'

    input:
    set file('ancestral.fa'), file('donor.fa'), alpha, file('input_tree'), oname from evo_sweep

    output:
    set file("${oname}.evo.fa") into descendents

    """
    scale_tree.py -a ${alpha} input_tree scaled_tree
    sgEvolver --indel-freq=${params.indel_freq} --small-ht-freq=${params.small_ht_freq} --large-ht-freq=${params.large_ht_freq} \
         --inversion-freq=${params.inversion_freq} --random-seed=${params.seed} scaled_tree \
         ancestral.fa donor.fa "${oname}.evo.aln" "${oname}.evo.fa"
    strip_semis.sh "${oname}.evo.fa"
    """
}

/**
 * Generate WGS read-pairs
 */
descendents = ChannelDuplicator.createFrom(descendents)

wgs_sweep = descendents.onCopy()
        .spread(tables.values())
        .spread(xfold)
        .map{ it += tuple(it[0].name[0..-8], it[1].name, it[2]).join(Globals.separator) }

process WGS_Reads {
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'false'

    input:
    set file('descendent.fa'), file('profile'), xf, oname from wgs_sweep

    output:
    set file("${oname}.wgs*.fq.gz"), oname into wgs_reads

    """
    metaART.py -C gzip -t $profile -M $xf -S ${params.seed} -z ${params.num_samples} -s ${params.wgs_ins_std} \
            -m ${params.wgs_ins_len} -l ${params.wgs_read_len} -n "${oname}.wgs" descendent.fa .
    """
}


/**
 * Map WGS read-pairs to reference
 */

wgs_reads = ChannelDuplicator.createFrom( wgs_reads.map { reads, oname -> [*reads, oname] } )

map_sweep = wgs_reads.onCopy()

process ReadMap {
    cpus 1
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'false'

    input:
    set file('ancestral.fa'), file('*.wgs*.fq.gz'), oname from map_sweep

    output:
    set file("*.bam") into wgs_bams

    """
    bwa index ancestral.fa
    for r in `ls *.wgs.*r1.fq.gz`
    do
        rbase=`basename \$r .r1.fq.gz`
        r2=\$rbase.r2.fq.gz
        bwa mem -t 1 ancestral.fa \$r \$r2 | samtools view -bS - | samtools sort -l 9 - \$rbase.bam
    done
    """
}

/**
 * Deconvolve the SNVs into strain genotypes
 */

decon_bams = ChannelDuplicator.createFrom( wgs_bams.map { bams, oname -> [*bams, oname] } )
decon_sweep = decon_bams.onCopy()

process Deconvolve {
    cpus 1
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'false'

    input:
    set file('*.bam'), oname from decon_sweep

    output:
    set file("decon.csv"), file("strains.tre") into deconvolution

    """
    snvbpnmft.py . 4 ancestral.fa ${oname}.*.bam
    """
}
