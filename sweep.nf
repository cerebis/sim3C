class Globals {
    static String sweep_separator = '%'
}

//sweep_separator = '%'
alpha_BL = [1]
xfold = [1]
nhic = [10000]
trees = file('test/trees/*nwk').collectEntries{ [it.name, it] }
tables = file('test/tables/*table').collectEntries{ [it.name, it] }
ancestor = file('test/NC*raw').collectEntries{ [it.name, it] }
out_path = 'out'

seq_len = 3000000
sg_scale = 1e-4
seed = 1234
wgs_read_len = 150
wgs_ins_len = 450
wgs_ins_std = 100
hic_inter_prob = 0.9
hic_read_len = 150

/**
 * Helper methods
 */

import groovyx.gpars.dataflow.DataflowQueue

class Helper {
    static def separators = /[ ,\t]/

    static int[] stringToInts(String str) {
        return str.split(separators).collect { elem -> elem as int }
    }

    static float[] stringToFloats(String str) {
        return str.split(separators).collect { elem -> elem as float }
    }

    static String[] stringToList(String str) {
        return str.split(Helper.separators)
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
        return name.split(Globals.sweep_separator)[0..-(n+1)].join(Globals.sweep_separator)
    }

    static String removeLevels(String name, int n) {
        return name.split(Globals.sweep_separator)[0..-(n+1)].join(Globals.sweep_separator)
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
        .from(ancestor.values())
        .spread(alpha_BL)
        .spread(trees.values())
        .map { it += it.collect { Helper.safeString(it) }.join(Globals.sweep_separator) }

process Evolve {
    cache 'deep'
    publishDir out_path, mode: 'copy', overwrite: 'false'

    input:
    set file('ancestral.fa'), alpha, file('input_tree'), oname from evo_sweep

    output:
    set file("${oname}.evo.fa") into descendents

    """
    prepParams.py --tree input_tree --seq ancestral.fa --seq-len $seq_len --tree-scale $alpha \
        --sg-scale $sg_scale --out-name "${oname}.evo.fa"
    simujobrun.pl ancestral.fa $seed
    """
}

/**
 * Generate WGS read-pairs
 */
descendents = ChannelDuplicator.createFrom(descendents)

wgs_sweep = descendents.onCopy()
        .spread(tables.values())
        .spread(xfold)
        .map{ it += tuple(it[0].name[0..-8], it[1].name, it[2]).join(Globals.sweep_separator) }

process WGS_Reads {
    cache 'deep'
    publishDir out_path, mode: 'copy', overwrite: 'false'

    input:
    set file('descendent.fa'), file('profile'), xf, oname from wgs_sweep

    output:
    set file("${oname}.wgs*.fq"), oname into wgs_reads

    """
    metaART.py -t $profile -M $xf -S $seed -s $wgs_ins_std -m $wgs_ins_len -l $wgs_read_len -n "${oname}.wgs" descendent.fa .
    """
}

/**
 * Generate 3C/HiC read-pairs
 */

hic_sweep = descendents.onCopy()
        .spread(tables.values())
        .spread(nhic)
        .map{ it += tuple(it[0].name[0..-8], it[1].name, it[2]).join(Globals.sweep_separator) }

process HIC_Reads {
    cache 'deep'
    publishDir out_path, mode: 'copy', overwrite: 'false'

    input:
    set file('descendent.fa'), file('profile'), nh, oname from hic_sweep

    output:
    set file("${oname}.hic.fa") into hic_reads

    """
    simForward.py -r $seed -n $nh -l $hic_read_len -p $hic_inter_prob -t $profile -s descendent.fa -o "${oname}.hic.fa"
    """
}

/**
 * Assemble WGS read-pairs
 */

wgs_reads = ChannelDuplicator.createFrom( wgs_reads.map { reads, oname -> [*reads, oname] } )

asm_sweep = wgs_reads.onCopy()

process Assemble {
    cpus 1
    cache 'deep'
    publishDir out_path, mode: 'copy', overwrite: 'false'

    input:
    set file('wgs_R1.fq'), file ('wgs_R2.fq'), oname from asm_sweep

    output:
    file("${oname}.contigs.fasta") into wgs_contigs

    """
    a5_pipeline.pl --threads=1 --metagenome wgs_R1.fq wgs_R2.fq ${oname}
    """
}

/**
 * Generate contig->descendent gold standard mapping
 */

wgs_contigs = ChannelDuplicator.createFrom(wgs_contigs)

a = wgs_contigs.onCopy().map { f -> [Helper.removeLevels(f.name,2), f, f.name[0..-15]] }
tr_sweep = descendents.onCopy()
        .spread(tables.values())
        .map { t -> [t[0].name[0..-8], t[0]] }
        .phase(a){ t -> t[0] }.map { t -> [*t[0], *t[1][1..2]] }

process Truth {
    cache 'deep'
    publishDir out_path, mode: 'copy', overwrite: 'false'

    input:
    set key, file('ref.fa'), file('contigs.fa'), oname from tr_sweep

    output:
    file("${oname}.truth.yaml") into truth_tables

    """
    if [ ! -e db.prj ]
    then
        lastdb db ref.fa
    fi

    lastal -P 1 db contigs.fa | maf-convert psl > ctg2ref.psl
    alignmentToTruth.py --afmt test --ofmt yaml ctg2ref.psl "${oname}.truth.yaml"
    """
}

/**
 * Map 3C/HiC read-pairs to assembly contigs
 */

a = wgs_contigs.onCopy()
        .map { f -> [Helper.removeLevels(f.name,1), f] }

hicmap_sweep = hic_reads
        .map { f -> [Helper.removeLevels(f.name, 1), f, f.name[0..-8]] }
        .phase(a){ t -> t[0] }.map { t -> [t[0][1], t[1][1], t[0][2]] }

process HiCMap {
    cache 'deep'
    publishDir out_path, mode: 'copy', overwrite: 'false'

    input:
    set file('hic.fa'), file('contigs.fa'), oname from hicmap_sweep

    output:
    file("${oname}.hic2ctg*") into hic2ctg_mapping


    """
    bwa index contigs.fa
    bwa mem -t 1 contigs.fa hic.fa | samtools view -bS - | samtools sort - "${oname}.hic2ctg"
    samtools index "${oname}.hic2ctg.bam"
    samtools idxstats "${oname}.hic2ctg.bam" > "${oname}.hic2ctg.idxstats"
    samtools flagstat "${oname}.hic2ctg.bam" > "${oname}.hic2ctg.flagstat"
    """
}

/**
 * Generate contig graphs from HiC mappings
 */

hic2ctg_mapping = ChannelDuplicator.createFrom(hic2ctg_mapping)

graph_sweep = hic2ctg_mapping.onCopy()
        .map {t -> [t[0], t[1], t[0].name[0..-13] ] }

process Graph {
    cache 'deep'
    publishDir out_path, mode: 'copy', overwrite: 'false'

    input:
    set file('hic2ctg.bam'), file('hic2ctg.bam.bai'), oname from graph_sweep

    output:
    file("${oname}.graphml") into graphs

    """
    bamToEdges_mod2.py --sim --afmt bam --strong 150 --graphml "${oname}.graphml" --merged hic2ctg.bam hic2ctg.e hic2ctg.n
    """
}

/**
 * Map WGS reads to contigs to infer read-depth
 */

a = wgs_contigs.onCopy()
        .map { f -> [f.name[0..-15], f] }

wgsmap_sweep = wgs_reads.onCopy()
        .map { t -> [t[2], t[0], t[1]] }
        .phase(a){ t -> t[0] }.map { t -> [*(t[0]), t[1][1]] }

process WGSMap {
    cache 'deep'
    publishDir out_path, mode: 'copy', overwrite: 'false'

    input:
    set oname, file('r1.fq'), file('r2.fq'), file('contigs.fa') from wgsmap_sweep

    output:
    set file("${oname}.wgs2ctg.bam"), oname into wgs2ctg_mapping

    """
    bwa index contigs.fa
    bwa mem -t 1 contigs.fa r1.fq r2.fq | samtools view -bS - | samtools sort - "${oname}.wgs2ctg"
    """
}

/**
 * Infer contig coverage from wgs mapping
 */

wgs2ctg_mapping = ChannelDuplicator.createFrom(wgs2ctg_mapping)

cov_sweep = wgs2ctg_mapping.onCopy()

process InferReadDepth {
    cache 'deep'
    publishDir out_path, mode: 'copy', overwrite: 'false'

    input:
    set file("wgs2ctg.bam"), oname from cov_sweep

    output:
    file("${oname}.wgs2ctg.cov") into wgs2ctg_coverage

    """
    bedtools genomecov -ibam wgs2ctg.bam | \
    awk '
    BEGIN{n=0}
    {
        # ignore whole genome records
        if (\$1 != "genome") {
            # store names as they appear in repository
            # we use this to preserve file order
            if (!(\$1 in seq_cov)) {
                name_repo[n++]=\$1
            }
            # sum uses relative weights from histogram
            seq_cov[\$1]+=\$2*\$3/\$4
        }
    }
    END{
        for (i=0; i<n; i++) {
            print i+1, name_repo[i], seq_cov[name_repo[i]]
        }
    }' > "${oname}.wgs2ctg.cov"
    """
}
