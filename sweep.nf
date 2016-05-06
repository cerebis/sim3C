sweep_sep = '%'
alpha_BL = [1]
xfold = [1]
nhic = [10000]
trees = file('nf/ref_data/trees/*nwk').collectEntries{ [it.name, it] }
tables = file('nf/ref_data/tables/*table').collectEntries{ [it.name, it] }
ancestor = file('nf/ref_data/NC*raw').collectEntries{ [it.name, it] }
out_path = 'junk-out'


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

String dropSuffix(str) {
    return str.lastIndexOf('.').with {it != -1 ? str[0..<it] : str}
}

String safeString(val) {
    String s = ""
    if (val instanceof Path) {
        s = val.name //dropSuffix(val.name)
    }
    else {
        s = val.toString()
    }
    return s.replaceAll(/[\\\/]/, "_")
}

String[] splitSampleName(Path path) {
    m = (path.name =~ /^(.*)_?([rR][0-9]+).*$/)
    return m[0][1..2]
}
String removeLevels(Path path, int n) {
    name = path.toAbsolutePath().toString()
    return name.split('%')[0..-(n+1)].join('%')
}

String removeLevels(String name, int n) {
    return name.split('%')[0..-(n+1)].join('%')
}

/**
 * Generate simulated communities
 */

sweep = Channel
        .from(ancestor.values())
        .spread(alpha_BL)
        .spread(trees.values())
        .map { it += it.collect { safeString(it) }.join(sweep_sep) }

process Evolve {
    cache 'deep'
    publishDir out_path, mode: 'copy', overwrite: 'false'

    input:
    set file('ancestral.fa'), alpha, file('input_tree'), oname from sweep

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

(a, descendents) = descendents.into(2)
(wgs_sweep, hic_sweep, tr_sweep) = a.spread(tables.values()).into(3)

wgs_sweep = wgs_sweep
        .spread(xfold)
        .map{ it += tuple(it[0].name[0..-8], it[1].name, it[2]).join(sweep_sep) }

hic_sweep = hic_sweep
        .spread(nhic)
        .map{ it += tuple(it[0].name[0..-8], it[1].name, it[2]).join(sweep_sep) }

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

(asm_sweep, wgs_reads) = wgs_reads.map { reads, oname -> [*reads, oname] }.into(2)

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

(a, wgs_contigs) = wgs_contigs.into(2)
a = a.map { f -> [removeLevels(f.name,2), f, f.name[0..-15]] }
b = tr_sweep.map { t -> [t[0].name[0..-8], t[0]] }
tr_sweep = b.phase(a){ t -> t[0] }.map { t -> [*t[0], *t[1][1..2]] }

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

(a, wgs_contigs) = wgs_contigs.into(2)
a = a.map { f -> [removeLevels(f.name,1), f] }
b = hic_reads.map { f -> [removeLevels(f.name, 1), f, f.name[0..-8]] }
hicmap_sweep = b.phase(a){ t -> t[0] }.map { t -> [t[0][1], t[1][1], t[0][2]] }

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

(a, hic2ctg_mapping) = hic2ctg_mapping.into(2)
graph_sweep = a.map {t -> [t[0], t[1], t[0].name[0..-13] ] }

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

(a, wgs_contigs) = wgs_contigs.into(2)
a = a.map { f -> [f.name[0..-15], f] }
(b, wgs_reads) = wgs_reads.into(2)
wgsmap_sweep = b.map { t -> [t[2], t[0], t[1]] }.phase(a){ t -> t[0] }.map { t -> [*(t[0]), t[1][1]] }

process WGSMap {
    cache 'deep'
    publishDir out_path, mode: 'copy', overwrite: 'false'

    input:
    set oname, file('r1.fq'), file('r2.fq'), file('contigs.fa') from wgsmap_sweep

    output:
    set file("${oname}.wgs2ctg.bam"), file("${oname}.wgs2ctg.bam.bai"), oname into wgs2ctg_mapping

    """
    bwa index contigs.fa
    bwa mem -t 1 contigs.fa r1.fq r2.fq | samtools view -bS - | samtools sort - "${oname}.wgs2ctg"
    """
}

/**
 * Infer per-contig coverage from wgs mapping
 */

process InferReadDepth {
    cache 'deep'
    publishDir out_path, mode: 'copy', overwrite: 'false'

    input:
    set file("wgs2ctg.bam"), file('wgs2ctg.bam.bai'), oname from wgs2ctg_mapping

    output:
    file("${oname}.wgs2ctg.cov") into wgs2ctg_coverage

    """
    calc_coverage.sh "${oname}.wgs2ctg.bam" > "${oname}.wgs2ctg.cov"
    """

}