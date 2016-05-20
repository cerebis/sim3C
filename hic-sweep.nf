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

/**
 * Include and initialise some utility classes
 */
class Globals {
    static String separator = '%'
}
helper = this.class.classLoader.parseClass(new File('Helper.groovy')).newInstance()
duplicator = this.class.classLoader.parseClass(new File('ChannelDuplicator.groovy')).newInstance()

trees = file(params.trees).collectEntries{ [it.name, it] }
tables = file(params.tables).collectEntries{ [it.name, it] }
ancestor = file(params.ancestor)
donor = file(params.donor)

alpha_BL = helper.stringToList(params.alpha)
xfold = helper.stringToList(params.xfold)
nhic = helper.stringToList(params.hic_pairs)

/**
 * Generate simulated communities
 */

evo_sweep = Channel
        .from([ancestor])
        .spread([donor])
        .spread(alpha_BL)
        .spread(trees.values())
        .map { it += it.collect { helper.safeString(it) }.join(Globals.separator) }

process Evolve {
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'true'

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
descendents = duplicator.createFrom(descendents)

wgs_sweep = descendents.onCopy()
        .spread(tables.values())
        .spread(xfold)
        .map{ it += tuple(it[0].name[0..-8], it[1].name, it[2]).join(Globals.separator) }

process WGS_Reads {
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'true'

    input:
    set file('descendent.fa'), file('profile'), xf, oname from wgs_sweep

    output:
    set file("${oname}.wgs.r*.fq.gz"), oname into wgs_reads

    """
    metaART.py -C gzip -t $profile -z 1 -M $xf -S ${params.seed} -s ${params.wgs_ins_std} \
            -m ${params.wgs_ins_len} -l ${params.wgs_read_len} -n "${oname}.wgs" descendent.fa .
    """
}

/**
 * Generate 3C/HiC read-pairs
 */

hic_sweep = descendents.onCopy()
        .spread(tables.values())
        .spread(nhic)
        .map{ it += tuple(it[0].name[0..-8], it[1].name, it[2]).join(Globals.separator) }

process HIC_Reads {
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'true'

    input:
    set file('descendent.fa'), file('profile'), nh, oname from hic_sweep

    output:
    set file("${oname}.hic.fa.gz") into hic_reads

    """
    simForward.py -C gzip -r ${params.seed} -n $nh -l ${params.hic_read_len} -p ${params.hic_inter_prob} \
           -t $profile descendent.fa "${oname}.hic.fa.gz"
    """
}

/**
 * Assemble WGS read-pairs
 */

wgs_reads = duplicator.createFrom( wgs_reads.map { reads, oname -> [*reads, oname] } )

asm_sweep = wgs_reads.onCopy()

process Assemble {
    cpus 1
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'true'

    input:
    set file('wgs_R1.fq.gz'), file ('wgs_R2.fq.gz'), oname from asm_sweep

    output:
    file("${oname}.contigs.fasta") into wgs_contigs

    """
    a5_pipeline.pl --threads=1 --metagenome wgs_R1.fq.gz wgs_R2.fq.gz ${oname}
    """
}

/**
 * Generate contig->descendent gold standard mapping
 */

wgs_contigs = duplicator.createFrom(wgs_contigs)

a = wgs_contigs.onCopy().map { f -> [helper.removeLevels(f.name,2), f, f.name[0..-15]] }
tr_sweep = descendents.onCopy()
        .spread(tables.values())
        .map { t -> [t[0].name[0..-8], t[0]] }
        .phase(a){ t -> t[0] }.map { t -> [*t[0], *t[1][1..2]] }

process Truth {
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'true'

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
        .map { f -> [helper.removeLevels(f.name, 1), f] }
        .groupTuple()

hicmap_sweep = hic_reads
        .map { f -> [helper.removeLevels(f.name, 1), f] }
        .groupTuple()
        .cross(a)
        .map { t -> [t[0][0], helper.product(t[0][1], t[1][1])] }
        .flatMap { t -> t[1] } 
        .map { t -> [*t, (t[1].name[0..-15].split(Globals.separator) + t[0].name[0..-8].split(Globals.separator)[-1]).join(Globals.separator)  ] }

process HiCMap {
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'true'

    input:
    set file('hic.fa.gz'), file('contigs.fa'), oname from hicmap_sweep

    output:
    file("${oname}.hic2ctg*") into hic2ctg_mapping


    """
    bwa index contigs.fa
    bwa mem -t 1 contigs.fa hic.fa.gz | samtools view -bS - | samtools sort -l 9 - "${oname}.hic2ctg"
    samtools index "${oname}.hic2ctg.bam"
    samtools idxstats "${oname}.hic2ctg.bam" > "${oname}.hic2ctg.idxstats"
    samtools flagstat "${oname}.hic2ctg.bam" > "${oname}.hic2ctg.flagstat"
    """
}

/**
 * Generate contig graphs from HiC mappings
 */

hic2ctg_mapping = duplicator.createFrom(hic2ctg_mapping)

graph_sweep = hic2ctg_mapping.onCopy()
        .map {t -> [t[0], t[1], t[0].name[0..-13] ] }

process Graph {
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'true'

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
    publishDir params.output, mode: 'symlink', overwrite: 'true'

    input:
    set oname, file('r1.fq.gz'), file('r2.fq.gz'), file('contigs.fa') from wgsmap_sweep

    output:
    set file("${oname}.wgs2ctg.bam"), oname into wgs2ctg_mapping

    """
    bwa index contigs.fa
    bwa mem -t 1 contigs.fa r1.fq.gz r2.fq.gz | samtools view -bS - | samtools sort -l 9 - "${oname}.wgs2ctg"
    """
}

/**
 * Infer contig coverage from wgs mapping
 */

wgs2ctg_mapping = duplicator.createFrom(wgs2ctg_mapping)

cov_sweep = wgs2ctg_mapping.onCopy()

process InferReadDepth {
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'true'

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
            # sum uses relative weights from histogram
            }
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
