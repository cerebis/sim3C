/*
 * meta-sweeper - for performing parametric sweeps of simulated
 * metagenomic sequencing experiments.
 * Copyright (C) 2016 "Aaron E Darling"
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
    publishDir params.output, mode: 'symlink', overwrite: 'false'

    input:
    set file('ancestral.fa'), file('donor.fa'), alpha, file('input_tree'), oname from evo_sweep

    output:
    set file("${oname}.evo.fa") into descendents
    val oname into evolvedoname

    """
    scale_tree.py -a ${alpha} input_tree scaled_tree
    sgEvolver --indel-freq=${params.indel_freq} --small-ht-freq=${params.small_ht_freq} --large-ht-freq=${params.large_ht_freq} \
         --inversion-freq=${params.inversion_freq} --random-seed=${params.seed} scaled_tree \
         ancestral.fa donor.fa "${oname}.evo.aln" "${oname}.evo.fa"
    strip_semis.sh "${oname}.evo.fa"
    """
}

ancestor_channel = duplicator.createFrom(Channel.from([ancestor]))

/**
 * Generate WGS read-pairs
 */
wgs_descendents = duplicator.createFrom(descendents)

wgs_sweep = wgs_descendents.onCopy()
        .spread(tables.values())
        .spread(xfold)
        .map{ it += tuple(it[0].name[0..-8], it[1].name, it[2]).join(Globals.separator) }

process WGS_Reads {
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'false'

    input:
    set file('descendent.fa'), file('profile'), xf, oname from wgs_sweep

    output:
    set file("${oname}.wgs.*.r1.fq.gz"), file("${oname}.wgs.*.r2.fq.gz"), oname into wgs_reads

    """
    metaART.py -C gzip -t $profile -M $xf -S ${params.seed} -z ${params.num_samples} -s ${params.wgs_ins_std} \
            -m ${params.wgs_ins_len} -l ${params.wgs_read_len} -n "${oname}.wgs" descendent.fa .
    """
}


/**
 * Map WGS read-pairs to reference
 */
map_ancestor = ancestor_channel.onCopy()
process ReadMap {
    cpus 1
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'false'

    input:
    file(ancestor) from map_ancestor
    set file(r1file), file(r2file), oname from wgs_reads

    output:
    set file("*.bam"), oname into map_bams

    """
    cp -L ancestor.fa ancestral.fa
    bwa index ancestral.fa
    for rr in `ls *.r1.fq.gz`
    do
        rbase=`basename \$rr .r1.fq.gz`
        r2=\$rbase.r2.fq.gz
        bwa mem -t 1 ancestral.fa \$rr \$r2 | samtools view -bS - | samtools sort -l 9 - \$rbase
    done
    """
}

/**
 * Deconvolve the SNVs into strain genotypes
 */
decon_ancestor = ancestor_channel.onCopy()
process Deconvolve {
    cpus 1
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'false'

    input:
    file(ancestor) from decon_ancestor
    set file('*.bam'), oname from map_bams

    output:
    set file("${oname}.decon.csv"), file("${oname}.snv_file.data.R"), file("${oname}.strains.tre"), oname into deconvolution

    """
    snvbpnmft.py . 4 ancestor.fa *.bam
    #java -Xmx1000m -jar \$JARPATH/beast.jar beast.xml 
    #java -jar \$JARPATH/treeannotator.jar -burnin 1000 -heights mean aln.trees strains.tre
    touch strains.tre decon.csv
    """
}


/**
 * Record the true strain genotypes
 */
truth_ancestor = ancestor_channel.onCopy()
process Truth {
    cache 'deep'
    publishDir params.output, mode: 'symlink', overwrite: 'false'

    input:
    file(ancestor) from truth_ancestor
    set file('evo.fa') from descendents
    val oname from evolvedoname

    output:
    set file("${oname}.truth.tsv") into truth

    """
    strain_truth.py --mauve-path=\$MAUVEPATH -o $oname.truth.tsv evo.fa ancestor.fa
    """
}

