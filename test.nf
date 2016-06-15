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

import static Helper.*


/*
 * Basic sweep collections
 *
 * These are set in the configuration file or on the command line.
 *
 * Explicit series must be quoted and delimited by commas, whereas wildcards
 * can be passed to filesystem related variables.
 *
 * Eg. --trees /path/to/trees/*.nwk --alpha "1,2,3,4"
 */

sweep = new Sweep()

sweep['ancestor'] = files(params.ancestor)
sweep['donor'] = files(params.donor)
sweep['alpha'] = stringToList(params.alpha)
sweep['tree'] = absPath(params.trees)
sweep['profile'] = absPath(params.profiles)
sweep['xfold'] = stringToList(params.xfold).collect{(int)it * 99}
sweep['n3c'] = stringToList(params.hic_pairs)


next = Channel.from(
            sweep.permuteNames('ancestor', 'donor', 'alpha', 'tree'))
        .map{ it.addKey() }

/**
 * Generate reference genomes from ancestor
 */

process Evolve {
    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    set key, ancestor, donor, alpha, tree from next

    output:
    set key, file("${key}.evo.fa") into evo_out

    """
    scale_tree.py -a $alpha $tree scaled_tree
    \$EXT_BIN/sgevolver/sgEvolver --indel-freq=${params.indel_freq} --small-ht-freq=${params.small_ht_freq} \
        --large-ht-freq=${params.large_ht_freq} --inversion-freq=${params.inversion_freq} \
        --random-seed=${params.seed} scaled_tree \
         $ancestor $donor "${key}.evo.aln" "${key}.evo.fa"
    strip_semis.sh "${key}.evo.fa"
    """

    /*"""
    echo "test" > "${key}.evo.fa"
    """*/
}

// add a name to new output
evo_out = evo_out.map { it.nameify(1, 'ref_seq') }

/**
 *  Make WGS reads
 */
(evo_out, next) = evo_out.into(2)
next = sweep.extend(next, 'profile', 'xfold')

process WGS_Reads {
    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    set key, ref_seq from next

    output:
    set key, file("${key}.wgs.r*.fq.gz"), ref_seq into wgs_out

    /*"""
    export PATH=\$EXT_BIN/art:\$PATH
    metaART.py -C gzip -t $profile -z 1 -M $xfold -S ${params.seed} -s ${params.wgs_ins_std} \
            -m ${params.wgs_ins_len} -l ${params.wgs_read_len} -n "${key}.wgs" $ref_seq .
    """*/

    """
    echo $key > ${key}.wgs.r1.fq.gz
    echo $key > ${key}.wgs.r2.fq.gz
    """
}

// add a name to new output
wgs_out = wgs_out.map { it.nameify(1, 'wgs_reads') }

/**
 * Make HiC reads
 */
(evo_out, next) = evo_out.into(2)
next = sweep.extend(next, 'profile', 'n3c')

process HIC_Reads {
    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    set key, ref_seq from next

    output:
    set key, file("${key}.hic.fa.gz"), ref_seq into hic_out

    /*"""
    simForward.py -C gzip -r ${params.seed} -n $n3c -l ${params.hic_read_len} -p ${params.hic_inter_prob} \
           -t $profile $ref_seq "${key}.hic.fa.gz"
    """*/

    """
    echo "${key['profile']}, ${key['n3c']}" > ${key}.hic.fa.gz
    """
}

// add a name to new output
hic_out = hic_out.map { it.nameify(1, 'hic_reads') }


/**
 * Assemble WGS reads
 */
(wgs_out, next) = wgs_out.into(2)

process Assemble {
    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    set key, reads, ref_seq from next

    output:
    set key, file("${key}.contigs.fasta"), reads, ref_seq into asm_out

    /*"""
    \$EXT_BIN/a5/bin/a5_pipeline.pl --threads=1 --metagenome $R1 $R2 ${key}
    """*/

    """
    echo $key $reads > ${key}.contigs.fasta
    """
}

// add a name to new output
asm_out = asm_out.map { it.nameify(1, 'contigs') }


/**
 * Make Truth Tables
 */
(asm_out, next) = asm_out.into(2)

process Truth {
    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    set key, contigs, reads, ref_seq from next

    output:
    set key, file("${key}.truth"), contigs, reads, ref_seq into truth_out

    /*"""
    if [ ! -e db.prj ]
    then
        lastdb db $ref_seq
    fi

    \$EXT_BIN/last/lastal -P 1 db $contigs | maf-convert psl > ctg2ref.psl
    alignmentToTruth.py --ofmt json ctg2ref.psl "${key}.truth"
    """*/

    """
    echo $key > ${key}.truth
    """
}

// add a name to new output
truth_out = truth_out.map { it.nameify(1, 'truth') }

/**
 * Map HiC reads to assembled contigs
 */
(asm_out, next) = asm_out.into(2)

next = hic_out.map { it -> [ (it[0]).id(5), *it[1..-1] ] }
        .cross( next.map { it -> [ (it[0]).id(5), *it[1..-1] ] } )
        .map { it.flatten().unique() }
        .map { it.pick('contigs', 'hic_reads') }


process HiCMap {
    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    set key, contigs, hic_reads from next

    output:
    set key, file("${key}.hic2ctg.bam"), contigs, hic_reads into hicmap_out

    /*"""
    export PATH=\$EXT_BIN/a5/bin:\$PATH
    bwa index $contigs
    bwa mem -t 1 $contigs $hic_reads | samtools view -bS - | samtools sort -l 9 - "${key}.hic2ctg"
    samtools index "${key}.hic2ctg.bam"
    samtools idxstats "${key}.hic2ctg.bam" > "${key}.hic2ctg.idxstats"
    samtools flagstat "${key}.hic2ctg.bam" > "${key}.hic2ctg.flagstat"
    """*/

    """
    echo $key > ${key}.hic2ctg.bam
    """
}

// add a name to new output
hicmap_out = hicmap_out.map { it.nameify(1, 'hic2ctg') }

/**
 * Generate contig graphs
 */
(hicmap_out, next) = hicmap_out.into(2)

process Graph {

    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    set key, hic2ctg, contigs, hic_reads from next

    output:
    set key, file("${key}.graphml"), contigs, hic_reads into graph_out

    /*"""
    bamToEdges_mod2.py --sim --afmt bam --strong 150 --graphml "${key}.graphml" --merged $hic_bam hic2ctg.e hic2ctg.n
    """*/

    """
    echo $key > ${key}.graphml
    """
}

// add a name to new output
graph_out = graph_out.map { it.nameify(1, 'graph') }

/**
 * Map WGS reads to contigs
 */

(asm_out, next) = asm_out.into(2)

process WGSMap {
    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    set key, contigs, reads, ref_seq from next

    output:
    set key, file("${key}.wgs2ctg.bam"), contigs, reads, ref_seq into wgsmap_out

    /*"""
    export PATH=\$EXT_BIN/a5/bin:\$PATH
    bwa index $contigs
    bwa mem -t 1 $contigs $R1 $R2 | samtools view -bS - | samtools sort -l 9 - "${key}.wgs2ctg"
    """*/

    """
    echo $key > ${key}.wgs2ctg.bam
    """
}

// add a name to new output
wgsmap_out = wgsmap_out.map { it.nameify(1, 'wgs2ctg') }

/**
 * Calculate assembly contig coverage
 */
(wgsmap_out, next) = wgsmap_out.into(2)


process InferReadDepth {
    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    set key, wgs2ctg, contigs, reads, ref_seq from next

    output:
    set key, file("${key}.wgs2ctg.cov"), wgs2ctg, contigs, reads, ref_seq into cov_out

    /*"""
    \$EXT_BIN/bedtools/bedtools genomecov -ibam $wgs_bam | \
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
    }' > "${key}.wgs2ctg.cov"
    """*/

    """
    echo $key > ${key}.wgs2ctg.cov
    """
}

// add a name to new output
cov_out = cov_out.map { it.nameify(1, 'coverage') }
