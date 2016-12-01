#!/usr/bin/env nextflow
/**
 * HiC data generation workflow
 *
 * Using a parametric sweep defined in hic.yaml, generate a set of simulated sequencing data.
 *
 * By defining the variable values in hic.yaml, a sweep will be performed as a permutation of
 * all variables. The results will be copied to the defined output folder (default: out).
 *
 * This folder is then used as input for subsequent workflows: hic-cluster.nf and hic-aggregate.nf
 *
 * Please note that the depth of the parametric sweep is such that choosing a wide range of values
 * at multiple levels can result in a very large final set and long computation time.
 *
 * Usage: hic-sweep.nf [--debug]
 */
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

import MetaSweeper

MetaSweeper ms = MetaSweeper.fromFile(new File('hic.yaml'))

/**
 * Generate phylogenetic trees for each clade within each community
 */

// Initial sweep begins with seeds and community's clades
gen_in = ms.createSweep()
        .withVariable('seed')
        .withVariable('community', true)
        .permute()

ms.describeSweep('Tree Generation')

process TreeGen {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, seed, clade from gen_in

    output:
    set key, file("${key}.nwk"), seed, clade into tree_out

    script:
    if (params.debug) {
        """
        echo $key > "${key}.nwk"
        """
    }
    else {
        """
        tree_generator.py --seed $seed --prefix ${clade.prefix} --suppress-rooting --mode random \
            --max-height 0.1 --birth-rate ${clade.tree.birth} --death-rate ${clade.tree.death} \
            --format newick --num-taxa ${clade.ntaxa} ${key}.nwk
        """
    }
}


/**
 * Generate evolved sequences for each clade from each community
 */

(tree_out, evo_in) = tree_out.into(2)

        // extend the sweep to include alpha
evo_in = ms.withVariable('alpha').extend(evo_in, 'alpha')

ms.describeSweep('Evolve Clades')

process Evolve {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, tree_file, seed, clade, alpha from evo_in

    output:
    set key, file("${key}.evo.fa"), seed, clade into evo_out

    script:
    if (params.debug) {
        """
        echo $key > "${key}.evo.fa"
        """
    }
    else {
        """
        scale_tree.py -a $alpha $tree_file scaled_tree
        \$EXT_BIN/sgevolver/sgEvolver --indel-freq=${ms.options.evo.indel_freq} \
            --small-ht-freq=${ms.options.evo.small_ht_freq} \
            --large-ht-freq=${ms.options.evo.large_ht_freq} \
            --inversion-freq=${ms.options.evo.inversion_freq} \
            --random-seed=$seed scaled_tree \
             $clade.ancestor $clade.donor "${key}.evo.aln" "${key}.evo.fa"
        strip_semis.sh "${key}.evo.fa"
        """
    }

}


/**
 * Generate abundance profiles for each clade within each community
 */

(evo_out, prof_in) = evo_out.into(2)

    // just the clade sequences, which are used to obtain taxon ids.
prof_in = prof_in.map { it.pick(1) }

process ProfileGen {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('clade_seq') from prof_in

    output:
    set key, file("${key}.prf") into prof_out

    script:
    if (params.debug) {
        """
        echo $key > "${key}.prf"
        """
    }
    else {
        def mu = key['community'].profile.mu
        def sigma = key['community'].profile.sigma
        """
        profile_generator.py --seed ${key['seed']} --dist lognormal --lognorm-mu $mu \
            --lognorm-sigma $sigma clade_seq ${key}.prf
        """
    }
}


/**
 * Merge clade abundance profiles into whole community profiles
 */
(prof_out, merge_prof_in) = prof_out.into(2)

        // group by a reduced key that is only the random seed and alpha
merge_prof_in = merge_prof_in.groupBy{ it[0].selectedKey('seed', 'alpha') }
        // convert the resulting map of sweep point results into table format and sort by file name
        .flatMap { it.collect { k, v -> [k, v.collect { vi -> vi[1] }.toSorted { a, b -> a.name <=> b.name }] } }

process ProfileMerge {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('clade_profile') from merge_prof_in

    output:
    set key, file("${key}.mprf") into merge_prof_out

    script:
    if (params.debug) {
        """
        echo $key > "${key}.mprf"
        """
    } else {
        """
        profile_merge.py clade_profile* ${key}.mprf
        """
    }
}

/**
 * Merge evolved sequences from the clades into whole communities
 */
(evo_out, merge_seq_in) = evo_out.into(2)

        // group by a reduced key that is only the random seed and alpha
merge_seq_in = merge_seq_in.groupBy { it.getKey().selectedKey('seed', 'alpha') }
        // convert the resulting map of sweep point results into table format and sort by file name
        .flatMap { it.collect { k, v -> [k, v.collect { vi -> vi[1] }.toSorted { a, b -> a.name <=> b.name }] } }

process MergeClades {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('clade_seq') from merge_seq_in

    output:
    set key, file("${key}.community.fa") into merge_seq_out

    script:
    if (params.debug) {
        """
        echo $key > "${key}.community.fa"
        """
    }
    else {
        """
        cat clade_seq* >> ${key}.community.fa
        """
    }
}

/**
 * Prepare a channel which is composed of paired community sequences and profiles
 */

(merge_seq_out, seq_prof) = merge_seq_out.into(2)
(merge_prof_out, tmp) = merge_prof_out.into(2)

        // select just the community sequences
seq_prof = seq_prof.map { it.pick(1) }
        // combine with their respective profiles, then flatten and simplify the rows
        .phase(tmp).map { it = it.flatten(); it.pick(1, 3) }


/**
 * Generate shotgun sequencing reads for for each whole community
 */
(seq_prof, wgs_in) = seq_prof.into(2)

        // Add WGS coverage to the sweep
wgs_in = ms.withVariable('xfold')
        // extend the channel
        .extend(wgs_in, 'xfold')

ms.describeSweep('WGS Read Generation')

process WGS_Reads {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(comm_seq), file(comm_prof), xfold from wgs_in

    output:
    set key, file("${key}.wgs.r1.fq.gz"), file("${key}.wgs.r2.fq.gz"), file(comm_seq) into wgs_out

    script:
    if (params.debug) {
        """
        echo "metaART.py -C gzip --profile $comm_prof -z 1 -M $xfold -S ${key['seed']} \
                -s ${ms.options['wgs']['ins_std']} -m ${ms.options['wgs']['ins_len']} \
                -l ${ms.options['wgs']['read_len']} -n ${key}.wgs $comm_seq ." > ${key}.wgs.r1.fq.gz

        echo "metaART.py -C gzip --profile $comm_prof -z 1 -M $xfold -S ${key['seed']} \
                -s ${ms.options['wgs']['ins_std']} -m ${ms.options['wgs']['ins_len']} \
                -l ${ms.options['wgs']['read_len']} -n ${key}.wgs $comm_seq ." > ${key}.wgs.r2.fq.gz
        """
    }
    else {
        """
        export PATH=\$EXT_BIN/art:\$PATH
        metaART.py -C gzip --profile $comm_prof -z 1 -M $xfold -S ${key['seed']} \
                -s ${ms.options['wgs']['ins_std']} -m ${ms.options['wgs']['ins_len']} \
                -l ${ms.options['wgs']['read_len']} -n "${key}.wgs" $comm_seq .
        wait_on_openfile.sh ${key}.wgs.r1.fq.gz
        wait_on_openfile.sh ${key}.wgs.r2.fq.gz
        """
    }
}


/**
 * Generate 3C reads for each whole community
 */
(seq_prof, hic_in) = seq_prof.into(2)

        // Add 3C coverage to the sweep
hic_in = ms.withVariable('n3c')
        // extend the channel
        .extend(hic_in, 'n3c')

ms.describeSweep('HiC Read Generation')

process HIC_Reads {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(comm_seq), file(comm_prof), n3c from hic_in

    output:
    set key, file("${key}.hic.fa.gz") into hic_out

    script:
    hic_sites = ms.options.n3c.hic_sites ? '--site-dup' : ''
    if (params.debug) {
        """
        echo "sim3C.py $hic_sites -C gzip -r ${key['seed']} -n $n3c -l ${ms.options['n3c']['read_len']} \
            --inter-prob ${ms.options['n3c']['inter_prob']} --profile $comm_prof $comm_seq \
            ${key}.hic.fa.gz" > ${key}.hic.fa.gz
        """
    }
    else {
        """
        sim3C.py $hic_sites -C gzip -r ${key['seed']} -n $n3c -l ${ms.options['n3c']['read_len']} \
            --inter-prob ${ms.options['n3c']['inter_prob']} --profile $comm_prof $comm_seq ${key}.hic.fa.gz
        wait_on_openfile.sh ${key}.hic.fa.gz
        """
    }
}


/**
 * Assemble WGS reads
 */
(wgs_out, asm_in) = wgs_out.into(2)

        // select just read-set files (R1, R2) and the community sequences
asm_in = asm_in.map { it.pick(1, 2, 3) }

process Assemble {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(reads1), file(reads2), file(comm_seq) from asm_in

    output:
    set key, file("${key}.contigs.fasta"), file(reads1), file(reads2), file(comm_seq) into asm_out

    script:
    if (params.debug) {
        """
        echo "\$EXT_BIN/a5/bin/a5_pipeline.pl --threads=1 --metagenome $reads1 $reads2 $key" > ${key}.contigs.fasta
        """
    }
    else {
        """
        export PATH=\$EXT_BIN/a5/bin:\$PATH
        a5_pipeline.pl --threads=1 --metagenome $reads1 $reads2 $key
        bwa index ${key}.contigs.fasta
        """
    }
}

/**
 * Infer Truth Tables for each community by mapping contigs to community references
 */
(asm_out, truth_in) = asm_out.into(2)

        // select just contigs and community sequences
truth_in = truth_in.map{ it.pick(1, 4) }

process Truth {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(contigs), file(comm_seq) from truth_in

    output:
    set key, file("${key}.truth"), file(contigs), file(comm_seq) into truth_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.truth
        """
    }
    else {
        """
        export PATH=\$EXT_BIN/last:\$PATH
        if [ ! -e db.prj ]
        then
            lastdb db $comm_seq
        fi

        lastal -P 1 db $contigs | maf-convert psl > ctg2ref.psl
        alignmentToTruth.py --ofmt json ctg2ref.psl "${key}.truth"
        """
    }

}

/**
 * Map HiC reads to assembled contigs
 */
(asm_out, hicmap_in) = asm_out.into(2)

        // join 3C reads and the results of assembly
hicmap_in = ms.joinChannels(hicmap_in, hic_out, 2)
        // select just contigs and 3C reads
        .map{ it.pick(1, -1) }

process HiCMap {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(contigs), file(hic_reads) from hicmap_in

    output:
    set key, file("${key}.hic2ctg.bam"), file(hic_reads), file(contigs) into hicmap_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.hic2ctg.bam
        """
    }
    else {
        """
        export PATH=\$EXT_BIN/a5/bin:\$PATH
        if [ ! -e "${contigs}.bwt" ]
        then
            bwa index $contigs
        fi
        bwa mem -t 1 $contigs $hic_reads | samtools view -bS - | samtools sort -l 9 - "${key}.hic2ctg"
        samtools index "${key}.hic2ctg.bam"
        samtools idxstats "${key}.hic2ctg.bam" > "${key}.hic2ctg.idxstats"
        samtools flagstat "${key}.hic2ctg.bam" > "${key}.hic2ctg.flagstat"
        """
    }
}

/**
 * Generate contig graphs
 */
(hicmap_out, graph_in) = hicmap_out.into(2)

        // select just hic bam, hic reads and contigs
graph_in = graph_in.map{ it.pick(1,2,3) }

process Graph {

    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(hic2ctg), file(hic_reads), file(contigs) from graph_in

    output:
    set key, file("${key}.graphml"), file(hic_reads), file(contigs) into graph_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.graphml
        """
    }
    else {
        """
        if [ ! -e "${hic2ctg}.bai" ]
        then
            samtools index $hic2ctg
        fi
        bamToEdges.py --strong 150 --preserve-zerodeg --merged $hic2ctg -o ${key}
        """
    }
}


/**
 * Map WGS reads to contigs
 */
(asm_out, wgsmap_in) = asm_out.into(2)

        // select just contigs, r1 and r2
wgsmap_in = wgsmap_in.map{ it.pick(1, 2, 3) }

process WGSMap {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(contigs), file(reads1), file(reads2) from wgsmap_in

    output:
    set key, file("${key}.wgs2ctg.bam"), file(contigs) into wgsmap_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.wgs2ctg.bam
        """
    }
    else {
        """
        export PATH=\$EXT_BIN/a5/bin:\$PATH
        if [ ! -e "${contigs}.bwt" ]
        then
            bwa index $contigs
        fi
        bwa mem -t 1 $contigs $reads1 $reads2 | samtools view -bS - | samtools sort -l 9 - "${key}.wgs2ctg"
        """
    }
}


/**
 * Calculate assembly contig coverage
 */
(wgsmap_out, cov_in) = wgsmap_out.into(2)

        // select wgs bam and contigs
cov_in = cov_in.map{ it.pick(1, 2) }

process InferReadDepth {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(wgs2ctg), file(contigs) from cov_in

    output:
    set key, file("${key}.wgs2ctg.cov"), file(wgs2ctg), file(contigs) into cov_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.wgs2ctg.cov
        """
    }
    else {
        """
        \$EXT_BIN/bedtools/bedtools genomecov -ibam $wgs2ctg | \
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
        """
    }
}
