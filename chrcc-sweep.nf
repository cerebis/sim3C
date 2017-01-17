#!/usr/bin/env nextflow
/**
 * Chromosome conforation capture (HiC/3C) data generation workflow
 *
 * Using a parametric sweep defined in chrcc.yaml, generate a set of simulated sequencing data.
 *
 * By defining the variable values in chrcc.yaml, a sweep will be performed as a permutation of
 * all variables. The results will be copied to the defined output folder (default: out).
 *
 * This folder is then used as input for subsequent workflows: chrcc-cluster.nf and chrcc-aggregate.nf
 *
 * Please note that the depth of the parametric sweep is such that choosing a wide range of values
 * at multiple levels can result in a very large final set and long computation time.
 *
 * Usage: chrcc-sweep.nf [--debug]
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

MetaSweeper ms = MetaSweeper.fromFile(new File('chrcc.yaml'))

//
// Generate phylogenetic trees for each clade within each community
//

def sweep = MetaSweeper.createSweep()
                .withVariable('seed', ms.variables.seed)
                .withVariable('clade', ms.variables.community.clades)
                .describe('Tree Generation')

// channel composed of the permutation of variables
gen_in = sweep.permutedChannel()

process TreeGen {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, seed, clade from gen_in

    output:
    set key, file("${key}.nwk"), seed, clade into tree_out

    script:
    if (params.debug) {
        """
        echo "$key ${clade.tree}" > "${key}.nwk"
        """
    }
    else {
        if (clade.isDefined()) {
            """
            echo "${clade.getDefined()}" | tree_scaler.py --max-height 0.1 - ${key}.nwk
            """
        }
        else {
            // presently we always assume birth-death
            assert clade.isSupportedAlgorithm() : 'Only birth_death is currently supported'
            """
            tree_generator.py --seed $seed --prefix ${clade.prefix} --suppress-rooting --mode random \
                --max-height 0.1 --birth-rate ${clade.tree.birth_rate} --death-rate ${clade.tree.death_rate} \
                --format newick --num-taxa ${clade.ntaxa} ${key}.nwk
            """
        }
    }
}


//
// Generate evolved sequences for each clade from each community
//

(tree_out, evo_in) = tree_out.into(2)

// add variation on alpha
sweep.withVariable('alpha', ms.variables.alpha)
        .describe('Evolve Clades')

// extend the channel to include new parameter
evo_in = sweep.extendChannel(evo_in, 'alpha')

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


//
// Merge evolved sequences from the clades into whole communities
//

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

//
// Generate abundance profiles for each clade within each community
//
(merge_seq_out, prof_in) = merge_seq_out.into(2)

community = Channel.value(ms.variables['community'])

process ProfileGen {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('community.fa') from prof_in
    val community

    output:
    set key, file("${key}.prf") into prof_out

    script:
    def mu = community.profile.mu
    def sigma = community.profile.sigma
    if (params.debug) {
        """
        echo "$key $mu $sigma" > "${key}.prf"
        """
    }
    else {
        """
        profile_generator.py --seed ${key['seed']} --dist lognormal --lognorm-mu $mu \
            --lognorm-sigma $sigma community.fa ${key}.prf
        """
    }
}


//
// Merge and pair community sequences and profiles
//

(merge_seq_out, seq_prof) = merge_seq_out.into(2)
(merge_prof_out, tmp) = prof_out.into(2)

// select just the community sequences
seq_prof = seq_prof.map { it.pick(1) }
        // combine with their respective profiles, then flatten and simplify the rows
        .phase(tmp).map { it = it.flatten(); it.pick(1, 3) }

//
// Generate shotgun sequencing reads for for each whole community
//
(seq_prof, wgs_in) = seq_prof.into(2)

// Add wgs coverage to sweep
sweep.withVariable('xfold', ms.variables.xfold)
        .describe('WGS Read Generation')

// extend the channel
wgs_in = sweep.extendChannel(wgs_in, 'xfold')

process WGS_Reads {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(comm_seq), file(comm_prof), xfold from wgs_in

    output:
    set key, file("${key}.wgs.r1.fq.gz"), file("${key}.wgs.r2.fq.gz"), file(comm_seq) into wgs_out

    script:
    opts = ms.options['wgs']
    if (params.debug) {
        """
        echo "metaART.py -C gzip --profile $comm_prof -z 1 -M $xfold -S ${key['seed']} \
                -s ${opts['insert_sd']} -m ${opts['insert_mean']} \
                -l ${opts['read_len']} -n ${key}.wgs $comm_seq ." > ${key}.wgs.r1.fq.gz

        echo "metaART.py -C gzip --profile $comm_prof -z 1 -M $xfold -S ${key['seed']} \
                -s ${opts['insert_sd']} -m ${opts['insert_mean']} \
                -l ${opts['read_len']} -n ${key}.wgs $comm_seq ." > ${key}.wgs.r2.fq.gz
        """
    }
    else {
        """
        export PATH=\$EXT_BIN/art:\$PATH
        metaART.py -C gzip --profile $comm_prof -z 1 -M $xfold -S ${key['seed']} \
                -s ${opts['insert_sd']} -m ${opts['insert_mean']} \
                -l ${opts['read_len']} -n "${key}.wgs" $comm_seq .
        wait_on_openfile.sh ${key}.wgs.r1.fq.gz
        wait_on_openfile.sh ${key}.wgs.r2.fq.gz
        """
    }
}


//
// Generate Conformation Capture (HiC/3C) reads for each whole community
//
(seq_prof, ccc_in) = seq_prof.into(2)

// Add 3C coverage to sweep
sweep.withVariable('num_3c', ms.variables.num_3c)
        .describe('Chromosome Conformation Capture (HiC/3C) Read Generation')

// extend the channel
ccc_in = sweep.extendChannel(ccc_in, 'num_3c')


process CCC_Reads {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(comm_seq), file(comm_prof), num_3c from ccc_in

    output:
    set key, file("${key}.ccc.fa.gz") into ccc_out

    script:
    opts = ms.options['ccc']
    if (params.debug) {
        """
        echo "sim3C.py -C gzip -m ${opts['method']} -r ${key['seed']} -n $num_3c -l ${opts['read_len']} -e ${opts['enzyme']} \
            --insert-mean ${opts['insert_mean']} --insert-sd ${opts['insert_sd']} --insert-max ${opts['insert_max']} \
            --machine-profile ${opts['machine_profile']} --profile $comm_prof $comm_seq \
            ${key}.ccc.fa.gz" > ${key}.ccc.fa.gz
        """
    }
    else {
        """
        sim3C.py -C gzip -m ${opts['method']} -r ${key['seed']} -n $num_3c -l ${opts['read_len']} -e ${opts['enzyme']} \
            --insert-mean ${opts['insert_mean']} --insert-sd ${opts['insert_sd']} --insert-max ${opts['insert_max']} \
            --machine-profile ${opts['machine_profile']} --profile $comm_prof $comm_seq ${key}.ccc.fa.gz

        wait_on_openfile.sh ${key}.ccc.fa.gz
        """
    }
}

//
// Assemble WGS reads
//
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

//
// Infer Truth Tables for each community by mapping contigs to community references
//
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

//
// Map CCC reads to assembled contigs
//
(asm_out, cccmap_in) = asm_out.into(2)

// join 3C reads and the results of assembly
cccmap_in = sweep.joinChannels(cccmap_in, ccc_out, 2)
        // select just contigs and 3C reads
        .map{ it.pick(1, -1) }

process CCCMap {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(contigs), file(ccc_reads) from cccmap_in

    output:
    set key, file("${key}.ccc2ctg.bam"), file(ccc_reads), file(contigs) into cccmap_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.ccc2ctg.bam
        """
    }
    else {
        """
        export PATH=\$EXT_BIN/a5/bin:\$PATH
        if [ ! -e "${contigs}.bwt" ]
        then
            bwa index $contigs
        fi
        bwa mem -t 1 $contigs $ccc_reads | samtools view -bS - | samtools sort -l 9 - "${key}.ccc2ctg"
        samtools index "${key}.ccc2ctg.bam"
        samtools idxstats "${key}.ccc2ctg.bam" > "${key}.ccc2ctg.idxstats"
        samtools flagstat "${key}.ccc2ctg.bam" > "${key}.ccc2ctg.flagstat"
        """
    }
}

//
// Generate contig graphs
//
(cccmap_out, graph_in) = cccmap_out.into(2)

// select just ChrCC bam, reads and contigs
graph_in = graph_in.map{ it.pick(1,2,3) }

process Graph {

    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(ccc2ctg), file(ccc_reads), file(contigs) from graph_in

    output:
    set key, file("${key}.graphml"), file(ccc_reads), file(contigs) into graph_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.graphml
        """
    }
    else {
        """
        if [ ! -e "${ccc2ctg}.bai" ]
        then
            samtools index $ccc2ctg
        fi
        bamToEdges.py --strong ${ms.options['ccc']['read_len']} --preserve-zerodeg --merged $ccc2ctg -o ${key}
        """
    }
}


//
// Map WGS reads to contigs
//
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


//
// Calculate assembly contig coverage
//
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
