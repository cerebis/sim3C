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

MetaSweeper ms = MetaSweeper.fromFile(new File('sweep.yaml'))

/**
 * Tree Generation
 */

// Initial sweep begins with seeds and community's clades
gen_in = ms.createSweep()
        .withVariable('seed')
        .withVariable('community', true)
        .permute()

ms.describeSweep('Tree Generation')

// a newick tree is generated for each seed and clade def.
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
        tree_generator.py --seed $seed --prefix ${clade.value.prefix} --suppress-rooting --mode random \
            --max-height 0.1 --birth-rate ${clade.value.tree.birth} --death-rate ${clade.value.tree.death} \
            --format newick --num-taxa ${clade.value.ntaxa} ${key}.nwk
        """
    }
}


/**
 * Evolve Clade Sequences
 */

(tree_out, evo_in) = tree_out.into(2)
// add alpha to the sweep and extend the channel
evo_in = ms.withVariable('alpha').extend(evo_in, 'alpha')

ms.describeSweep('Evolve Clades')

// sequences are produced for all taxa in the clade
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
             $clade.value.ancestor $clade.value.donor "${key}.evo.aln" "${key}.evo.fa"
        strip_semis.sh "${key}.evo.fa"
        """
    }

}


/**
 * Profile Generation
 */

(evo_out, prof_in) = evo_out.into(2)
prof_in = prof_in.map{[it[0], it[1]]}

// generate an abundance profile for each clade
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
        def mu = key['community'].value.profile.mu
        def sigma = key['community'].value.profile.sigma
        """
        profile_generator.py --seed ${key['seed']} --dist lognormal --lognorm-mu $mu \
            --lognorm-sigma $sigma clade_seq ${key}.prf
        """
    }
}


/**
 * Profile Merging
 */
(prof_out, merge_prof_in) = prof_out.into(2)
merge_prof_in = merge_prof_in.groupBy{ it[0].selectedKey('seed', 'alpha') }
        .flatMap{ it.collect { k,v -> [k] +  v.flatten() } }
        .map{ [it[0], [it[2],it[4]]] }

// merge the clade profiles together
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
 * Merge clades into communities
 */
(evo_out, merge_seq_in) = evo_out.into(2)
merge_seq_in = merge_seq_in.groupBy{ it[0].selectedKey('seed', 'alpha') }
        .flatMap{it.collect { k,v -> [k] +  v.collect{cl -> [cl[1],cl[3]]}.flatten() } }
        .map{ [it[0], [it[1],it[3]]] }

// the sequences for all clades are concatenated together
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
 * Make WGS reads
 */

// join the results of both sequence and profile generation
// pick out just the required variables.
(merge_seq_out, wgs_in) = merge_seq_out.into(2)
(merge_prof_out, tmp) = merge_prof_out.into(2)
wgs_in = wgs_in.map{ [it[0], it[1]] }.phase(tmp).map{ [*it[0], it[1][1]]}
// add WGS coverage to the sweep, extend our channel
wgs_in = ms.withVariable('xfold').extend(wgs_in, 'xfold')
ms.describeSweep('WGS Read Generation')

// from community sequences, generate simulated WGS reads
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
 * Make 3C reads
 */
// join the results of both sequence and profile generation
// pick out just the required variables.
(merge_seq_out, hic_in) = merge_seq_out.into(2)
(merge_prof_out, tmp) = merge_prof_out.into(2)
hic_in = hic_in.map{ [it[0], it[1]]}.phase(tmp).map{ [*it[0], it[1][1]]}
// add 3C coverage to the sweep, extend our channel
hic_in = ms.withVariable('n3c').extend(hic_in, 'n3c')
ms.describeSweep('HiC Read Generation')

// from community sequences, generate 3C sequence
process HIC_Reads {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(comm_seq), file(comm_prof), n3c from hic_in

    output:
    set key, file("${key}.hic.fa.gz") into hic_out

    script:
    if (params.debug) {
        """
        echo "sim3C.py -C gzip -r ${key['seed']} -n $n3c -l ${ms.options['n3c']['read_len']} \
            --inter-prob ${ms.options['n3c']['inter_prob']} --profile $comm_prof $comm_seq \
            ${key}.hic.fa.gz" > ${key}.hic.fa.gz
        """
    }
    else {
        """
        sim3C.py -C gzip -r ${key['seed']} -n $n3c -l ${ms.options['n3c']['read_len']} \
            --inter-prob ${ms.options['n3c']['inter_prob']} --profile $comm_prof $comm_seq ${key}.hic.fa.gz
        wait_on_openfile.sh ${key}.hic.fa.gz
        """
    }
}


/**
 * Assemble WGS reads
 */
(wgs_out, asm_in) = wgs_out.into(2)
asm_in = asm_in.map{[it[0], it[1], it[2], it[3]]}

// assemble each metagenome
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
 * Make Truth Tables
 */
(asm_out, truth_in) = asm_out.into(2)
truth_in = truth_in.map{ [it[0], it[1], it[4]] }

// generate ground-truth tables by aligning assembly contigs to community references
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
// join 3C reads and assembly results
hicmap_in = ms.joinChannels(hicmap_in, hic_out, 2)
        .map{ [it[0], it[1], it[-1]] }

// map 3C reads to assembly contigs
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
graph_in = graph_in.map{[it[0], it[1], it[2], it[3]]}

// contig graphs are generated from 3C mappings
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
wgsmap_in = wgsmap_in.map{ [it[0], it[1], it[2], it[3]]}

// map WGS reads to contigs
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
cov_in = cov_in.map{[it[0], it[1], it[2]]}

// depth inferred from WGS2CTG mapping
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
