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

gen_in = ms.createSweep()
        .withVariable('seed')
        .withVariable('community',true)
        .permute()

/**
 * Generate reference genomes from ancestor
 */

process TreeGen {
    publishDir params.output, mode: 'copy', overwrite: 'true'

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
            --max-height 0.1 --birth-rate ${clade.value.tree['birth']} --death-rate ${clade.value.tree['death']} \
            --format newick --num-taxa ${clade.value.ntaxa} ${key}.nwk
        """
    }
}
tree_out = tree_out.map{ it.nameify(1, 'tree_file') }

(tree_out, evo_in) = tree_out.into(2)
// add abundance profile and wgs coverage to initial sweep
evo_in = ms.withVariable('alpha')
            .extend(evo_in, 'alpha')

process Evolve {
    publishDir params.output, mode: 'copy', overwrite: 'true'

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
        \$EXT_BIN/sgevolver/sgEvolver --indel-freq=${ms.options['evo']['indel_freq']} \
            --small-ht-freq=${ms.options['evo']['small_ht_freq']} \
            --large-ht-freq=${ms.options['evo']['large_ht_freq']} \
            --inversion-freq=${ms.options['evo']['inversion_freq']} \
            --random-seed=$seed scaled_tree \
             $clade.value.ancestor $clade.value.donor "${key}.evo.aln" "${key}.evo.fa"
        strip_semis.sh "${key}.evo.fa"
        """
    }

}

// add a name to new output
evo_out = evo_out.map { it.nameify(1, 'clade_seq') }

(evo_out, prof_in) = evo_out.into(2)
prof_in = prof_in.map{[it[0], it[1].value]}

process ProfileGen {
    publishDir params.output, mode: 'copy', overwrite: 'true'

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

prof_out = prof_out.map { it.nameify(1, 'clade_profile') }

(prof_out, merge_prof_in) = prof_out.into(2)
merge_prof_in = merge_prof_in.groupBy{ it[0].selectedKey('seed','alpha') }
        .flatMap{ it.collect { k,v -> [k] +  v.flatten() } }
        .map{ [it[0], [it[2],it[4]]*.value] }

process ProfileMerge {
    publishDir params.output, mode: 'copy', overwrite: 'true'

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

(evo_out, merge_seq_in) = evo_out.into(2)

// TODO
// - this needs to be cleaned up. We have to expose the underlying Path
// type to the Processor. Currently NamedValue gets in the way of this
// happening and we're forced to unwrap the objects.
// Either figure work-around or we should consider removing NamedValue

merge_seq_in = merge_seq_in.groupBy{ it[0].selectedKey('seed','alpha') }
        .flatMap{it.collect { k,v -> [k] +  v.collect{cl -> [cl[1],cl[3]]}.flatten() } }
        .map{ [it[0], [it[1],it[3]]*.value] }

process MergeClades {
    publishDir params.output, mode: 'copy', overwrite: 'true'

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

merge_seq_out = merge_seq_out.map { it.nameify(1, 'comm_seq') }

//
//Make WGS reads
//

(merge_seq_out, wgs_in) = merge_seq_out.into(2)
(merge_prof_out, tmp) = merge_prof_out.into(2)

wgs_in = wgs_in.map{ [it[0], it[1].value] }.phase(tmp).map{ [*it[0], it[1][1]]}
wgs_in = ms.withVariable('xfold').extend(wgs_in, 'xfold')

process WGS_Reads {
    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(comm_seq), file(comm_prof), xfold from wgs_in

    output:
    set key, file("${key}.wgs.r*.fq.gz"), file(comm_seq) into wgs_out

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


// add a name to new output
wgs_out = wgs_out.map { it.nameify(1, 'wgs_reads') }



// Make HiC reads
//

// add abundance profile and 3c depth to initial sweep
(merge_seq_out, hic_in) = merge_seq_out.into(2)
(merge_prof_out, tmp) = merge_prof_out.into(2)

hic_in = hic_in.map{ [it[0], it[1].value]}.phase(tmp).map{ [*it[0], it[1][1]]}
hic_in = ms.withVariable('n3c').extend(hic_in, 'n3c')

process HIC_Reads {
    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(comm_seq), file(comm_prof), n3c from hic_in

    output:
    set key, file("${key}.hic.fa.gz") into hic_out

    script:
    if (params.debug) {
        """
        echo "simForward.py -C gzip -r ${key['seed']} -n $n3c -l ${ms.options['n3c']['read_len']} \
            -p ${ms.options['n3c']['inter_prob']} --profile $comm_prof $comm_seq ${key}.hic.fa.gz" > ${key}.hic.fa.gz
        """
    }
    else {
        """
        simForward.py -C gzip -r ${key['seed']} -n $n3c -l ${ms.options['n3c']['read_len']} \
            -p ${ms.options['n3c']['inter_prob']} --profile $comm_prof $comm_seq "${key}.hic.fa.gz"
        wait_on_openfile.sh ${key}.hic.fa.gz
        """
    }
}


// add a name to new output
hic_out = hic_out.map { it.nameify(1, 'hic_reads') }


//
// Assemble WGS reads
//
(wgs_out, asm_in) = wgs_out.into(2)
asm_in = asm_in.map{[it[0], it[1].value.sort(), it[2]]}

process Assemble {
    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(reads), file(comm_seq) from asm_in

    output:
    set key, file("${key}.contigs.fasta"), reads, file(comm_seq) into asm_out

    script:
    if (params.debug) {
        """
        echo "\$EXT_BIN/a5/bin/a5_pipeline.pl --threads=1 --metagenome ${reads[0]} ${reads[1]} $key" > ${key}.contigs.fasta
        """
    }
    else {
        """
        export PATH=\$EXT_BIN/a5/bin:\$PATH
        a5_pipeline.pl --threads=1 --metagenome ${reads[0]} ${reads[1]} $key
        bwa index ${key}.contigs.fasta
        """
    }
}

// add a name to new output
asm_out = asm_out.map { it.nameify(1, 'contigs') }


//
// Make Truth Tables
//
(asm_out, truth_in) = asm_out.into(2)

truth_in = truth_in.map{ [it[0], it[1].value, it[3]] }

process Truth {
    publishDir params.output, mode: 'copy', overwrite: 'true'

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


// add a name to new output
truth_out = truth_out.map { it.nameify(1, 'truth') }

//
// Map HiC reads to assembled contigs
//
(asm_out, hicmap_in) = asm_out.into(2)
// combine results of hic and assembly processes, reduce to unique columns and select those relevant
hicmap_in = ms.sweep.joinChannels(hic_out, hicmap_in, 2)
        .map{ [it[0], it[1].value, it[2].value] }

process HiCMap {
    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(hic_reads), file(contigs) from hicmap_in

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
        bwa mem -t 1 $contigs $hic_reads | samtools view -bS - | samtools sort -l 9 - "${key}.hic2ctg"
        samtools index "${key}.hic2ctg.bam"
        samtools idxstats "${key}.hic2ctg.bam" > "${key}.hic2ctg.idxstats"
        samtools flagstat "${key}.hic2ctg.bam" > "${key}.hic2ctg.flagstat"
        """
    }
}

// add a name to new output
hicmap_out = hicmap_out.map { it.nameify(1, 'hic2ctg') }

//
// Generate contig graphs
//
(hicmap_out, graph_in) = hicmap_out.into(2)
graph_in = graph_in.map{[it[0], it[1].value, it[2], it[3]]}

process Graph {

    publishDir params.output, mode: 'copy', overwrite: 'true'

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
        bamToEdges_mod2.py --sim --afmt bam --strong 150 --graphml "${key}.graphml" --merged $hic2ctg hic2ctg.e hic2ctg.n
        """
    }
}

// add a name to new output
graph_out = graph_out.map { it.nameify(1, 'graph') }


//
// Map WGS reads to contigs
//
(asm_out, wgsmap_in) = asm_out.into(2)
wgsmap_in = wgsmap_in.map{ [it[0], it[1].value, it[2]]}

process WGSMap {
    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(contigs), file(reads) from wgsmap_in

    output:
    set key, file("${key}.wgs2ctg.bam"), file(contigs), reads into wgsmap_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.wgs2ctg.bam
        """
    }
    else {
        """
        export PATH=\$EXT_BIN/a5/bin:\$PATH
        bwa mem -t 1 $contigs ${reads[0]} ${reads[1]} | samtools view -bS - | samtools sort -l 9 - "${key}.wgs2ctg"
        """
    }
}

// add a name to new output
wgsmap_out = wgsmap_out.map { it.nameify(1, 'wgs2ctg') }


//
// Calculate assembly contig coverage
//
(wgsmap_out, cov_in) = wgsmap_out.into(2)


process InferReadDepth {
    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(wgs2ctg), file(contigs), file(reads) from cov_in

    output:
    set key, file("${key}.wgs2ctg.cov"), file(wgs2ctg), file(contigs), reads into cov_out

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

// add a name to new output
cov_out = cov_out.map { it.nameify(1, 'coverage') }
