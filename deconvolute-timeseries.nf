#!/usr/bin/env nextflow
/* deconvolution sweep.
** test metagenomic deconvolution over a range of simulated strain evolution parameters
**
** Usage: deconvolute-timeseries.nf
*/
import MetaSweeper
// TODO this import and its dependent method below will eventually be moved to MetaSweeper
import groovy.util.GroovyCollections

MetaSweeper ms = MetaSweeper.fromFile(new File('timeseries.yaml'))

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
// convert the resulting map of sweep point results into table format
        .flatMap{ it.collect { k,v -> [k] +  v.flatten() } }
// include just the clade abundance profiles from each sweep point
        .map { it.pick([2, 4]) }

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
// convert the resulting map of sweep point results into table format
        .flatMap { it.collect { k, v -> [k] +  v.collect { cl -> [cl[1], cl[3]] }.flatten() } }
// include just the clade sequences from each sweep point
        .map { it.pick([1, 3]) }

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
    set key, file("${key}.wgs.*.r1.fq.gz"), file("${key}.wgs.*.r2.fq.gz"), file(comm_seq), file("${key}.cov") into wgs_out

    script:
    if (params.debug) {
        """
        for ((n=1; n<=$ms.options.num_samples; n++))
        do
            echo "metaART.py -C gzip --profile $comm_prof -z $ms.options.num_samples -M $xfold -S ${key['seed']} \
                    -s $ms.options.wgs.ins_std -m $ms.options.wgs.ins_len -l $ms.options.wgs.read_len
                    --coverage-out ${key}.cov -n ${key}.wgs $comm_seq ." > ${key}.wgs.\$n.r1.fq.gz

            echo "metaART.py -C gzip --profile $comm_prof -z $ms.options.num_samples -M $xfold -S ${key['seed']} \
                    -s $ms.options.wgs.ins_std -m $ms.options.wgs.ins_len -l $ms.options.wgs.read_len \
                    --coverage-out ${key}.cov -n ${key}.wgs $comm_seq ." > ${key}.wgs.\$n.r2.fq.gz
        done

        touch ${key}.cov
        """
    }
    else {
        """
        export PATH=\$EXT_BIN/art:\$PATH
        metaART.py -C gzip --profile $comm_prof -z $ms.options.num_samples -M $xfold -S ${key['seed']} \
                -s $ms.options.wgs.ins_std -m $ms.options.wgs.ins_len -l $ms.options.wgs.read_len \
                --coverage-out ${key}.cov -n ${key}.wgs $comm_seq .
        #wait_on_openfile.sh ${key}.wgs.*.r1.fq.gz
        #wait_on_openfile.sh ${key}.wgs.*.r2.fq.gz
        """
    }
}

/**
 * Map WGS reads to reference sequences
 */
// ancestral sequence for community
ancestor_in = Channel.value(file(ms.variables.community.clades[0].ancestor))

(wgs_out, map_in) = wgs_out.into(2)

        // ref, R1, R2
map_in = map_in.map{ it.pick(1, 2) }
        // pair the R1s and R2s -- might move this to a method in MetaSweeper
        // then flatten nested lists to become one row per read-pair
        .flatMap{ GroovyCollections.transpose(it[1].sort(), it[2].sort())
                    // extend the key to include sample number, extracting it from read-pair file names
                    .collect{ pair ->
                        def nsamp = pair.collect { ri -> (ri =~ /wgs\.(\d+)\.r[12].fq.*$/)[0][1] }
                        assert nsamp.size() == 2 && nsamp[0] == nsamp[1] : 'Error: read-pairs do not share the same sample index'
                        [ms.extendKey(it.getKey(), 'nsamp', nsamp[0]), *pair] } }

process WGSMap {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(reads1), file(reads2) from map_in
    file(anc) from ancestor_in

    output:
    set key, file("${key}.wgs2ref.bam"), file(reads1), file(reads2) into map_out

    script:
    if (params.debug) {
        """
        if [ ! -e $anc ]; then echo "no ancestral sequence found"; exit 1; fi
        echo $key > ${key}.wgs2ref.bam
        """
    }
    else {
        """
        export PATH=\$EXT_BIN/a5/bin:\$PATH
        if [ ! -e "${anc}.bwt" ]
        then
            bwa index $anc
        fi
        bwa mem -t 1 $anc $reads1 $reads2 | samtools view -bS - | samtools sort -l 9 - ${key}.wgs2ref
        """
    }
}

/**
 * Deconvolve the SNVs into strain genotypes
 */
// ancestral sequence for community
ancestor_in = Channel.value(file(ms.variables.community.clades[0].ancestor))

(map_out, deconv_in) = map_out.into(2)
            // just the bam files
deconv_in = deconv_in.map { it.pick(1) }
                // remove the sample number from key
                .map { [it.getKey().popLevels(1), *it.dropKey()] }
                // group each samples time-series bams on the new reduced key and sort by base filename
                .groupTuple(sort: {it.name})

process Deconvolve {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('tp*.bam') from deconv_in
    file(anc) from ancestor_in

    output:
    set key, file("${key}.decon.csv"), file("${key}.snv_file.data.R"), file("${key}.strains.tre") into deconv_out

    script:
    if (params.debug) {
        """
        if [ ! -e $anc]; then exit 1; fi
        echo $key > ${key}.decon.csv
        echo $key > ${key}.snv_file.data.R
        echo $key > ${key}.strains.tre
        """
    }
    else {
        """
        export PATH=\$EXT_BIN/lofreq_star:\$PATH
        snvbpnmft.py . $anc *.bam
        mv decon.csv ${key}.decon.csv
        mv elbos.csv ${key}.elbos.csv
        mv snv_file.data.R ${key}.snv_file.data.R
        java -Xmx1000m \$JNI_LIB_PATH -jar \$EXT_BIN/beast/beast.jar beast.xml
        java -jar \$EXT_BIN/treeanno/treeannotator.jar -burnin 1000 -heights mean aln.trees ${key}.strains.tre
        """
    }
}

/**
 * Record the true strain genotypes
 */
(seq_prof, truth_in) = seq_prof.into(2)

// ancestral sequence for community
ancestor_in = Channel.value(file(ms.variables.community.clades[0].ancestor))

// just the community reference sequence
truth_in = truth_in.map{ it.pick(1) }

process Truth {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(ref) from truth_in
    file(anc) from ancestor_in

    output:
    set key, file("${key}.truth.tsv") into truth_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.truth.tsv
        """
    }
    else {
        """
        strain_truth.py --mauve-path=\$MAUVEPATH -o ${key}.truth.tsv $ref $anc
        """
    }
}

/**
 * Measure accuracy of strain genotypes
 */

// join truth and deconv outputs at sweep depth of 2.
accuracy_in = ms.joinChannels(truth_out, deconv_out, 2)

process Accuracy {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(truthfile), file(snvbpnmf), file(snv_file), file(tree_file) from accuracy_in

    output:
    set file("${key}.truth.report.txt") into accuracy_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.truth.report.txt
        """
    }
    else {
        """
        measure_accuracy.py --bpnmf=${snvbpnmf} --truth=${truthfile} --sites=${snv_file} > ${key}.truth.report.txt
        """
    }
}
