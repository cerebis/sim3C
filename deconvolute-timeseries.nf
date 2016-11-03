#!/usr/bin/env nextflow
/* deconvolution sweep.
** test metagenomic deconvolution over a range of simulated strain evolution parameters
**
** Usage: deconvolute-timeseries.nf -c sweep.config
*/
import static Helper.*

sweep = new Sweep()

sweep['ancestor'] = files(params.ancestor)
sweep['donor'] = files(params.donor)
sweep['alpha'] = stringToList(params.alpha)
sweep['tree'] = absPath(params.trees)
sweep['timepoints'] = stringToList(params.num_samples)
sweep['profile'] = absPath(params.profiles)
sweep['xfold'] = stringToList(params.xfold)

println sweep.description()

evo_in = sweep.permutedChannel('ancestor', 'donor', 'alpha', 'tree')

process Evolve {
    cache 'deep'

    input:
    set key, ancestor, donor, alpha, tree from evo_in

    output:
    set key, file("${key}.evo.fa") into evo_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.evo.fa
        """
    }
    else {
        """
        scale_tree.py -a $alpha $tree scaled_tree
        \$EXT_BIN/sgevolver/sgEvolver --indel-freq=${params.indel_freq} --small-ht-freq=${params.small_ht_freq} \
            --large-ht-freq=${params.large_ht_freq} --inversion-freq=${params.inversion_freq} \
            --random-seed=${params.seed} scaled_tree \
             $ancestor $donor ${key}.evo.aln ${key}.evo.fa
        strip_semis.sh ${key}.evo.fa
        """
    }
}

// add a name to new output
evo_out = evo_out.map { it.nameify(1, 'ref_seq') }

/**
 *  Make WGS reads
 */
(evo_out, wgs_in) = evo_out.into(2)
// add abundance profile, wgs coverage, and timepoints to initial sweep
wgs_in = sweep.extendChannel(wgs_in, 'profile', 'xfold', 'timepoints')

process WGS_Reads {
    cache 'deep'
    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    set key, ref_seq from wgs_in

    output:
    set key, file("${key}.wgs.*.r1.fq.gz"), file("${key}.wgs.*.r2.fq.gz"), file("${key}.coverage.tsv"), ref_seq into wgs_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.wgs.1.r1.fq.gz
        echo $key > ${key}.wgs.1.r2.fq.gz
        """
    }
    else {
        """
        export PATH=\$EXT_BIN/art:\$PATH
        metaART.py -C gzip -t ${key['profile']} -z ${key['timepoints']} -M ${key['xfold']} -S ${params.seed} -s ${params.wgs_ins_std} \
                -m ${params.wgs_ins_len} -l ${params.wgs_read_len} -n "${key}.wgs" $ref_seq .
	    wait_on_openfile.sh ${key}.wgs.r1.fq.gz
	    wait_on_openfile.sh ${key}.wgs.r2.fq.gz
        mv coverage.tsv ${key}.coverage.tsv
        """
    }
}

// add a name to new output
//wgs_out = wgs_out.map { it.nameify(1, 'wgs_reads') }


/**
 * Map WGS reads to contigs
 */

(wgs_out, map_in) = wgs_out.into(2)

process ReadMap {
    cache 'deep'
    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(reads1), file(reads2), coverage_file, ref_seq from map_in

    output:
    set key, file("*.bam"), reads1, reads2, ref_seq into map_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.bam
        """
    }
    else {
        fff = files(params.ancestor)[0]
        """
        export PATH=\$EXT_BIN/a5/bin:\$PATH

	    for rr in `ls *.r1.fq.gz`
        do
        	rbase=`basename \$rr .r1.fq.gz`
        	r2=\$rbase.r2.fq.gz
        	bwa mem -t 1 $fff \$rr \$r2 | samtools view -bS - | samtools sort -l 9 - \$rbase
        done
        """
    }
}


/**
 * Deconvolve the SNVs into strain genotypes
 */
 
 (map_out, deconv_in) = map_out.into(2)

process Deconvolve {
    cache 'deep'
    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(bams), reads1, reads2, ref_seq from deconv_in

    output:
    set key, file("${key}.decon.csv"), file("${key}.snv_file.data.R"), file("${key}.strains.tre"), bams, reads1, reads2, ref_seq into deconv_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.decon.csv
        echo $key > ${key}.snv_file.data.R
        echo $key > ${key}.strains.tre
        """
    }
    else {
        fff = files(params.ancestor)[0]
        """
        export PATH=\$EXT_BIN/lofreq_star:\$PATH
        snvbpnmft.py . $fff *.bam
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
(evo_out, truth_in) = evo_out.into(2)

process Truth {
    cache 'deep'
    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    set key, ref_seq from truth_in

    output:
    set key, file("${key}.truth.tsv") into truth_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.truth.tsv
        """
    }
    else {
        fff = files(params.ancestor)[0]
        """
        strain_truth.py --mauve-path=\$MAUVEPATH -o ${key}.truth.tsv ${ref_seq} ${fff}
        """
    }
}

(truth_out, accuracy_in) = truth_out.into(2)
(deconv_out, acc_deconv_in) = deconv_out.into(2)

accuracy_in = accuracy_in
        .cross(acc_deconv_in) { it -> ((Key)it[0]).id(4) }
        .map{ t -> [t[1][0], t[0][1]] + t[1][1..3] }

/**
 * Measure accuracy of strain genotypes
 */
process Accuracy {
    cache 'deep'
    publishDir params.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file(truthfile), file(snvbpnmf), file(snv_file), file(tree_file) from accuracy_in

    output:
    set file("${key}.truth.report.txt") into truth

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
