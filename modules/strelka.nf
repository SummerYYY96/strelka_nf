process STRELKA {
    
    publishDir "${params.outputDir}/variants/${caller}/raw", mode: 'copy'

    input:
    // A single channel emitting a tuple of three items from the samplesheet
    tuple val(sample), val(tumor), val(normal), path(tumorBam), path(tumorBai), path(normalBam), path(normalBai)
    tuple path(ref_fa), path(ref_fai), path(ref_dict)
    tuple path(targets_bgz), path(targets_tbi)

    output:
    tuple val("${caller}"), val("snvs"), val(sample), val(tumor), val(normal), path("${somatic_snvs}"), emit: strelka_snvs
    tuple val("${caller}"), val("indel"), val(sample), val(tumor), val(normal), path("${somatic_indels}"), emit: strelka_indels
    
    script:
    caller = "Strelka"
    prefix = "${sample}.${caller}"
    runDir = "${prefix}.Strelka"
    somatic_indels = "${prefix}.somatic.indels.vcf"
    somatic_indels_gz = "${prefix}.somatic.indels.vcf.gz"
    somatic_indels_tbi = "${prefix}.somatic.indels.vcf.gz.tbi"
    somatic_snvs = "${prefix}.somatic.snvs.vcf"
    somatic_snvs_gz = "${prefix}.somatic.snvs.vcf.gz"
    somatic_snvs_tbi = "${prefix}.somatic.snvs.vcf.gz.tbi"

    """
    # check points
    echo "Running Strelka process for comparisonID: ${sample}"
    echo "Tumor BAM: ${tumorBam}"
    echo "Reference Fasta: ${ref_fa}"

    configureStrelkaSomaticWorkflow.py \
    --normalBam "${normalBam}" \
    --tumorBam "${tumorBam}" \
    --referenceFasta "${ref_fa}" \
    --runDir ${runDir} \
    --callRegions "${targets_bgz}" \
    --exome

    python ${runDir}/runWorkflow.py \
    -m local \
    -j \${NSLOTS:-\${NTHREADS:-1}}

    mv ${runDir}/results/variants/somatic.indels.vcf.gz \
    "${somatic_indels_gz}"
    gunzip -c "${somatic_indels_gz}" > "${somatic_indels}"

    mv ${runDir}/results/variants/somatic.snvs.vcf.gz \
    "${somatic_snvs_gz}"
    gunzip -c "${somatic_snvs_gz}" > "${somatic_snvs}"
    """
}
