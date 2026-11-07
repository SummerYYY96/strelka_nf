process NORMALIZE_VCF_PAIRS {

    publishDir "${params.outputDir}/variants/${caller}/normalized", mode: 'copy'

    input:
    tuple val(caller), val(callerType), val(sample), val(tumor), val(normal), file(vcf)
    tuple path(ref_fa), path(ref_fai), path(ref_dict)

    output:
    tuple val(caller), val(callerType), val(sample), val(tumor), val(normal), file("${norm_vcf}"), emit: norm_vcfs_pairs

    script:
    prefix = "${sample}.${caller}.${callerType}"
    norm_vcf = "${prefix}.norm.vcf"

    """
    cat ${vcf} | \
    bcftools norm --multiallelics -both --output-type v - | \
    bcftools norm --fasta-ref "${ref_fa}" --output-type v - > \
    "${norm_vcf}"
    """
}