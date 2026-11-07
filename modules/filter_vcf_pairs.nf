process FILTER_VCF_PAIRS {
    // filter the .vcf for tumor-normal pairs
    // https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php
    tag "${caller}"
    publishDir "${params.outputDir}/variants/${caller}/filtered", mode: 'copy', pattern: "*${filtered_vcf}"

    input:
    tuple val(caller), val(callerType), val(sample), val(tumor), val(normal), path(norm_vcf)

    output:
    tuple val(caller), val(callerType), val(sample), val(tumor), val(normal), path(filtered_vcf), emit:filtered_vcf_pairs

    script:
    prefix = "${sample}.${caller}.${callerType}"
    filtered_vcf = "${prefix}.filtered.vcf"
    if( caller == 'Strelka' )
        """
        # Keep VCF header lines and records whose FILTER column is exactly PASS
        awk '/^#/ || \$7 == "PASS"' "${norm_vcf}" > "${filtered_vcf}"
        """
    else
        error "Invalid caller: ${caller}"
}