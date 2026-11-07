process VCF_TO_TSV_PAIRS {

    tag "${caller}"
    publishDir "${params.outputDir}/variants/${caller}/tsv", mode: 'copy', pattern: "*${tsv_file}"
    publishDir "${params.outputDir}/variants/${caller}/tsv", mode: 'copy', pattern: "*${reformat_tsv}"

    input:
    tuple val(caller), val(callerType), val(sample), val(tumor), val(normal), path(vcf) 
    tuple path(ref_fa), path(ref_fai), path(ref_dict)

    output:
    tuple val(caller), val(callerType), val(sample), val(tumor), val(normal), path(vcf), path(reformat_tsv), emit: vcf_tsv_pairs

    script:
    prefix = "${sample}.${caller}.${callerType}"
    tsv_file = "${prefix}.tsv"
    reformat_tsv = "${prefix}.reformat.tsv"
    if( caller == 'Strelka' )
        if ( callerType == "snvs" )
            """
            # convert to tsv format
            # NOTE: automatically filters for only PASS entries
            gatk.sh -T VariantsToTable \
            -R "${ref_fa}" \
            -V "${vcf}" \
            -F CHROM \
            -F POS \
            -F ID \
            -F REF \
            -F ALT \
            -F FILTER \
            -F DP \
            -F SOMATIC \
            -F QSS \
            -F MQ \
            -F SNVSB \
            -F SomaticEVS \
            -GF DP \
            -GF AU \
            -GF TU \
            -GF CU \
            -GF GU \
            -o "${tsv_file}"

            # vcf header example for Strelka SNVs
            ##fileformat=VCFv4.1
            ##FILTER=<ID=PASS,Description="All filters passed">
            ##source=strelka
            ##source_version=2.9.10
            ##content=strelka somatic snv calls
            ##priorSomaticSnvRate=0.0001
            ##INFO=<ID=QSS,Number=1,Type=Integer,Description="Quality score for any somatic snv, ie. for the ALT allele to be present at a significantly different frequency in the tumor and normal">
            ##INFO=<ID=TQSS,Number=1,Type=Integer,Description="Data tier used to compute QSS">
            ##INFO=<ID=NT,Number=1,Type=String,Description="Genotype of the normal in all data tiers, as used to classify somatic variants. One of {ref,het,hom,conflict}.">
            ##INFO=<ID=QSS_NT,Number=1,Type=Integer,Description="Quality score reflecting the joint probability of a somatic variant and NT">
            ##INFO=<ID=TQSS_NT,Number=1,Type=Integer,Description="Data tier used to compute QSS_NT">
            ##INFO=<ID=SGT,Number=1,Type=String,Description="Most likely somatic genotype excluding normal noise states">
            ##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation">
            ##INFO=<ID=DP,Number=1,Type=Integer,Description="Combined depth across samples">
            ##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
            ##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
            ##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref read-position in the tumor">
            ##INFO=<ID=SNVSB,Number=1,Type=Float,Description="Somatic SNV site strand bias">
            ##INFO=<ID=PNOISE,Number=1,Type=Float,Description="Fraction of panel containing non-reference noise at this site">
            ##INFO=<ID=PNOISE2,Number=1,Type=Float,Description="Fraction of panel containing more than one non-reference noise obs at this site">
            ##INFO=<ID=SomaticEVS,Number=1,Type=Float,Description="Somatic Empirical Variant Score (EVS) expressing the phred-scaled probability of the call being a false positive observation.">
            ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth for tier1 (used+filtered)">
            ##FORMAT=<ID=FDP,Number=1,Type=Integer,Description="Number of basecalls filtered from original read depth for tier1">
            ##FORMAT=<ID=SDP,Number=1,Type=Integer,Description="Number of reads with deletions spanning this site at tier1">
            ##FORMAT=<ID=SUBDP,Number=1,Type=Integer,Description="Number of reads below tier1 mapping quality threshold aligned across this site">
            ##FORMAT=<ID=AU,Number=2,Type=Integer,Description="Number of 'A' alleles used in tiers 1,2">
            ##FORMAT=<ID=CU,Number=2,Type=Integer,Description="Number of 'C' alleles used in tiers 1,2">
            ##FORMAT=<ID=GU,Number=2,Type=Integer,Description="Number of 'G' alleles used in tiers 1,2">
            ##FORMAT=<ID=TU,Number=2,Type=Integer,Description="Number of 'T' alleles used in tiers 1,2">
            ##FILTER=<ID=LowEVS,Description="Somatic Empirical Variant Score (SomaticEVS) is below threshold">
            ##FILTER=<ID=LowDepth,Description="Tumor or normal sample read depth at this locus is below 2">

            # reformat and adjust the TSV table for consistency downstream
            # add extra columns to the VCF TSV file for downstream
            reformat-vcf-table.py -c StrelkaSomaticSNV -s "${tumor}" -i "${tsv_file}" | \
            paste-col.py --header "Sample" -v "${tumor}" | \
            paste-col.py --header "Tumor" -v "${tumor}" | \
            paste-col.py --header "Normal" -v "${normal}" | \
            paste-col.py --header "VariantCallerType" -v "${callerType}"  | \
            paste-col.py --header "VariantCaller" -v "${caller}" > \
            "${reformat_tsv}"

            # make sure that the input and output tables have the same number of rows
            tsv_lines="\$(wc -l < "${tsv_file}" )"
            reformat_lines="\$( wc -l < "${reformat_tsv}" )"
            if [ "\$tsv_lines" -ne "\$reformat_lines" ]; then echo "ERROR: reformat table has different number of rows!"; exit 1; fi

            vcf_lines="\$(grep -v '^##' ${vcf} | wc -l)"
            if [ "\$vcf_lines" -ne "\$reformat_lines" ]; then echo "ERROR: reformat table has different number of entries than vcf file!"; exit 1; fi
            """
        else if( callerType == 'indel' )
            """
            # convert to tsv format
            # NOTE: automatically filters for only PASS entries
            gatk.sh -T VariantsToTable \
            -R "${ref_fa}" \
            -V "${vcf}" \
            -F CHROM \
            -F POS \
            -F ID \
            -F REF \
            -F ALT \
            -F FILTER \
            -F SOMATIC \
            -F MQ \
            -F SomaticEVS \
            -F QSI \
            -GF DP \
            -GF TAR \
            -GF TIR \
            -GF TOR \
            -o "${tsv_file}"

            # vcf header example for Strelka indels
            ##fileformat=VCFv4.1
            ##FILTER=<ID=PASS,Description="All filters passed">
            ##source=strelka
            ##source_version=2.9.10
            ##content=strelka somatic indel calls
            ##priorSomaticIndelRate=1e-06
            ##INFO=<ID=QSI,Number=1,Type=Integer,Description="Quality score for any somatic variant, ie. for the ALT haplotype to be present at a significantly different frequency in the tumor and normal">
            ##INFO=<ID=TQSI,Number=1,Type=Integer,Description="Data tier used to compute QSI">
            ##INFO=<ID=NT,Number=1,Type=String,Description="Genotype of the normal in all data tiers, as used to classify somatic variants. One of {ref,het,hom,conflict}.">
            ##INFO=<ID=QSI_NT,Number=1,Type=Integer,Description="Quality score reflecting the joint probability of a somatic variant and NT">
            ##INFO=<ID=TQSI_NT,Number=1,Type=Integer,Description="Data tier used to compute QSI_NT">
            ##INFO=<ID=SGT,Number=1,Type=String,Description="Most likely somatic genotype excluding normal noise states">
            ##INFO=<ID=RU,Number=1,Type=String,Description="Smallest repeating sequence unit in inserted or deleted sequence">
            ##INFO=<ID=RC,Number=1,Type=Integer,Description="Number of times RU repeats in the reference allele">
            ##INFO=<ID=IC,Number=1,Type=Integer,Description="Number of times RU repeats in the indel allele">
            ##INFO=<ID=IHP,Number=1,Type=Integer,Description="Largest reference interrupted homopolymer length intersecting with the indel">
            ##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
            ##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
            ##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation">
            ##INFO=<ID=OVERLAP,Number=0,Type=Flag,Description="Somatic indel possibly overlaps a second indel.">
            ##INFO=<ID=SomaticEVS,Number=1,Type=Float,Description="Somatic Empirical Variant Score (EVS) expressing the phred-scaled probability of the call being a false positive observation.">
            ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth for tier1">
            ##FORMAT=<ID=DP2,Number=1,Type=Integer,Description="Read depth for tier2">
            ##FORMAT=<ID=TAR,Number=2,Type=Integer,Description="Reads strongly supporting alternate allele for tiers 1,2">
            ##FORMAT=<ID=TIR,Number=2,Type=Integer,Description="Reads strongly supporting indel allele for tiers 1,2">
            ##FORMAT=<ID=TOR,Number=2,Type=Integer,Description="Other reads (weak support or insufficient indel breakpoint overlap) for tiers 1,2">
            ##FORMAT=<ID=DP50,Number=1,Type=Float,Description="Average tier1 read depth within 50 bases">
            ##FORMAT=<ID=FDP50,Number=1,Type=Float,Description="Average tier1 number of basecalls filtered from original read depth within 50 bases">
            ##FORMAT=<ID=SUBDP50,Number=1,Type=Float,Description="Average number of reads below tier1 mapping quality threshold aligned across sites within 50 bases">
            ##FORMAT=<ID=BCN50,Number=1,Type=Float,Description="Fraction of filtered reads within 50 bases of the indel.">
            ##FILTER=<ID=LowEVS,Description="Somatic Empirical Variant Score (SomaticEVS) is below threshold">
            ##FILTER=<ID=LowDepth,Description="Tumor or normal sample read depth at this locus is below 2">

            reformat-vcf-table.py -c StrelkaSomaticIndel -s "${tumor}" -i "${tsv_file}" | \
            paste-col.py --header "Sample" -v "${tumor}" | \
            paste-col.py --header "Tumor" -v "${tumor}" | \
            paste-col.py --header "Normal" -v "${normal}" | \
            paste-col.py --header "VariantCallerType" -v "${callerType}"  | \
            paste-col.py --header "VariantCaller" -v "${caller}" > \
            "${reformat_tsv}"

            # make sure that the input and output tables have the same number of rows
            if [ "\$(wc -l < "${tsv_file}" )" -ne "\$( wc -l < "${reformat_tsv}" )" ]; then echo "ERROR: reformat table has different number of rows!"; exit 1; fi
            """
        else
            error "Invalid Strelka callerType: ${callerType}"
    else
        error "Invalid caller: ${caller}"
}
