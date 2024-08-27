params.ref_dir = "/gpfs/data/molecpathlab/ref"
// params.gatk_bundle_dir = "${params.ref_dir}/gatk-bundle"

params.targetbed = "/gpfs/data/molecpathlab/ref/targets/targets.629.bed"
params.ref_fa = "${params.ref_dir}/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
params.ref_fai = "${params.ref_dir}/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai"
params.ref_dict = "${params.ref_dir}/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.dict"

// params.germline_resource_gz = "${params.gatk_bundle_dir}/af-only-gnomad.raw.sites.hg19.vcf.gz"
// params.germline_resource_gz_tbi = "${params.gatk_bundle_dir}/af-only-gnomad.raw.sites.hg19.vcf.gz.tbi"

params.samplesheet = "sample.pairs.tsv"
params.outputDir = "output"

// Set up input data channels for targets, ref and germline vcf files
Channel.fromPath( file(params.targetbed) ).into { targets_bed; targets_bed2 }

Channel.fromPath( file(params.ref_fa) ).into { ref_fasta; ref_fasta2; ref_fasta3; ref_fasta4; ref_fasta5 }
Channel.fromPath( file(params.ref_fai) ).into { ref_fai; ref_fai2; ref_fai3; ref_fai4; ref_fai5 }
Channel.fromPath( file(params.ref_dict) ).into { ref_dict; ref_dict2; ref_dict3; ref_dict4; ref_dict5 }

// Channel.fromPath( file(params.germline_resource_gz) ).set { germline_resource_gz }
// Channel.fromPath( file(params.germline_resource_gz_tbi) ).set { germline_resource_gz_tbi }

// Read in sample paris file with sample id and bams
Channel.fromPath( file(params.samplesheet) )
       .splitCsv(header: true, sep: '\t')
       .map{row ->
         def comparisonID = row['Sample']
         def tumorID = row['Tumor']
         def normalID = row['Normal']
         def tumorBam = row['Tumor_Bam'].tokenize( ',' ).collect { file(it) }
         def tumorBai = row['Tumor_Bai'].tokenize( ',' ).collect { file(it) }
         def normalBam = row['Normal_Bam'].tokenize( ',' ).collect { file(it) }
         def normalBai = row['Normal_Bai'].tokenize( ',' ).collect { file(it) }
         return [ comparisonID, tumorID, normalID, tumorBam, tumorBai, normalBam, normalBai ]
       }
       .tap { samples_bam_bai;  samples_bam_bai2}

samples_bam_bai.combine(ref_fasta)
               .combine(ref_fai)
               .combine(ref_dict)
              //  .combine(targets_bed)
              //  .combine(germline_resource_gz)
              //  .combine(germline_resource_gz_tbi)
               .set { sample_pairs_ref }


process targets_zip {
    input:
    file(targets_bed) from targets_bed

    output:
    set file("${output_bgz}"), file("${output_index}") into targets_zipped, targets_zipped2

    script:
    output_bgz = "targets.bed.bgz"
    output_index = "targets.bed.bgz.tbi"
    """
    sort -V -k1,1 -k2,2 "${targets_bed}" > targets.sorted.bed
    bgzip -c targets.sorted.bed > "${output_bgz}"
    tabix -p bed "${output_bgz}"
    """
}

process manta {
    publishDir "${params.outputDir}/variants/${caller}/raw", mode: 'copy'

    input:
    set val(comparisonID), val(tumorID), val(normalID), file(tumorBam), file(tumorBai), file(normalBam), file(normalBai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bgz), file(targets_tbi) from sample_pairs_ref.combine(targets_zipped)

    output:
    set val("${caller}"), val("${callerType}"), val(comparisonID), val(tumorID), val(normalID), val("${chunkLabel}"), file("${candidateSmallIndels_gz}"), file("${candidateSmallIndels_tbi}") into mantaToStrelka
    file("${candidateSV}")
    file("${diploidSV}")
    file("${somaticSV}")
    file("${candidateSmallIndels}")
    // set val("${caller}"), val("${callerType}"), val(comparisonID), val(tumorID), val(normalID), val("${chunkLabel}"), file("${candidateSmallIndels}"), file("${candidateSmallIndels_tbi}"), file("${candidateSV}"), file("${candidateSV_tbi}"), file("${diploidSV}"), file("${diploidSV_tbi}"), file("${somaticSV}"), file("${somaticSV_tbi}") into mantaOutput

    script:
    caller = "Manta"
    chunkLabel = "NA"
    callerType = "NA"
    prefix = "${comparisonID}.${caller}.${callerType}.${chunkLabel}"
    runDir = "${prefix}.Manta"
    candidateSmallIndels = "${prefix}.candidateSmallIndels.vcf"
    candidateSmallIndels_gz = "${prefix}.candidateSmallIndels.vcf.gz"
    candidateSmallIndels_tbi = "${prefix}.candidateSmallIndels.vcf.gz.tbi"
    candidateSV = "${prefix}.candidateSV.vcf"
    candidateSV_gz = "${prefix}.candidateSV.vcf.gz"
    candidateSV_tbi = "${prefix}.candidateSV.vcf.gz.tbi"
    diploidSV = "${prefix}.diploidSV.vcf"
    diploidSV_gz = "${prefix}.diploidSV.vcf.gz"
    diploidSV_tbi = "${prefix}.diploidSV.vcf.gz.tbi"
    somaticSV = "${prefix}.somaticSV.vcf"
    somaticSV_gz = "${prefix}.somaticSV.vcf.gz"
    somaticSV_tbi = "${prefix}.somaticSV.vcf.gz.tbi"
    """
    # Debugging
    echo "Tumor BAM: ${tumorBam}"
    echo "Tumor BAI: ${tumorBai}"
    echo "Normal BAM: ${normalBam}"
    echo "Normal BAI: ${normalBai}"
    echo "Reference Fasta: ${ref_fasta}"
    echo "Reference Fai: ${ref_fai}"
    echo "Reference Dict: ${ref_dict}"
    echo "Targets BGZ: ${targets_bgz}"
    echo "Targets TBI: ${targets_tbi}"
    echo "Output directory is: ${params.outputDir}"

    configManta.py \
    --normalBam "${normalBam}" \
    --tumorBam "${tumorBam}" \
    --referenceFasta "${ref_fasta}" \
    --runDir "${runDir}" \
    --callRegions "${targets_bgz}" \
    --exome

    python ${runDir}/runWorkflow.py \
    -m local \
    -j \${NSLOTS:-\${NTHREADS:-1}}

    # needed for Strelka
    mv ${runDir}/results/variants/candidateSmallIndels.vcf.gz \
    "${candidateSmallIndels_gz}"
    gunzip -c "${candidateSmallIndels_gz}" > "${candidateSmallIndels}"
    mv ${runDir}/results/variants/candidateSmallIndels.vcf.gz.tbi \
    "${candidateSmallIndels_tbi}"

    mv ${runDir}/results/variants/candidateSV.vcf.gz \
    "${candidateSV_gz}"
    gunzip -c "${candidateSV_gz}" > "${candidateSV}"

    mv ${runDir}/results/variants/diploidSV.vcf.gz \
    "${diploidSV_gz}"
    gunzip -c "${diploidSV_gz}" > "${diploidSV}"

    mv ${runDir}/results/variants/somaticSV.vcf.gz \
    "${somaticSV_gz}"
    gunzip -c "${somaticSV_gz}" > "${somaticSV}"
    """
}

samples_bam_bai2.combine(mantaToStrelka)
.filter { items ->
    def comparisonID = items[0]
    def tumorID = items[1]
    def tumorBam = items[2]
    def tumorBai = items[3]
    def normalID = items[4]
    def normalBam = items[5]
    def normalBai = items[6]
    def MantaCaller = items[7]
    def MantaCallerType = items[8]
    def MantaComparisonID = items[9]
    def MantaTumorID = items[10]
    def MantaNormalID = items[11]
    def MantaChunkLabel = items[12]
    def MantaCandidateSmallIndels = items[13]
    def MantaCandidateSmallIndels_tbi = items[14]

    def comparisonID_match = comparisonID == MantaComparisonID
    def tumorID_match = tumorID == MantaTumorID
    def normalID_match = normalID == MantaNormalID
    def all_matches = [ comparisonID_match, tumorID_match, normalID_match ]
    for ( match in all_matches ){
        if ( match == false ){
            return false
        }
    }
    return true
}
.map { items ->
    def comparisonID = items[0]
    def tumorID = items[1]
    def tumorBam = items[2]
    def tumorBai = items[3]
    def normalID = items[4]
    def normalBam = items[5]
    def normalBai = items[6]
    def MantaCaller = items[7]
    def MantaCallerType = items[8]
    def MantaComparisonID = items[9]
    def MantaTumorID = items[10]
    def MantaNormalID = items[11]
    def MantaChunkLabel = items[12]
    def MantaCandidateSmallIndels = items[13]
    def MantaCandidateSmallIndels_tbi = items[14]

    return([comparisonID, tumorID, tumorBam, tumorBai, normalID, normalBam, normalBai, MantaCandidateSmallIndels, MantaCandidateSmallIndels_tbi])
}
.combine(ref_fasta2) // add reference genome and targets
.combine(ref_fai2)
.combine(ref_dict2)
.tap { samples_bam_bai_manta }

process strelka {
    publishDir "${params.outputDir}/variants/${caller}/raw", mode: 'copy'

    input:
    set val(comparisonID), val(tumorID), file(tumorBam), file(tumorBai), val(normalID), file(normalBam), file(normalBai), file(small_indels), file(small_indels_tbi), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bgz), file(targets_tbi) from samples_bam_bai_manta.combine(targets_zipped2)

    output:
    set val("${caller}"), val("snvs"), val(comparisonID), val(tumorID), val(normalID), val("${chunkLabel}"), file("${somatic_snvs}") into strelka_snvs
    set val("${caller}"), val("indel"), val(comparisonID), val(tumorID), val(normalID), val("${chunkLabel}"), file("${somatic_indels}") into strelka_indels

    script:
    caller = "Strelka"
    chunkLabel = "NA"
    callerType = "NA"
    prefix = "${comparisonID}.${caller}.${callerType}.${chunkLabel}"
    runDir = "${prefix}.Strelka"
    somatic_indels = "${prefix}.somatic.indels.vcf"
    somatic_indels_gz = "${prefix}.somatic.indels.vcf.gz"
    somatic_indels_tbi = "${prefix}.somatic.indels.vcf.gz.tbi"
    somatic_snvs = "${prefix}.somatic.snvs.vcf"
    somatic_snvs_gz = "${prefix}.somatic.snvs.vcf.gz"
    somatic_snvs_tbi = "${prefix}.somatic.snvs.vcf.gz.tbi"
    """
    configureStrelkaSomaticWorkflow.py \
    --normalBam "${normalBam}" \
    --tumorBam "${tumorBam}" \
    --referenceFasta "${ref_fasta}" \
    --indelCandidates "${small_indels}" \
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

strelka_snvs.mix(strelka_indels).set{ raw_vcfs_pairs }
process normalize_vcfs_pairs {
    publishDir "${params.outputDir}/variants/${caller}/normalized", mode: 'copy'

    input:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file(vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from raw_vcfs_pairs.combine(ref_fasta3).combine(ref_fai3).combine(ref_dict3)

    output:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file("${norm_vcf}") into norm_vcfs_pairs

    script:
    prefix = "${comparisonID}.${caller}.${callerType}.${chunkLabel}"
    norm_vcf = "${prefix}.norm.vcf"
        """
        cat ${vcf} | \
        bcftools norm --multiallelics -both --output-type v - | \
        bcftools norm --fasta-ref "${ref_fasta}" --output-type v - > \
        "${norm_vcf}"
        """
}

// get all the paired sample vcfs for downstream processing
// vcfs_mutect2.mix(
//     vcfs_lofreq_somatic_snvs_vcf_norm,
//     vcfs_lofreq_somatic_indels_vcf_norm,
//     // vcfs_lofreq_somatic_snvs_minus_dbsnp_vcf_norm,
//     // vcfs_lofreq_somatic_indels_minus_dbsnp_vcf_norm,
//     norm_vcfs_pairs).set { vcfs_pairs }

process filter_vcf_pairs {
    // filter the .vcf for tumor-normal pairs
    // https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php
    tag "${caller}.${chunkLabel}"
    publishDir "${params.outputDir}/variants/${caller}/filtered", mode: 'copy', pattern: "*${filtered_vcf}"

    input:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file(vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from norm_vcfs_pairs.combine(ref_fasta4).combine(ref_fai4).combine(ref_dict4)

    output:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file("${filtered_vcf}") into filtered_vcf_pairs // to tsv

    script:
    prefix = "${comparisonID}.${caller}.${callerType}.${chunkLabel}"
    filtered_vcf = "${prefix}.filtered.vcf"
    if( caller == 'MuTect2' )
        """
        # filter VCF
        # report if:
        # only keep 'PASS' entries

        # get the header
        # grep '^#' "${vcf}" > "${filtered_vcf}"
        # get the 'PASS' entries
        # grep -v '^#' "${vcf}" | grep 'PASS' >> "${filtered_vcf}" || :

        # old method
        # gatk.sh -T SelectVariants \
        # -R "${ref_fasta}" \
        # -V "${vcf}" \
        # -select 'vc.isNotFiltered()' \
        # > "${filtered_vcf}"

        # other criteria:

        # T frequency is more than 3%
        # ( Tumor Allelic depth alt / (Tumor Allelic depth ref + Tumor Allelic depth alt ) )  > 0.03
        # -select "(vc.getGenotype('TUMOR').getAD().1 / (vc.getGenotype('TUMOR').getAD().0 + vc.getGenotype('TUMOR').getAD().1) )  > 0.03" \

        # N frequency is less than 5%
        # ( Normal Allelic depth alt / ( Normal Allelic depth ref + Normal Allelic depth alt ) )  < 0.05
        # -select "(vc.getGenotype('NORMAL').getAD().1 / (vc.getGenotype('NORMAL').getAD().0 + vc.getGenotype('NORMAL').getAD().1) )  < 0.05" \

        # at least 5 variant call supporting reads
        # Tumor Allelic depth alt > 5
        # -select "vc.getGenotype('TUMOR').getAD().1 > 5" \

        # T frequency is sufficiently higher (5x) than N frequency
        # "we recommend applying post-processing filters, e.g. by hard-filtering calls with low minor allele frequencies"
        # ( Tumor Allelic depth alt / ( Tumor Allelic depth ref + Tumor Allelic depth alt ) ) > ( Normal Allelic depth alt / ( Normal Allelic depth ref + Normal Allelic depth alt ) ) * 5
        # -select "(vc.getGenotype('TUMOR').getAD().1 / (vc.getGenotype('TUMOR').getAD().0 + vc.getGenotype('TUMOR').getAD().1) ) > (vc.getGenotype('NORMAL').getAD().1 / (vc.getGenotype('NORMAL').getAD().0 + vc.getGenotype('NORMAL').getAD().1) ) * 5" \

        # vc.getGenotype('TUMOR').getAD().0 ; Tumor Allelic depth ref
        # vc.getGenotype('TUMOR').getAD().1 ; Tumor Allelic depth alt
        # vc.getGenotype('NORMAL').getAD().0 ; Normal Allelic depth ref
        # vc.getGenotype('NORMAL').getAD().1 ; Normal Allelic depth alt

        # example variant format:
        ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
        # CHROM	chr7
        # POS	2946342
        # ID	.
        # REF	A
        # ALT	G
        # QUAL	.
        # FILTER	clustered_events;homologous_mapping_event
        # INFO	ECNT=16;HCNT=3;MAX_ED=65;MIN_ED=5;NLOD=33.41;TLOD=23.04
        # FORMAT	GT:AD:AF:ALT_F1R2:ALT_F2R1:FOXOG:PGT:PID:QSS:REF_F1R2:REF_F2R1
        # TUMOR	0/1:1333,17:0.013:7:10:0.588:0|1:2946342_A_G:40125,535:641:689
        # NORMAL	0/0:137,0:0.00:0:0:.:0|1:2946342_A_G:3959,0:53:80

        # get the header
        grep '^#' "${vcf}" > "${filtered_vcf}"
        # get the 'PASS' entries
        grep -v '^#' "${vcf}" | grep 'PASS' >> "${filtered_vcf}" || :
        """
    else if( caller == 'LoFreqSomatic' )
        """
        # do not report if:
        # - frequency is less than 1%, greater than 99%
        # - depth less than 200
        # gatk.sh -T SelectVariants \
        # -R "${ref_fasta}" \
        # -V "${vcf}" \
        # -select "AF > 0.01"  \
        # -select "AF < 0.99"  \
        # -select "DP > 100"  \
        # > "${filtered_vcf}"

        # get the header
        grep '^#' "${vcf}" > "${filtered_vcf}"
        # get the 'PASS' entries
        grep -v '^#' "${vcf}" | grep 'PASS' >> "${filtered_vcf}" || :
        """
    else if( caller == 'Strelka' )
        """
        # only keep 'PASS' entries
        # filter out TQSS_NT=2 https://github.com/Illumina/strelka/issues/65
        # gatk.sh -T SelectVariants \
        # -R "${ref_fasta}" \
        # -V "${vcf}" \
        # -select "TQSS_NT != 2"  \
        # --excludeFiltered \
        # > "${filtered_vcf}"

        # get the header
        grep '^#' "${vcf}" > "${filtered_vcf}"
        # get the 'PASS' entries
        grep -v '^#' "${vcf}" | grep 'PASS' >> "${filtered_vcf}" || :
        """
        // ##FILTER=<ID=PASS,Description="All filters passed">
        // ##INFO=<ID=QSS,Number=1,Type=Integer,Description="Quality score for any somatic snv, ie. for the ALT allele to be present at a significantly different frequency in the tumor and normal">
        // ##INFO=<ID=TQSS,Number=1,Type=Integer,Description="Data tier used to compute QSS">
        // ##INFO=<ID=NT,Number=1,Type=String,Description="Genotype of the normal in all data tiers, as used to classify somatic variants. One of {ref,het,hom,conflict}.">
        // ##INFO=<ID=QSS_NT,Number=1,Type=Integer,Description="Quality score reflecting the joint probability of a somatic variant and NT">
        // ##INFO=<ID=TQSS_NT,Number=1,Type=Integer,Description="Data tier used to compute QSS_NT">
        // ##INFO=<ID=SGT,Number=1,Type=String,Description="Most likely somatic genotype excluding normal noise states">
        // ##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation">
        // ##INFO=<ID=DP,Number=1,Type=Integer,Description="Combined depth across samples">
        // ##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
        // ##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
        // ##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref read-position in the tumor">
        // ##INFO=<ID=SNVSB,Number=1,Type=Float,Description="Somatic SNV site strand bias">
        // ##INFO=<ID=PNOISE,Number=1,Type=Float,Description="Fraction of panel containing non-reference noise at this site">
        // ##INFO=<ID=PNOISE2,Number=1,Type=Float,Description="Fraction of panel containing more than one non-reference noise obs at this site">
        // ##INFO=<ID=SomaticEVS,Number=1,Type=Float,Description="Somatic Empirical Variant Score (EVS) expressing the phred-scaled probability of the call being a false positive observation.">
        // ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth for tier1 (used+filtered)">
        // ##FORMAT=<ID=FDP,Number=1,Type=Integer,Description="Number of basecalls filtered from original read depth for tier1">
        // ##FORMAT=<ID=SDP,Number=1,Type=Integer,Description="Number of reads with deletions spanning this site at tier1">
        // ##FORMAT=<ID=SUBDP,Number=1,Type=Integer,Description="Number of reads below tier1 mapping quality threshold aligned across this site">
        // ##FORMAT=<ID=AU,Number=2,Type=Integer,Description="Number of 'A' alleles used in tiers 1,2">
        // ##FORMAT=<ID=CU,Number=2,Type=Integer,Description="Number of 'C' alleles used in tiers 1,2">
        // ##FORMAT=<ID=GU,Number=2,Type=Integer,Description="Number of 'G' alleles used in tiers 1,2">
        // ##FORMAT=<ID=TU,Number=2,Type=Integer,Description="Number of 'T' alleles used in tiers 1,2">
        // ##FILTER=<ID=LowEVS,Description="Somatic Empirical Variant Score (SomaticEVS) is below threshold">
        // ##FILTER=<ID=LowDepth,Description="Tumor or normal sample read depth at this locus is below 2">

        // ##FILTER=<ID=PASS,Description="All filters passed">
        // ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
        // ##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
        // ##INFO=<ID=PF,Number=1,Type=Integer,Description="The number of samples carry the variant">
        // ##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
        // ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
        // ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
        // ##INFO=<ID=NTLEN,Number=.,Type=Integer,Description="Number of bases inserted in place of deleted code">
        // ##FORMAT=<ID=PL,Number=3,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
        // ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        // ##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Reference depth, how many reads support the reference">
        // ##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allele depth, how many reads support this allele">
    else
        error "Invalid caller: ${caller}"
}

process vcf_to_tsv_pairs {
    tag "${caller}.${chunkLabel}"
    publishDir "${params.outputDir}/variants/${caller}/tsv", mode: 'copy', pattern: "*${tsv_file}"
    publishDir "${params.outputDir}/variants/${caller}/tsv", mode: 'copy', pattern: "*${reformat_tsv}"

    input:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file(vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from filtered_vcf_pairs.combine(ref_fasta5).combine(ref_fai5).combine(ref_dict5)

    output:
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file(vcf), file("${reformat_tsv}") into vcf_tsv_pairs // to annotation
    set val(caller), val(callerType), val(comparisonID), val(tumorID), val(normalID), val(chunkLabel), file("${reformat_tsv}") into vcf_tsv_pairs2

    script:
    prefix = "${comparisonID}.${caller}.${callerType}.${chunkLabel}"
    tsv_file = "${prefix}.tsv"
    reformat_tsv = "${prefix}.reformat.tsv"
    if( caller == 'MuTect2' )
        """
        # convert VCF to TSV
        # NOTE: automatically filters for only PASS entries
        gatk.sh -T VariantsToTable \
        -R "${ref_fasta}" \
        -V "${vcf}" \
        -F CHROM -F POS -F ID -F REF -F ALT -F FILTER -F QUAL -F AC -F AN -F NLOD -F TLOD \
        -GF AD -GF DP -GF AF \
        -o "${tsv_file}"

        # .vcf field descriptions:
        ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
        ##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele fraction of the event in the tumor">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
        ##INFO=<ID=NLOD,Number=1,Type=String,Description="Normal LOD score">
        ##INFO=<ID=TLOD,Number=1,Type=String,Description="Tumor LOD score">

        # reformat and adjust the TSV table for consistency downstream
        # add extra columns to the VCF TSV file for downstream
        reformat-vcf-table.py -c MuTect2 -s "${tumorID}" -i "${tsv_file}" | \
        paste-col.py --header "Sample" -v "${tumorID}"  | \
        paste-col.py --header "Tumor" -v "${tumorID}"  | \
        paste-col.py --header "Normal" -v "${normalID}"  | \
        paste-col.py --header "VariantCallerType" -v "${callerType}"  | \
        paste-col.py --header "VariantCaller" -v "${caller}" > \
        "${reformat_tsv}"

        # make sure that the input and output tables have the same number of rows
        if [ "\$(wc -l < "${tsv_file}" )" -ne "\$( wc -l < "${reformat_tsv}" )" ]; then echo "ERROR: reformat table has different number of rows!"; exit 1; fi
        """
    else if( caller == 'LoFreqSomatic' )
        """
        # convert to tsv format
        # NOTE: automatically filters for only PASS entries
        gatk.sh -T VariantsToTable \
        -R "${ref_fasta}" \
        -V "${vcf}" \
        -F CHROM \
        -F POS \
        -F ID \
        -F REF \
        -F ALT \
        -F QUAL \
        -F FILTER \
        -F DP \
        -F DP4 \
        -F AF \
        -F SB \
        -F UQ \
        -F CONSVAR \
        -F AN \
        -F AC \
        -F HRUN \
        -F INDEL \
        -F UNIQ \
        -F SOMATIC \
        -o "${tsv_file}"

        # .vcf field descriptions:
        ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
        ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
        ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
        ##INFO=<ID=CONSVAR,Number=0,Type=Flag,Description="Indicates that the variant is a consensus variant (as opposed to a low frequency variant).">
        ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
        ##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">
        ##INFO=<ID=DP4,Number=4,Type=Integer,Description="Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
        ##INFO=<ID=HRUN,Number=1,Type=Integer,Description="Homopolymer length to the right of report indel position">
        ##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
        ##INFO=<ID=SB,Number=1,Type=Float,Description="Strand Bias">
        ##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic event">
        ##INFO=<ID=UNIQ,Number=0,Type=Flag,Description="Unique, i.e. not detectable in paired sample">
        ##INFO=<ID=UQ,Number=1,Type=Integer,Description="Phred-scaled uniq score at this position">

        # reformat and adjust the TSV table for consistency downstream
        # add extra columns to the VCF TSV file for downstream
        reformat-vcf-table.py -c LoFreq -s "${tumorID}" -i "${tsv_file}" | \
        paste-col.py --header "Sample" -v "${tumorID}" | \
        paste-col.py --header "Tumor" -v "${tumorID}" | \
        paste-col.py --header "Normal" -v "${normalID}" | \
        paste-col.py --header "VariantCallerType" -v "${callerType}"  | \
        paste-col.py --header "VariantCaller" -v "${caller}" > \
        "${reformat_tsv}"

        # make sure that the input and output tables have the same number of rows
        if [ "\$(wc -l < "${tsv_file}" )" -ne "\$( wc -l < "${reformat_tsv}" )" ]; then echo "ERROR: reformat table has different number of rows!"; exit 1; fi
        """
    else if( caller == 'Strelka' )
        if ( callerType == "snvs" )
            """
            # convert to tsv format
            # NOTE: automatically filters for only PASS entries
            gatk.sh -T VariantsToTable \
            -R "${ref_fasta}" \
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
            reformat-vcf-table.py -c StrelkaSomaticSNV -s "${tumorID}" -i "${tsv_file}" | \
            paste-col.py --header "Sample" -v "${tumorID}" | \
            paste-col.py --header "Tumor" -v "${tumorID}" | \
            paste-col.py --header "Normal" -v "${normalID}" | \
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
            -R "${ref_fasta}" \
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

            reformat-vcf-table.py -c StrelkaSomaticIndel -s "${tumorID}" -i "${tsv_file}" | \
            paste-col.py --header "Sample" -v "${tumorID}" | \
            paste-col.py --header "Tumor" -v "${tumorID}" | \
            paste-col.py --header "Normal" -v "${normalID}" | \
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
// //sample_bam_pairs_ref_germlinevcf.subscribe { println "value $it" }
// process mutect2 {

//     publishDir "${params.outputdir}/mutect2", mode: 'copy', overwrite: true
//     echo true

//     input:
//     set val(sampleID), val(tumorID), val(normalID), file(tumorBam), file(tumorBai), file(normalBam), file(normalBai), file(ref_fasta), file(ref_fai), file(ref_dict), file(targets_bed), file(germline_resource_gz), file(germline_resource_gz_tbi) from sample_bam_pairs_ref_germlinevcf

//     output:
//     file("${vcf_file}")

//     script:
//     prefix = "${tumorID}.${normalID}"
//     vcf_file = "${prefix}.vcf"

//     """
//     gatk --java-options \"-Xms8G -Xmx8G\" Mutect2 \
//     --seconds-between-progress-updates 600 \
//     --native-pair-hmm-threads 4 \
//     --reference "${ref_fasta}" \
//     --germline-resource "${germline_resource_gz}" \
//     --dont-use-soft-clipped-bases \
//     --max-reads-per-alignment-start 100 \
//     --intervals "${targets_bed}" \
//     --interval-padding 10 \
//     --input "${tumorBam}" \
//     --input "${normalBam}" \
//     --tumor-sample "${tumorID}" \
//     --normal-sample "${normalID}" \
//     --output "${vcf_file}"
//     """
// }
