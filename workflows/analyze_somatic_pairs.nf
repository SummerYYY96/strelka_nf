include { INPUTS }      from '../modules/inputs.nf'
include { TARGETS_ZIP } from '../modules/targets_zip.nf'
include {STRELKA} from '../modules/strelka.nf'
include {NORMALIZE_VCF_PAIRS} from '../modules/normalize_vcf_pairs.nf'
include {FILTER_VCF_PAIRS} from '../modules/filter_vcf_pairs.nf'
include {VCF_TO_TSV_PAIRS} from '../modules/vcf_to_tsv_pairs.nf'

workflow SOMATIC_PAIRS {

    // Step 1: parse and emit standardized input channels
    INPUTS()

    // Step 2: preprocess targets bed into bgzip/tabix pair
    TARGETS_ZIP(INPUTS.out.targets)

    ch_bams = INPUTS.out.samples

    ch_refs = INPUTS.out.refs
    ch_targets = TARGETS_ZIP.out.zipped_targets.first()

    // Step 3: call strelka reverse
    STRELKA(ch_bams, ch_refs, ch_targets)

    raw_vcfs_pairs = STRELKA.out.strelka_snvs.mix(STRELKA.out.strelka_indels)

    // Step 4: normalize VCF pairs 
    NORMALIZE_VCF_PAIRS(raw_vcfs_pairs, ch_refs)

    ch_norm_vcf_pairs = NORMALIZE_VCF_PAIRS.out.norm_vcfs_pairs

    // Step 5: filter vcf
    FILTER_VCF_PAIRS(ch_norm_vcf_pairs)

    ch_filter_vcf_pairs = FILTER_VCF_PAIRS.out.filtered_vcf_pairs
    
    // Step 6: transform vcf to tsv 
    VCF_TO_TSV_PAIRS(ch_filter_vcf_pairs, ch_refs)
}
