include { INPUTS }      from '../modules/inputs.nf'
include { TARGETS_ZIP } from '../modules/targets_zip.nf'
include {STRELKA} from '../modules/strelka.nf'
include {NORMALIZE_VCF_PAIRS} from '../modules/normalize_vcf_pairs.nf'

workflow SOMATIC_PAIRS {

    // Step 1: parse and emit standardized input channels
    INPUTS()

    // Step 2: preprocess targets bed into bgzip/tabix pair
    TARGETS_ZIP(INPUTS.out.targets)

    ch_bams = INPUTS.out.samples

    ch_refs = INPUTS.out.refs
    ch_targets = TARGETS_ZIP.out.zipped_targets.first()

    STRELKA(ch_bams, ch_refs, ch_targets)

    raw_vcfs_pairs = STRELKA.out.strelka_snvs.mix(STRELKA.out.strelka_indels)

    NORMALIZE_VCF_PAIRS(raw_vcfs_pairs, ch_refs)
}
