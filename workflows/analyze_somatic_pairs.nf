include { INPUTS }      from '../modules/inputs.nf'
include { TARGETS_ZIP } from '../modules/targets_zip.nf'

workflow SOMATIC_PAIRS {

    // Step 1: parse and emit standardized input channels
    INPUTS()

    // Step 2: preprocess targets bed into bgzip/tabix pair
    zipped_targets = TARGETS_ZIP(INPUTS.out.targets)

    emit:
    samples        = INPUTS.out.samples
    refs           = INPUTS.out.refs
    targets        = INPUTS.out.targets
    zipped_targets = zipped_targets
}
