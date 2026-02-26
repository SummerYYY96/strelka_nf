nextflow.enable.dsl = 2

include { SOMATIC_PAIRS } from './workflows/analyze_somatic_pairs.nf'

workflow {
    SOMATIC_PAIRS()
}
