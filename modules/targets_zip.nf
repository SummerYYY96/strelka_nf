process TARGETS_ZIP {
    input:
    file targets_bed

    output:
    tuple(
        file("targets.bed.bgz"), file("targets.bed.bgz.tbi")
    )

    script:
    """
    sort -V -k1,1 -k2,2 "${targets_bed}" > targets.sorted.bed
    bgzip -c targets.sorted.bed > targets.bed.bgz
    tabix -p bed targets.bed.bgz
    """
}