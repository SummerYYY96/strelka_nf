workflow INPUTS {

    main:
    samples = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.Sample,
                row.Tumor,
                row.Normal,
                file(row.Tumor_Bam),
                file(row.Tumor_Bai),
                file(row.Normal_Bam),
                file(row.Normal_Bai)
            )
        }

    refs = Channel.value([
        file(params.ref_fa),
        file(params.ref_fai),
        file(params.ref_dict)
    ])

    targets = Channel.fromPath(params.targetbed)

    emit:
    samples
    refs
    targets
}
