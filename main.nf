#!/usr/bin/env nextflow

/*
* Define workflow
*/
workflow {

// read in the sample sheet:
    Channel
        .fromPath(params.input)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [id: row.sample, condition: row.condition]
            [meta, file(row.path, checkIfExists: true), file(row.path + ".bai", checkIfExists: true)]
        }
        .set { ch_bam }

// collect meta data & bam/bai files:
    ch_bam
        .collect { meta, bam, bai -> [meta, bam, bai] }
        .map { items ->
            def triples = items.collate(3)
            def metas   = triples.collect { it[0] }
            def bams    = triples.collect { it[1] }
            def bais    = triples.collect { it[2] }
            [metas, bams, bais]
        }
        .set { ch_all_bams }

// actual steps:
// 1 - deeptools, make counts/bin txt file
    MULTIBAMSUMMARY(ch_all_bams)

// 2 - normalize counts/bin txt file by BPM
    NORMALIZE_BINS(
        MULTIBAMSUMMARY.out.bed.map { metas, bed -> bed },
        MULTIBAMSUMMARY.out.bed.map { metas, bed -> metas }
    )

// workflow now splits: step 3 & 4 should run in parallel
    
// Step 3: make one chromDict per sample in parallel
    NORMALIZE_BINS.out.metas
        .flatten()
        .combine(NORMALIZE_BINS.out.bpm)
        | R_CHROMDICT_SAMPLE

// Step 4: make one chromDict per condition in sequence
    R_CHROMDICT_CONDITION(
        NORMALIZE_BINS.out.bpm,
        NORMALIZE_BINS.out.metas
    )
}

/*
* Workflow step 1
*/
process MULTIBAMSUMMARY {
    tag "all_samples"
    label 'process_high'
    container 'quay.io/biocontainers/deeptools:3.5.4--pyhdfd78af_1'

    publishDir "${params.outputPath}", mode: 'copy', pattern: "bed/*"

    input:
    tuple val(metas), path(bams), path(bais)

    output:
    tuple val(metas), path("bed/bam_summary_rawCounts.bed"), emit: bed  // carry metas forward

    script:
    def bam_files     = bams.join(' ')
    def blacklist_arg = params.blacklist ? "--blackListFileName ${params.blacklist}" : ""

    """
    mkdir -p bed

    multiBamSummary bins \\
        --bamfiles ${bam_files} \\
        -p ${task.cpus} \\
        --outRawCounts bed/bam_summary_rawCounts.bed \\
        -bs ${params.binSize} \\
        --minMappingQuality ${params.mapQuality} \\
        --ignoreDuplicates \\
        ${blacklist_arg}
    """
}


/*
* Workflow step 2
*/
process NORMALIZE_BINS {
    tag "all_samples"
    label 'process_medium'
    container 'docker://rocker/tidyverse:4.4.1'

    publishDir "${params.outputPath}/bed", mode: 'copy'

    input:
    path bed
    val metas

    output:
    path "bam_summary_BPM.tsv", emit: bpm
    val metas, emit: metas

    script:
    def samples = metas.collect { it.id }.join(',')

    """
    Rscript ${projectDir}/bin/normalize_bins.R \\
        "${bed}" \\
        "${samples}" \\
        "bam_summary_BPM.tsv"
    """
}

/*
* Workflow step 3
*/
process R_CHROMDICT_SAMPLE {
    tag "${meta.id}"
    label 'process_medium'
    container 'docker://rocker/tidyverse:4.4.1'

    publishDir "${params.outputPath}/perSample_chromDicts", mode: 'copy'

    input:
    tuple val(meta), path(bpm)

    output:
    path "*.rds", emit: per_sample

    script:
    """
    Rscript ${projectDir}/bin/chromdict_sample.R \\
        "${meta.id}" \\
        "${bpm}"
    """
}

/*
* Workflow step 4
*/
process R_CHROMDICT_CONDITION {
    tag "all_samples"
    label 'process_medium'
    container 'docker://rocker/tidyverse:4.4.1'

    publishDir "${params.outputPath}/perCondition_chromDicts", mode: 'copy'

    input:
    path bpm
    val  metas

    output:
    path "*.rds", emit: per_condition

    script:
    def samples = metas.collect { it.id }.join(',')
    def conditions = metas.collect { it.condition }.join(',')

    """
    Rscript ${projectDir}/bin/chromdict_condition.R \\
        "${samples}" \\
        "${conditions}" \\
        "${bpm}"
    """
}

