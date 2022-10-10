//
// Check input samplesheet and get read channel
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:

    if ( samplesheet.extension == 'tsv' ){
        ch_splitsamplesheet_for_branch = SAMPLESHEET_CHECK ( samplesheet )
            .csv
            .splitCsv ( header:true, sep:'\t' )
    } else {
        ch_splitsamplesheet_for_branch = SAMPLESHEET_CHECK ( samplesheet )
            .csv
            .splitCsv ( header:true )
    }

    ch_splitsamplesheet_for_branch
        .branch {
            bam: it.bam.toString().endsWith(".bam")
            fastq: it
        }.set { reads }

    fastqs = reads.fastq.map { create_fastq_channel(it) }
    bams = reads.bam.map { create_bam_channel(it) }

    emit:
    fastqs    // channel: [ val(meta), [ reads ] ]
    bams      // channel  [ val(mea), bam ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    // TODO create spanning main metadata
    meta.id                 = [ row.sample_id, row.library_id, "L" + row.lane ].join("_").trim()

    meta.sample_id          = row.sample_id
    meta.library_id         = row.library_id
    meta.lane               = row.lane
    meta.colour_chemistry   = row.colour_chemistry
    meta.single_end         = row.pairment == 'single'
    meta.strandedness       = row.strandedness
    meta.damage_treatment   = row.damage_treatment

    def array = []
    if (!file(row.r1).exists()) {
        exit 1, "[nf-core/eager] error: Please check input samplesheet. Read 1 FASTQ file does not exist! File: ${row.r1}"
    }
    if ( meta.single_end) {
        array = [ meta, [ file(row.r1) ] ]
    } else {
        if (!file(row.r2).exists()) {
            exit 1, "[nf-core/eager] error: Please check input samplesheet. Read 2 FASTQ file does not exist! File: ${row.r2}"
        }
        array = [ meta, [ file(row.r1), file(row.r2) ] ]
    }
    return array
}

def create_bam_channel(LinkedHashMap row) {
    def meta = [:]
    // TODO create spanning main metadata
    meta.id                 = [ row.sample_id, row.library_id ].join("_").trim()

    meta.sample_id          = row.sample_id
    meta.library_id         = row.library_id
    meta.strandedness       = row.strandedness
    meta.damage_treatment   = row.damage_treatment

    def array = []
    if (!file(row.bam).exists()) {
        exit 1, "[nf-core/eager] error: Please check input samplesheet. BAM file does not exist!\nFile: ${row.bam}"
    } else {
        array = [ meta, file(row.bam) ]
    }
    return array
}
