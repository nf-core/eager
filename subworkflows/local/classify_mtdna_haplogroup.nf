// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { addNewMetaFromAttributes      } from '../../subworkflows/local/utils_nfcore_eager_pipeline/main'
include { HAPLOGREP3_CLASSIFY           } from '../../modules/nf-core/haplogrep3/classify/main'

workflow CLASSIFY_MTDNA_HAPLOGROUP {

    take:
    ch_mtdna_vcf

    main:
    ch_versions      = Channel.empty()
    ch_haplogroups   = Channel.empty()

    ch_input_haplogrep3 = ch_mtdna_vcf

    HAPLOGREP3_CLASSIFY(ch_input_haplogrep3)
    ch_haplogroups = HAPLOGREP3_CLASSIFY.out.txt
    ch_versions    = ch_versions.mix(HAPLOGREP3_CLASSIFY.out.versions)

    emit:
    haplogroups    = ch_haplogroups
    versions       = ch_versions
}
