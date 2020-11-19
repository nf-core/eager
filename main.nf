#!/usr/bin/env nextflow
/*

    // You wouldn't want to make UDG treated reads even shorter, so skip trimming if UDG.
    // We assume same trim amount for both non-UDG/UDG half as could trim a bit more off half-UDG to match non-UDG if needed, with minimal effect
    // Note: Trimming of e.g. adapters are sequencing artefacts and should be removed before mapping, so we don't account for this here.
    ch_bamutils_decision = ch_rmdup_for_bamutils.branch{
        totrim: it[6] == 'none' || it[6] == 'half'
        notrim: it[6] == 'full'
    }

} else {

    ch_bamutils_decision = ch_rmdup_for_bamutils.branch{
        totrim: it[6] == "dummy"
        notrim: it[6] == 'full' || it[6] == 'none' || it[6] == 'half'
    }

}

process bam_trim {
    label 'mc_small'
    tag "${libraryid}"
    publishDir "${params.outdir}/trimmed_bam", mode: params.publish_dir_mode

    when: params.run_trim_bam

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(bam), path(bai) from ch_bamutils_decision.totrim

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.trimmed.bam"), file("*.trimmed.bam.{bai,csi}") into ch_trimmed_from_bamutils

    script:
    def softclip = params.bamutils_softclip ? '-c' : ''
    def size = params.large_ref ? '-c' : ''
    def left_clipping = udg == "half" ? "${params.bamutils_clip_half_udg_left}" : "${params.bamutils_clip_none_udg_left}"
    def right_clipping = udg == "half" ? "${params.bamutils_clip_half_udg_right}" : "${params.bamutils_clip_none_udg_right}"
    """
    bam trimBam $bam tmp.bam -L ${left_clipping} -R ${right_clipping} ${softclip}
    samtools sort -@ ${task.cpus} tmp.bam -o ${libraryid}.trimmed.bam
    samtools index ${libraryid}.trimmed.bam ${size}
    """
}

// Post trimming merging of libraries to single samples, except for SS/DS
// libraries as they should be genotyped separately, because we will assume
// that if trimming is turned on, 'lab-removed' libraries can be combined with
// merged with 'in-silico damage removed' libraries to improve genotyping

ch_trimmed_formerge = ch_bamutils_decision.notrim
  .mix(ch_trimmed_from_bamutils)
  .groupTuple(by:[0,4,5])
  .map{
        def samplename = it[0]
        def libraryid  = it[1]
        def lane = it[2]
        def seqtype = it[3]
        def organism = it[4]
        def strandedness = it[5]
        def udg = it[6]
        def bam = it[7].flatten()
        def bai = it[8].flatten()

      [samplename, libraryid, lane, seqtype, organism, strandedness, udg, bam, bai ]
  }
  .branch{
    skip_merging: it[7].size() == 1
    merge_me: it[7].size() > 1
  }

//////////////////////////////////////////////////////////////////////////
/* --    POST aDNA BAM MODIFICATION AND FINAL MAPPING STATISTICS     -- */
//////////////////////////////////////////////////////////////////////////

process additional_library_merge {
  label 'sc_tiny'
  tag "${samplename}"
  publishDir "${params.outdir}/merged_bams/additional", mode: params.publish_dir_mode

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(bam), path(bai) from ch_trimmed_formerge.merge_me

  output:
  tuple samplename, val("${samplename}_libmerged"), lane, seqtype, organism, strandedness, udg, path("*_libmerged_rg_add.bam"), path("*_libmerged_rg_add.bam.{bai,csi}") into ch_output_from_trimmerge

  script:
  def size = params.large_ref ? '-c' : ''
  """
  samtools merge ${samplename}_libmerged_add.bam ${bam}
  picard AddOrReplaceReadGroups I=${samplename}_libmerged_add.bam O=${samplename}_libmerged_rg_add.bam RGID=1 RGLB="${samplename}_additionalmerged" RGPL=illumina RGPU=4410 RGSM="${samplename}_additionalmerged" VALIDATION_STRINGENCY=LENIENT
  samtools index ${samplename}_libmerged_rg_add.bam ${size}
  """
}

ch_trimmed_formerge.skip_merging
  .mix(ch_output_from_trimmerge)
  .into{ ch_output_from_bamutils; ch_addlibmerge_for_qualimap; ch_for_sexdeterrmine }

  // General mapping quality statistics for whole reference sequence - e.g. X and % coverage

process qualimap {
    label 'mc_small'
    tag "${samplename}"
    publishDir "${params.outdir}/qualimap", mode: params.publish_dir_mode

    when:
    !params.skip_qualimap

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(bam), path(bai) from ch_addlibmerge_for_qualimap
    file fasta from ch_fasta_for_qualimap.collect()

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*") into ch_qualimap_results

    script:
    def snpcap = params.snpcapture_bed != '' ? "-gff ${params.snpcapture_bed}" : ''
    """
    qualimap bamqc -bam $bam -nt ${task.cpus} -outdir . -outformat "HTML" ${snpcap}
    """
}

/////////////////////////////
/* --    GENOTYPING     -- */
/////////////////////////////

// Reroute files for genotyping; we have to ensure to select lib-merged BAMs, as input channel will also contain the un-merged ones resulting in unwanted multi-sample VCFs
if ( params.run_genotyping && params.genotyping_source == 'raw' ) {
    ch_output_from_bamutils
      .into { ch_damagemanipulation_for_skipgenotyping; ch_damagemanipulation_for_genotyping_ug; ch_damagemanipulation_for_genotyping_hc; ch_damagemanipulation_for_genotyping_freebayes; ch_damagemanipulation_for_genotyping_pileupcaller; ch_damagemanipulation_for_genotyping_angsd }

} else if ( params.run_genotyping && params.genotyping_source == "trimmed" && !params.run_trim_bam )  {
    exit 1, "[nf-core/eager] error: Cannot run genotyping with 'trimmed' source without running BAM trimming (--run_trim_bam)! Please check input parameters."

} else if ( params.run_genotyping && params.genotyping_source == "trimmed" && params.run_trim_bam )  {
    ch_output_from_bamutils
        .into { ch_damagemanipulation_for_skipgenotyping; ch_damagemanipulation_for_genotyping_ug; ch_damagemanipulation_for_genotyping_hc; ch_damagemanipulation_for_genotyping_freebayes; ch_damagemanipulation_for_genotyping_pileupcaller; ch_damagemanipulation_for_genotyping_angsd }

} else if ( params.run_genotyping && params.genotyping_source == "pmd" && !params.run_pmdtools )  {
    exit 1, "[nf-core/eager] error: Cannot run genotyping with 'pmd' source without running pmtools (--run_pmdtools)! Please check input parameters."

} else if ( params.run_genotyping && params.genotyping_source == "pmd" && params.run_pmdtools )  {
  ch_output_from_pmdtools
    .into { ch_damagemanipulation_for_skipgenotyping; ch_damagemanipulation_for_genotyping_ug; ch_damagemanipulation_for_genotyping_hc; ch_damagemanipulation_for_genotyping_freebayes; ch_damagemanipulation_for_genotyping_pileupcaller; ch_damagemanipulation_for_genotyping_angsd }

} else if ( !params.run_genotyping && !params.run_trim_bam && !params.run_pmdtools )  {
    ch_rmdup_for_skipdamagemanipulation
    .into { ch_damagemanipulation_for_skipgenotyping; ch_damagemanipulation_for_genotyping_ug; ch_damagemanipulation_for_genotyping_hc; ch_damagemanipulation_for_genotyping_freebayes; ch_damagemanipulation_for_genotyping_pileupcaller; ch_damagemanipulation_for_genotyping_angsd }

} else if ( !params.run_genotyping && !params.run_trim_bam && params.run_pmdtools )  {
    ch_rmdup_for_skipdamagemanipulation
    .into { ch_damagemanipulation_for_skipgenotyping; ch_damagemanipulation_for_genotyping_ug; ch_damagemanipulation_for_genotyping_hc; ch_damagemanipulation_for_genotyping_freebayes; ch_damagemanipulation_for_genotyping_pileupcaller; ch_damagemanipulation_for_genotyping_angsd }

} else if ( !params.run_genotyping && params.run_trim_bam && !params.run_pmdtools )  {
    ch_rmdup_for_skipdamagemanipulation
    .into { ch_damagemanipulation_for_skipgenotyping; ch_damagemanipulation_for_genotyping_ug; ch_damagemanipulation_for_genotyping_hc; ch_damagemanipulation_for_genotyping_freebayes; ch_damagemanipulation_for_genotyping_pileupcaller; ch_damagemanipulation_for_genotyping_angsd }

}

// Unified Genotyper - although not-supported, better for aDNA (because HC does de novo assembly which requires higher coverages), and needed for MultiVCFAnalyzer

process genotyping_ug {
  label 'mc_small'
  tag "${samplename}"
  publishDir "${params.outdir}/genotyping", mode: params.publish_dir_mode, pattern: '*{.vcf.gz,.realign.bam,realign.bai}'

  when:
  params.run_genotyping && params.genotyping_tool == 'ug'

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_damagemanipulation_for_genotyping_ug
  file fasta from ch_fasta_for_genotyping_ug.collect()
  file fai from ch_fai_for_ug.collect()
  file dict from ch_dict_for_ug.collect()

  output:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*vcf.gz") into ch_ug_for_multivcfanalyzer,ch_ug_for_vcf2genome
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.realign.{bam,bai}") optional true

  script:
  def defaultbasequalities = params.gatk_ug_defaultbasequalities == '' ? '' : " --defaultBaseQualities ${params.gatk_ug_defaultbasequalities}"
  def keep_realign = params.gatk_ug_keep_realign_bam ? "samtools index ${samplename}.realign.bam" : "rm ${samplename}.realign.{bam,bai}"
  if (params.gatk_dbsnp == '')
    """
    samtools index -b ${bam}
    gatk3 -T RealignerTargetCreator -R ${fasta} -I ${bam} -nt ${task.cpus} -o ${samplename}.intervals ${defaultbasequalities}
    gatk3 -T IndelRealigner -R ${fasta} -I ${bam} -targetIntervals ${samplename}.intervals -o ${samplename}.realign.bam ${defaultbasequalities}
    gatk3 -T UnifiedGenotyper -R ${fasta} -I ${samplename}.realign.bam -o ${samplename}.unifiedgenotyper.vcf -nt ${task.cpus} --genotype_likelihoods_model ${params.gatk_ug_genotype_model} -stand_call_conf ${params.gatk_call_conf} --sample_ploidy ${params.gatk_ploidy} -dcov ${params.gatk_downsample} --output_mode ${params.gatk_ug_out_mode} ${defaultbasequalities}

    $keep_realign

    pigz -p ${task.cpus} ${samplename}.unifiedgenotyper.vcf
    """
  else if (params.gatk_dbsnp != '')
    """
    samtools index ${bam}
    gatk3 -T RealignerTargetCreator -R ${fasta} -I ${bam} -nt ${task.cpus} -o ${samplename}.intervals ${defaultbasequalities}
    gatk3 -T IndelRealigner -R ${fasta} -I ${bam} -targetIntervals ${samplenane}.intervals -o ${samplename}.realign.bam ${defaultbasequalities}
    gatk3 -T UnifiedGenotyper -R ${fasta} -I ${samplename}.realign.bam -o ${samplename}.unifiedgenotyper.vcf -nt ${task.cpus} --dbsnp ${params.gatk_dbsnp} --genotype_likelihoods_model ${params.gatk_ug_genotype_model} -stand_call_conf ${params.gatk_call_conf} --sample_ploidy ${params.gatk_ploidy} -dcov ${params.gatk_downsample} --output_mode ${params.gatk_ug_out_mode} ${defaultbasequalities}

    $keep_realign

    pigz -p ${task.cpus} ${samplename}.unifiedgenotyper.vcf
    """
}

 // HaplotypeCaller as 'best practise' tool for human DNA in particular

process genotyping_hc {
  label 'mc_small'
  tag "${samplename}"
  publishDir "${params.outdir}/genotyping", mode: params.publish_dir_mode

  when:
  params.run_genotyping && params.genotyping_tool == 'hc'

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_damagemanipulation_for_genotyping_hc
  file fasta from ch_fasta_for_genotyping_hc.collect()
  file fai from ch_fai_for_hc.collect()
  file dict from ch_dict_for_hc.collect()

  output:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*vcf.gz")

  script:
  if (params.gatk_dbsnp == '')
    """
    gatk HaplotypeCaller -R ${fasta} -I ${bam} -O ${samplename}.haplotypecaller.vcf -stand-call-conf ${params.gatk_call_conf} --sample-ploidy ${params.gatk_ploidy} --output-mode ${params.gatk_hc_out_mode} --emit-ref-confidence ${params.gatk_hc_emitrefconf}
    pigz -p ${task.cpus} ${samplename}.haplotypecaller.vcf
    """

  else if (params.gatk_dbsnp != '')
    """
    gatk HaplotypeCaller -R ${fasta} -I ${bam} -O ${samplename}.haplotypecaller.vcf --dbsnp ${params.gatk_dbsnp} -stand-call-conf ${params.gatk_call_conf} --sample_ploidy ${params.gatk_ploidy} --output_mode ${params.gatk_hc_out_mode} --emit-ref-confidence ${params.gatk_hc_emitrefconf}
    pigz -p ${task.cpus} ${samplename}.haplotypecaller.vcf
    """
}

 // Freebayes for 'more efficient/simple' and more generic genotyping (vs HC)

process genotyping_freebayes {
  label 'mc_small'
  tag "${samplename}"
  publishDir "${params.outdir}/genotyping", mode: params.publish_dir_mode

  when:
  params.run_genotyping && params.genotyping_tool == 'freebayes'

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_damagemanipulation_for_genotyping_freebayes
  file fasta from ch_fasta_for_genotyping_freebayes.collect()
  file fai from ch_fai_for_freebayes.collect()
  file dict from ch_dict_for_freebayes.collect()

  output:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*vcf.gz")

  script:
  def skip_coverage = "${params.freebayes_g}" == 0 ? "" : "-g ${params.freebayes_g}"
  """
  freebayes -f ${fasta} -p ${params.freebayes_p} -C ${params.freebayes_C} ${skip_coverage} ${bam} > ${samplename}.freebayes.vcf
  pigz -p ${task.cpus} ${samplename}.freebayes.vcf
  """
}


 // Branch channel by strandedness
ch_damagemanipulation_for_genotyping_pileupcaller
  .branch{
      singleStranded: it[5] == "single"
      doubleStranded: it[5] == "double"
  }
  .set{ch_input_for_genotyping_pileupcaller}

 // Create pileupcaller input tuples
ch_input_for_genotyping_pileupcaller.singleStranded
  .groupTuple(by:[5])
  .set {ch_prepped_for_pileupcaller_single}

ch_input_for_genotyping_pileupcaller.doubleStranded
  .groupTuple(by:[5])
  .set {ch_prepped_for_pileupcaller_double}

process genotyping_pileupcaller {
  label 'mc_small'
  tag "${strandedness}"
  publishDir "${params.outdir}/genotyping", mode: params.publish_dir_mode

  when:
  params.run_genotyping && params.genotyping_tool == 'pileupcaller'

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, bam, bai from ch_prepped_for_pileupcaller_double.mix(ch_prepped_for_pileupcaller_single)
  file fasta from ch_fasta_for_genotyping_pileupcaller.collect()
  file fai from ch_fai_for_pileupcaller.collect()
  file dict from ch_dict_for_pileupcaller.collect()
  path(bed) from ch_bed_for_pileupcaller.collect()
  path(snp) from ch_snp_for_pileupcaller.collect().dump(tag: "Pileupcaller SNP file")

  output:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("pileupcaller.${strandedness}.*") into ch_for_eigenstrat_snp_coverage

  script:
  def use_bed = bed.getName() != 'nf-core_eager_dummy.txt' ? "-l ${bed}" : ''
  def use_snp = snp.getName() != 'nf-core_eager_dummy2.txt' ? "-f ${snp}" : ''

  def transitions_mode = strandedness == "single" ? "" : "${params.pileupcaller_transitions_mode}" == 'SkipTransitions' ? "--skipTransitions" : "${params.pileupcaller_transitions_mode}" == 'TransitionsMissing' ? "--transitionsMissing" : ""
  def caller = "--${params.pileupcaller_method}"
  def ssmode = strandedness == "single" ? "--singleStrandMode" : ""
  def bam_list = bam.flatten().join(" ")
  def sample_names = samplename.flatten().join(",")
  """
  samtools mpileup -B -q 30 -Q 30 ${use_bed} -f ${fasta} ${bam_list} | pileupCaller ${caller} ${ssmode} ${transitions_mode} --sampleNames ${sample_names} ${use_snp} -e pileupcaller.${strandedness}
  """
 }

process eigenstrat_snp_coverage {
   label 'mc_tiny'
   tag "${strandedness}"
   publishDir "${params.outdir}/genotyping", mode: params.publish_dir_mode

   when:
   params.run_genotyping && params.genotyping_tool == 'pileupcaller'

   input:
   tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*") from ch_for_eigenstrat_snp_coverage.dump()

   output:
   tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*.json") into ch_eigenstrat_snp_cov_for_multiqc
   path("*_eigenstrat_coverage.txt")

   script:
   /*
   The following code block can be swapped in once the eigenstratdatabasetools MultiQC module becomes available.
   """
   eigenstrat_snp_coverage -i pileupcaller.${strandedness} -s ".txt" >${strandedness}_eigenstrat_coverage.txt -j ${strandedness}_eigenstrat_coverage_mqc.json
   """
   */
   """
   eigenstrat_snp_coverage -i pileupcaller.${strandedness} -s ".txt" >${strandedness}_eigenstrat_coverage.txt
   parse_snp_cov.py ${strandedness}_eigenstrat_coverage.txt
   """
 }

process genotyping_angsd {
  label 'mc_small'
  tag "${samplename}"
  publishDir "${params.outdir}/genotyping", mode: params.publish_dir_mode

  when:
  params.run_genotyping && params.genotyping_tool == 'angsd'

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_damagemanipulation_for_genotyping_angsd
  file fasta from ch_fasta_for_genotyping_angsd.collect()
  file fai from ch_fai_for_angsd.collect()
  file dict from ch_dict_for_angsd.collect()

  output:
  path("${samplename}*")

  script:
  switch ( "${params.angsd_glmodel}" ) {
    case "samtools":
    angsd_glmodel = "1"; break
    case "gatk":
    angsd_glmodel = "2"; break
    case "soapsnp":
    angsd_glmodel = "3"; break
    case "syk":
    angsd_glmodel = "4"; break
  }

  switch ( "${params.angsd_glformat}" ) {
    case "text":
    angsd_glformat = "4"; break
    case "binary":
    angsd_glformat = "1"; break
    case "beagle":
    angsd_glformat = "2"; break
    case "binary_three":
    angsd_glformat = "3"; break
  }

  def angsd_fasta = !params.angsd_createfasta ? '' : params.angsd_fastamethod == 'random' ? '-doFasta 1 -doCounts 1' : '-doFasta 2 -doCounts 1'
  """
  echo ${bam} > bam.filelist
  mkdir angsd
  angsd -bam bam.filelist -nThreads ${task.cpus} -GL ${angsd_glmodel} -doGlF ${angsd_glformat} ${angsd_fasta} -out ${samplename}.angsd
  """
}

////////////////////////////////////
/* --    CONSENSUS CALLING     -- */
////////////////////////////////////

// Generate a simple consensus-called FASTA file based on genotype VCF

process vcf2genome {
  label  'mc_small'
  tag "${samplename}"
  publishDir "${params.outdir}/consensus_sequence", mode: params.publish_dir_mode

  when:
  params.run_vcf2genome

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(vcf) from ch_ug_for_vcf2genome
  file fasta from ch_fasta_for_vcf2genome.collect()

  output:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*.fasta.gz")

  script:
  def out = "${params.vcf2genome_outfile}" == '' ? "${samplename}.fasta" : "${params.vcf2genome_outfile}"
  def fasta_head = "${params.vcf2genome_header}" == '' ? "${samplename}" : "${params.vcf2genome_header}"
  """
  pigz -f -d -p ${task.cpus} *.vcf.gz
  vcf2genome -draft ${out}.fasta -draftname "${fasta_head}" -in ${vcf.baseName} -minc ${params.vcf2genome_minc} -minfreq ${params.vcf2genome_minfreq} -minq ${params.vcf2genome_minq} -ref ${fasta} -refMod ${out}_refmod.fasta -uncertain ${out}_uncertainy.fasta
  pigz -p ${task.cpus} *.fasta
  pigz -p ${task.cpus} *.vcf
  """
}

// More complex consensus caller with additional filtering functionality (e.g. for heterozygous calls) to generate SNP tables and other things sometimes used in aDNA bacteria studies

// Create input channel for MultiVCFAnalyzer, possibly mixing with pre-made VCFs.
if (params.additional_vcf_files == '') {
    ch_vcfs_for_multivcfanalyzer = ch_ug_for_multivcfanalyzer.map{ it[7] }.collect()
} else {
    ch_vcfs_for_multivcfanalyzer = ch_ug_for_multivcfanalyzer.map{ it [7] }.collect().mix(ch_extravcfs_for_multivcfanalyzer)
}

process multivcfanalyzer {
  label  'mc_small'
  publishDir "${params.outdir}/multivcfanalyzer", mode: params.publish_dir_mode

  when:
  params.genotyping_tool == 'ug' && params.run_multivcfanalyzer && params.gatk_ploidy == '2'

  input:
  file vcf from ch_vcfs_for_multivcfanalyzer.collect()
  file fasta from ch_fasta_for_multivcfanalyzer.collect()

  output:
  file('fullAlignment.fasta.gz')
  file('info.txt.gz')
  file('snpAlignment.fasta.gz')
  file('snpAlignmentIncludingRefGenome.fasta.gz')
  file('snpStatistics.tsv.gz')
  file('snpTable.tsv.gz')
  file('snpTableForSnpEff.tsv.gz')
  file('snpTableWithUncertaintyCalls.tsv.gz')
  file('structureGenotypes.tsv.gz')
  file('structureGenotypes_noMissingData-Columns.tsv.gz')
  file('MultiVCFAnalyzer.json') optional true into ch_multivcfanalyzer_for_multiqc

  script:
  def write_freqs = params.write_allele_frequencies ? "T" : "F"
  """
  gunzip -f *.vcf.gz
  multivcfanalyzer ${params.snp_eff_results} ${fasta} ${params.reference_gff_annotations} . ${write_freqs} ${params.min_genotype_quality} ${params.min_base_coverage} ${params.min_allele_freq_hom} ${params.min_allele_freq_het} ${params.reference_gff_exclude} *.vcf
  pigz -p ${task.cpus} *.tsv *.txt snpAlignment.fasta snpAlignmentIncludingRefGenome.fasta fullAlignment.fasta
  """
 }

////////////////////////////////////////////////////////////
/* --    HUMAN DNA SPECIFIC ADDITIONAL INFORMATION     -- */
////////////////////////////////////////////////////////////

// Mitochondrial to nuclear ratio helps to evaluate quality of tissue sampled

 process mtnucratio {
  tag "${samplename}"
  publishDir "${params.outdir}/mtnucratio", mode: params.publish_dir_mode

  when:
  params.run_mtnucratio

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(bam), path(bai) from ch_rmdup_formtnucratio

  output:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*.mtnucratio")
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*.json") into ch_mtnucratio_for_multiqc

  script:
  """
  mtnucratio ${bam} "${params.mtnucratio_header}"
  """
 }

// Human biological sex estimation

// As we collect all files for a single sex_deterrmine run, we DO NOT use the normal input/output tuple
process sex_deterrmine {
    label 'sc_small'
    publishDir "${params.outdir}/sex_determination", mode: params.publish_dir_mode

    input:
    path bam from ch_for_sexdeterrmine.map { it[7] }.collect()
    path(bed) from ch_bed_for_sexdeterrmine

    output:
    file "SexDet.txt"
    file "*.json" into ch_sexdet_for_multiqc

    when:
    params.run_sexdeterrmine

    script:
    def filter = bed.getName() != 'nf-core_eager_dummy.txt' ? "-b $bed" : ''
    """

    for i in *.bam; do
        echo \$i >> bamlist.txt
    done

    samtools depth -aa -q30 -Q30 $filter -f bamlist.txt | sexdeterrmine -f bamlist.txt > SexDet.txt
    """
}

// Human DNA nuclear contamination estimation

 process nuclear_contamination{
    label 'sc_small'
    tag "${samplename}"
    publishDir "${params.outdir}/nuclear_contamination", mode: params.publish_dir_mode

    when:
    params.run_nuclear_contamination

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(input), path(bai) from ch_for_nuclear_contamination

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path('*.X.contamination.out') into ch_from_nuclear_contamination

    script:
    """
    samtools index ${input}
    angsd -i ${input} -r ${params.contamination_chrom_name}:5000000-154900000 -doCounts 1 -iCounts 1 -minMapQ 30 -minQ 30 -out ${libraryid}.doCounts
    contamination -a ${libraryid}.doCounts.icnts.gz -h ${baseDir}/assets/angsd_resources/HapMapChrX.gz 2> ${libraryid}.X.contamination.out
    """
 }

// As we collect all files for a single print_nuclear_contamination run, we DO NOT use the normal input/output tuple
process print_nuclear_contamination{
    label 'sc_tiny'
    publishDir "${params.outdir}/nuclear_contamination", mode: params.publish_dir_mode

    when:
    params.run_nuclear_contamination

    input:
    val 'Contam' from ch_from_nuclear_contamination.map { it[7] }.collect()

    output:
    file 'nuclear_contamination.txt'
    file 'nuclear_contamination_mqc.json' into ch_nuclear_contamination_for_multiqc

    script:
    """
    print_x_contamination.py ${Contam.join(' ')}
    """
 }

/////////////////////////////////////////////////////////
/* --    METAGENOMICS-SPECIFIC ADDITIONAL STEPS     -- */
/////////////////////////////////////////////////////////

// MALT is a super-fast BLAST replacement typically used for pathogen detection or microbiome profiling against large databases, here using off-target reads from mapping

// As we collect all files for a all metagenomic runs, we DO NOT use the normal input/output tuple!
if (params.metagenomic_tool == 'malt') {
  ch_bam_filtering_for_metagenomic
    .set {ch_bam_filtering_for_metagenomic_malt}

  ch_bam_filtering_for_metagenomic_kraken = Channel.empty()
} else if (params.metagenomic_tool == 'kraken') {
  ch_bam_filtering_for_metagenomic
    .set {ch_bam_filtering_for_metagenomic_kraken}

  ch_bam_filtering_for_metagenomic_malt = Channel.empty()
} else if ( params.metagenomic_tool == '' ) {
  ch_bam_filtering_for_metagenomic_malt = Channel.empty()
  ch_bam_filtering_for_metagenomic_kraken = Channel.empty()

}

// As we collect all files for a single MALT run, we DO NOT use the normal input/output tuple
process malt {
  label 'mc_small'
  publishDir "${params.outdir}/metagenomic_classification/malt", mode: params.publish_dir_mode

  when:
  params.run_metagenomic_screening && params.run_bam_filtering && params.bam_unmapped_type == 'fastq' && params.metagenomic_tool == 'malt'

  input:
  file fastqs from ch_bam_filtering_for_metagenomic_malt.map { it[7] }.collect()
  file db from ch_db_for_malt

  output:
  path("*.rma6") into ch_rma_for_maltExtract
  path("*.sam.gz") optional true
  path("malt.log") into ch_malt_for_multiqc

  script:
  if ( "${params.malt_min_support_mode}" == "percent" ) {
    min_supp = "-supp ${params.malt_min_support_percent}"
  } else if ( "${params.malt_min_support_mode}" == "reads" ) {
    min_supp = "-sup ${params.metagenomic_min_support_reads}"
  }
  def sam_out = params.malt_sam_output ? "-a . -f SAM" : ""
  """
  malt-run \
  -J-Xmx${task.memory.toGiga()}g \
  -t ${task.cpus} \
  -v \
  -o . \
  -d ${db} \
  ${sam_out} \
  -id ${params.percent_identity} \
  -m ${params.malt_mode} \
  -at ${params.malt_alignment_mode} \
  -top ${params.malt_top_percent} \
  ${min_supp} \
  -mq ${params.malt_max_queries} \
  --memoryMode ${params.malt_memory_mode} \
  -i ${fastqs.join(' ')} |&tee malt.log
  """
}

// MaltExtract performs aDNA evaluation from the output of MALT (damage patterns, read lengths etc.)

// As we collect all files for a single MALT extract run, we DO NOT use the normal input/output tuple
process maltextract {
  label 'mc_medium'
  publishDir "${params.outdir}/maltextract/", mode: params.publish_dir_mode

  when:
  params.run_maltextract && params.metagenomic_tool == 'malt'

  input:
  file rma6 from ch_rma_for_maltExtract.collect()
  file taxon_list from ch_taxonlist_for_maltextract
  file ncbifiles from ch_ncbifiles_for_maltextract

  output:
  path "results/" type('dir')
  file "results/*_Wevid.json" optional true into ch_hops_for_multiqc

  script:
  def destack = params.maltextract_destackingoff ? "--destackingOff" : ""
  def downsam = params.maltextract_downsamplingoff ? "--downSampOff" : ""
  def dupremo = params.maltextract_duplicateremovaloff ? "--dupRemOff" : ""
  def matches = params.maltextract_matches ? "--matches" : ""
  def megsum = params.maltextract_megansummary ? "--meganSummary" : ""
  def topaln = params.maltextract_topalignment ?  "--useTopAlignment" : ""
  def ss = params.single_stranded ? "--singleStranded" : ""
  """
  MaltExtract \
  -Xmx${task.memory.toGiga()}g \
  -t ${taxon_list} \
  -i ${rma6.join(' ')} \
  -o results/ \
  -r ${ncbifiles} \
  -p ${task.cpus} \
  -f ${params.maltextract_filter} \
  -a ${params.maltextract_toppercent} \
  --minPI ${params.maltextract_percentidentity} \
  ${destack} \
  ${downsam} \
  ${dupremo} \
  ${matches} \
  ${megsum} \
  ${topaln} \
  ${ss}

  postprocessing.AMPS.r -r results/ -m ${params.maltextract_filter} -t ${task.cpus} -n ${taxon_list} -j
  """
}

// Kraken is offered as a replacement for MALT as MALT is _very_ resource hungry

if (params.run_metagenomic_screening && params.database.endsWith(".tar.gz") && params.metagenomic_tool == 'kraken'){
  comp_kraken = file(params.database)

  process decomp_kraken {
    input:
    path(ckdb) from comp_kraken

    output:
    path(dbname) into ch_krakendb

    script:
    dbname = params.database.tokenize("/")[-1].tokenize(".")[0]
    """
    tar xvzf $ckdb
    """
  }

} else if (! params.database.endsWith(".tar.gz") && params.run_metagenomic_screening && params.metagenomic_tool == 'kraken') {
    ch_krakendb = path(params.database)
} else {
    ch_krakendb = Channel.empty()
}

process kraken {
  tag "$prefix"
  label 'mc_huge'
  publishDir "${params.outdir}/metagenomic_classification/kraken", mode: params.publish_dir_mode

  when:
  params.run_metagenomic_screening && params.run_bam_filtering && params.bam_unmapped_type == 'fastq' && params.metagenomic_tool == 'kraken'

  input:
  path(fastq) from ch_bam_filtering_for_metagenomic_kraken.map { it[7] }
  path(krakendb) from ch_krakendb

  output:
  file "*.kraken.out" into ch_kraken_out
  tuple prefix, path("*.kreport") into ch_kraken_report, ch_kraken_for_multiqc

  script:
  prefix = fastq.toString().tokenize('.')[0]
  out = prefix+".kraken.out"
  kreport = prefix+".kreport"

  """
  kraken2 --db ${krakendb} --threads ${task.cpus} --output $out --report $kreport $fastq
  """
}

process kraken_parse {
  tag "$name"
  errorStrategy 'ignore'

  input:
  tuple val(name), path(kraken_r) from ch_kraken_report

  output:
  tuple val(name), path('*.kraken_parsed.csv') into ch_kraken_parsed

  script:
  out = name+".kraken_parsed.csv"
  """
  kraken_parse.py -c ${params.metagenomic_min_support_reads} -o $out $kraken_r
  """
}

process kraken_merge {
  publishDir "${params.outdir}/metagenomic_classification/kraken", mode: params.publish_dir_mode

  input:
  file csv_count from ch_kraken_parsed.map{ it[1] }.collect()

  output:
  path('kraken_count_table.csv')

  script:
  out = "kraken_count_table.csv"
  """
  merge_kraken_res.py -o $out
  """
}

//////////////////////////////////////
/* --    PIPELINE COMPLETION     -- */
//////////////////////////////////////

// Pipeline documentation for on-server guidance

process output_documentation {
    label 'sc_tiny'
    publishDir "${params.outdir}/documentation", mode: params.publish_dir_mode

    input:
    file output_docs from ch_output_docs
    file images from ch_output_docs_images

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

// Collect all software versions for inclusion in MultiQC report

process get_software_versions {
  label 'sc_tiny'
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version &> v_pipeline.txt
    echo $workflow.nextflow.version &> v_nextflow.txt

    fastqc --version &> v_fastqc.txt 2>&1 || true
    AdapterRemoval --version  &> v_adapterremoval.txt 2>&1 || true
    fastp --version &> v_fastp.txt 2>&1 || true
    bwa &> v_bwa.txt 2>&1 || true
    circulargenerator --help | head -n 1 &> v_circulargenerator.txt 2>&1 || true
    samtools --version &> v_samtools.txt 2>&1 || true
    dedup -v &> v_dedup.txt 2>&1 || true
    ## bioconda recipe of picard is incorrectly set up and extra warning made with stderr, this ugly command ensures only version exported
    ( exec 7>&1; picard MarkDuplicates --version 2>&1 >&7 | grep -v '/' >&2 ) 2> v_markduplicates.txt || true
    qualimap --version &> v_qualimap.txt 2>&1 || true
    preseq &> v_preseq.txt 2>&1 || true
    gatk --version 2>&1 | head -n 1 > v_gatk.txt 2>&1 || true
    freebayes --version &> v_freebayes.txt 2>&1 || true
    bedtools --version &> v_bedtools.txt 2>&1 || true
    damageprofiler --version &> v_damageprofiler.txt 2>&1 || true
    bam --version &> v_bamutil.txt 2>&1 || true
    pmdtools --version &> v_pmdtools.txt 2>&1 || true
    angsd -h |& head -n 1 | cut -d ' ' -f3-4 &> v_angsd.txt 2>&1 || true
    multivcfanalyzer --help | head -n 1 &> v_multivcfanalyzer.txt 2>&1 || true
    malt-run --help |& tail -n 3 | head -n 1 | cut -f 2 -d'(' | cut -f 1 -d ',' &> v_malt.txt 2>&1 || true
    MaltExtract --help | head -n 2 | tail -n 1 &> v_maltextract.txt 2>&1 || true
    multiqc --version &> v_multiqc.txt 2>&1 || true
    vcf2genome -h |& head -n 1 &> v_vcf2genome.txt || true
    mtnucratio --help &> v_mtnucratiocalculator.txt || true
    sexdeterrmine --version &> v_sexdeterrmine.txt || true
    kraken2 --version | head -n 1 &> v_kraken.txt || true
    endorS.py --version &> v_endorSpy.txt || true
    pileupCaller --version &> v_sequencetools.txt 2>&1 || true
    bowtie2 --version | grep -a 'bowtie2-.* -fdebug' > v_bowtie2.txt || true
    eigenstrat_snp_coverage --version | cut -d ' ' -f2 >v_eigenstrat_snp_coverage.txt || true

    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

// MultiQC file generation for pipeline report

process multiqc {
    label 'sc_medium'

    publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode

    input:
    file multiqc_config from ch_multiqc_config
    file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
    file software_versions_mqc from software_versions_yaml.collect().ifEmpty([])
    file logo from ch_eager_logo
    file ('fastqc_raw/*') from ch_prefastqc_for_multiqc.collect().ifEmpty([])
    path('fastqc/*') from ch_fastqc_after_clipping.collect().ifEmpty([])
    file ('adapter_removal/*') from ch_adapterremoval_logs.collect().ifEmpty([])
    file ('mapping/bt2/*') from ch_bt2_for_multiqc.collect().ifEmpty([])
    file ('flagstat/*') from ch_flagstat_for_multiqc.collect().ifEmpty([])
    file ('flagstat_filtered/*') from ch_bam_filtered_flagstat_for_multiqc.collect().ifEmpty([])
    file ('preseq/*') from ch_preseq_for_multiqc.collect().ifEmpty([])
    file ('damageprofiler/dmgprof*/*') from ch_damageprofiler_results.collect().ifEmpty([])
    file ('qualimap/qualimap*/*') from ch_qualimap_results.collect().ifEmpty([])
    file ('markdup/*') from ch_markdup_results_for_multiqc.collect().ifEmpty([])
    file ('dedup*/*') from ch_dedup_results_for_multiqc.collect().ifEmpty([])
    file ('fastp/*') from ch_fastp_for_multiqc.collect().ifEmpty([])
    file ('sexdeterrmine/*') from ch_sexdet_for_multiqc.collect().ifEmpty([])
    file ('mutnucratio/*') from ch_mtnucratio_for_multiqc.collect().ifEmpty([])
    file ('endorspy/*') from ch_endorspy_for_multiqc.collect().ifEmpty([])
    file ('multivcfanalyzer/*') from ch_multivcfanalyzer_for_multiqc.collect().ifEmpty([])
    file ('malt/*') from ch_malt_for_multiqc.collect().ifEmpty([])
    file ('kraken/*') from ch_kraken_for_multiqc.collect().ifEmpty([])
    file ('hops/*') from ch_hops_for_multiqc.collect().ifEmpty([])
    file ('nuclear_contamination/*') from ch_nuclear_contamination_for_multiqc.collect().ifEmpty([])
    file ('genotyping/*') from ch_eigenstrat_snp_cov_for_multiqc.collect().ifEmpty([])

    file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

    output:
    file "*multiqc_report.html" into ch_multiqc_report
    file "*_data"

    script:
    def rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    def rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    def custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc -f $rtitle $rfilename $multiqc_config $custom_config_file .
    """
}

// Send completion emails if requested, so user knows data is ready

workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/eager] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/eager] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report instanceof ArrayList) {
                log.warn "[nf-core/eager] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/eager] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/eager] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if (mqc_report == NULL) {
                log.warn "[nf-core/eager] Could not attach MultiQC report to summary email"
            } else if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
                mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/eager] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/eager]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/eager]${c_red} Pipeline completed with errors${c_reset}-"
    }
}

/////////////////////////////////////
/* --    AUXILARY FUNCTIONS     -- */
/////////////////////////////////////

def nfcoreHeader() {
    // Log colours ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/eager v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}

// Channelling the TSV file containing FASTQ or BAM
def extract_data(tsvFile) {
    Channel.fromPath(tsvFile)
        .splitCsv(header: true, sep: '\t')
        .dump()
        .map { row ->

            def expected_keys = ['Sample_Name', 'Library_ID', 'Lane', 'Colour_Chemistry', 'SeqType', 'Organism', 'Strandedness', 'UDG_Treatment', 'R1', 'R2', 'BAM']
            if ( !row.keySet().containsAll(expected_keys) ) exit 1, "[nf-core/eager] error: Invalid TSV input - malformed column names. Please check input TSV. Column names should be: ${expected_keys.join(", ")}"

            checkNumberOfItem(row, 11)

            def samplename = row.Sample_Name
            def libraryid  = row.Library_ID
            def lane = row.Lane
            def colour = row.Colour_Chemistry
            def seqtype = row.SeqType
            def organism = row.Organism
            def strandedness = row.Strandedness
            def udg = row.UDG_Treatment
            def r1 = row.R1.matches('NA') ? 'NA' : return_file(row.R1)
            def r2 = row.R2.matches('NA') ? 'NA' : return_file(row.R2)
            def bam = row.BAM.matches('NA') ? 'NA' : return_file(row.BAM)

            // check no empty metadata fields
            if (samplename == '' || libraryid == '' || lane == '' || colour == '' || seqtype == '' || seqtype == '' || udg == '' || r1 == '' || r2 == '') exit 1, "[nf-core/eager] error: a field does not contain any information. Ensure all cells are filled or contain 'NA' for optional fields. Check row:\n ${row}"

            // Check no 'empty' rows
            if (r1.matches('NA') && r2.matches('NA') && bam.matches('NA') && bai.matches('NA')) exit 1, "[nf-core/eager] error: A row in your TSV appears to have all files defined as NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check row for: ${samplename}"

            // Ensure BAMs aren't submitted with PE
            if (!bam.matches('NA') && seqtype.matches('PE')) exit 1, "[nf-core/eager] error: BAM input rows in TSV cannot be set as PE, only SE. See '--help' flag and documentation under 'running the pipeline' for more information. Check row for: ${samplename}"

            // Check valid UDG treatment
            if (!udg.matches('none') && !udg.matches('half') && !udg.matches('full')) exit 1, "[nf-core/eager] error: UDG treatment can only be 'none', 'half' or 'full'. See '--help' flag and documentation under 'running the pipeline' for more information. You have '${udg}'"

            // Check valid colour chemistry
            if (!colour == 2 && !colour == 4) exit 1, "[nf-core/eager] error: Colour chemistry in TSV can either be 2 (e.g. NextSeq/NovaSeq) or 4 (e.g. HiSeq/MiSeq)"

            //  Ensure that we do not accept incompatible chemistry setup
            if (!seqtype.matches('PE') && !seqtype.matches('SE')) exit 1, "[nf-core/eager] error:  SeqType for one or more rows in TSV is neither SE nor PE! see '--help' flag and documentation under 'running the pipeline' for more information. You have: '${seqtype}'"

           // So we don't accept existing files that are wrong format: e.g. fasta or sam
            if ( !r1.matches('NA') && !has_extension(r1, "fastq.gz") && !has_extension(r1, "fq.gz") && !has_extension(r1, "fastq") && !has_extension(r1, "fq")) exit 1, "[nf-core/eager] error: A specified R1 file either has a non-recognizable FASTQ extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${r1}"
            if ( !r2.matches('NA') && !has_extension(r2, "fastq.gz") && !has_extension(r2, "fq.gz") && !has_extension(r2, "fastq") && !has_extension(r2, "fq")) exit 1, "[nf-core/eager] error: A specified R2 file either has a non-recognizable FASTQ extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${r2}"
            if ( !bam.matches('NA') && !has_extension(bam, "bam")) exit 1, "[nf-core/eager] error: A specified R1 file either has a non-recognizable BAM extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${bam}"

            [ samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, r1, r2, bam ]

        }

    }

// Check if a row has the expected number of item
def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "[nf-core/eager] error:  Invalid TSV input - malformed row (e.g. missing column) in ${row}, see '--help' flag and documentation under 'running the pipeline' for more information"
    return true
}

// Return file if it exists
def return_file(it) {
    if (!file(it).exists()) exit 1, "[nf-core/eager] error: Cannot find supplied FASTQ or BAM input file. If using input method TSV set to NA if no file required. See '--help' flag and documentation under 'running the pipeline' for more information. Check file: ${it}"
    return file(it)
}

// Check file extension
def has_extension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Extract FastQs from Path
// Create a channel of FASTQs from a directory pattern: "my_samples/*/"
// All FASTQ files in subdirectories are collected and emitted;
// they must have _R1_ and/or _R2_ in their names.
def retrieve_input_paths(input, colour_chem, pe_se, ds_ss, udg_treat, bam_in) {

  if ( !bam_in ) {
        if( pe_se ) {
            log.info "Generating single-end FASTQ data TSV"
            Channel
                .fromFilePairs( input, size: 1 )
                .filter { it =~/.*.fastq.gz|.*.fq.gz|.*.fastq|.*.fq/ }
                .ifEmpty { exit 1, "[nf-core/eager] error:  Your specified FASTQ read files did not end in: '.fastq.gz', '.fq.gz', '.fastq', or '.fq'. Did you forget --bam?" }
                .map { row -> [ row[0], [ row[1][0] ] ] }
                .ifEmpty { exit 1, "[nf-core/eager] error:  --input was empty - no input files supplied!" }
                .into { ch_reads_for_faketsv; ch_reads_for_validate }

                // Check we don't have any duplicated sample names due to fromFilePairs behaviour of calculating sample name from anything before R1/R2 glob
                ch_reads_for_validate
                  .groupTuple()
                  .map{
                    if ( validate_size(it[1], 1) ) { null } else { exit 1, "[nf-core/eager] error: You have supplied non-unique sample names (text before R1/R2 indication). Did you accidentally supply paired-end data?  see '--help' flag and documentation under 'running the pipeline' for more information. Check duplicates of: ${it[0]}" }
                  }

        } else if (!pe_se ){
            log.info "Generating paired-end FASTQ data TSV"

            Channel
                .fromFilePairs( input )
                .filter { it =~/.*.fastq.gz|.*.fq.gz|.*.fastq|.*.fq/ }
                .ifEmpty { exit 1, "[nf-core/eager] error: Files could not be found. Do the specified FASTQ read files end in: '.fastq.gz', '.fq.gz', '.fastq', or '.fq'? Did you forget --single_end?" }
                .map { row -> [ row[0], [ row[1][0], row[1][1] ] ] }
                .ifEmpty { exit 1, "[nf-core/eager] error: --input was empty - no input files supplied!" }
                .into { ch_reads_for_faketsv; ch_reads_for_validate }

                // Check we don't have any duplicated sample names due to fromFilePairs behaviour of calculating sample name from anything before R1/R2 glob
                ch_reads_for_validate
                  .groupTuple()
                  .map{
                    if ( validate_size(it[1], 1) ) { null } else { exit 1, "[nf-core/eager] error: You have supplied non-unique sample names (text before R1/R2 indication). See '--help' flag and documentation under 'running the pipeline' for more information. Check duplicates of: ${it[0]}" }
                  }

        }

    } else if ( bam_in ) {
              log.info "Generating BAM data TSV"

         Channel
            .fromFilePairs( input, size: 1 )
            .filter { it =~/.*.bam/ }
            .map { row -> [ row[0], [ row[1][0] ] ] }
            .ifEmpty { exit 1, "[nf-core/eager] error: Cannot find any bam file matching: ${input}" }
            .set { ch_reads_for_faketsv }

    }

ch_reads_for_faketsv
  .map{

      def samplename = it[0]
      def libraryid  = it[0]
      def lane = 0
      def colour = "${colour_chem}"
      def seqtype = pe_se ? 'SE' : 'PE'
      def organism = 'NA'
      def strandedness = ds_ss ? 'single' : 'double'
      def udg = udg_treat
      def r1 = !bam_in ? return_file(it[1][0]) : 'NA'
      def r2 = !bam_in && !pe_se ? return_file(it[1][1]) : 'NA'
      def bam = bam_in && pe_se ? return_file(it[1][0]) : 'NA'

      [ samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, r1, r2, bam ]
  }
  .ifEmpty {exit 1, "[nf-core/eager] error: Invalid file paths with --input"}

}

// Function to check length of collection in a channel closure is as expected (e.g. with .map())
def validate_size(collection, size){
    if ( collection.size() != size ) { return false } else { return true }
}
