#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    "nf-core/eager": ["v_pipeline.txt", r"(\S+)"],
    "Nextflow": ["v_nextflow.txt", r"(\S+)"],
    "FastQC": ["v_fastqc.txt", r"FastQC v(\S+)"],
    "MultiQC": ["v_multiqc.txt", r"multiqc, version (\S+)"],
    'AdapterRemoval':['v_adapterremoval.txt', r"AdapterRemoval ver. (\S+)"],
    'Picard MarkDuplicates': ['v_markduplicates.txt', r"(\S+)"],
    'Samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'Preseq': ['v_preseq.txt', r"Version: (\S+)"],
    'BWA': ['v_bwa.txt', r"Version: (\S+)"], 
    'Bowtie2': ['v_bowtie2.txt', r"bowtie2-([0-9]+\.[0-9]+\.[0-9]+) -fdebug"],
    'Qualimap': ['v_qualimap.txt', r"QualiMap v.(\S+)"],
    'GATK HaplotypeCaller': ['v_gatk.txt', r" v(\S+)"],
    #'GATK UnifiedGenotyper': ['v_gatk3_5.txt', r"version (\S+)"],
    'bamUtil' : ['v_bamutil.txt', r"Version: (\S+);"],
    'fastP': ['v_fastp.txt', r"([\d\.]+)"],
    'DamageProfiler' : ['v_damageprofiler.txt', r"DamageProfiler v(\S+)"],
    'angsd':['v_angsd.txt',r"version: (\S+)"],
    'bedtools':['v_bedtools.txt',r"bedtools v(\S+)"],
    'circulargenerator':['v_circulargenerator.txt',r"CircularGeneratorv(\S+)"],
    'DeDup':['v_dedup.txt',r"DeDup v(\S+)"],
    'freebayes':['v_freebayes.txt',r"v([0-9]\S+)"],
    'sequenceTools':['v_sequencetools.txt',r"(\S+)"],
    'maltextract':['v_maltextract.txt', r"version(\S+)"],
    'malt':['v_malt.txt',r"version (\S+)"],
    'multivcfanalyzer':['v_multivcfanalyzer.txt', r"MultiVCFAnalyzer - (\S+)"],
    'pmdtools':['v_pmdtools.txt',r"pmdtools v(\S+)"],
    'sexdeterrmine':['v_sexdeterrmine.txt',r"(\S+)"],
    'MTNucRatioCalculator':['v_mtnucratiocalculator.txt',r"Version: (\S+)"],
    'VCF2genome':['v_vcf2genome.txt', r"VCF2Genome \(v. ([0-9].[0-9]+) "],
    'endorS.py':['v_endorSpy.txt', r"endorS.py (\S+)"],
    'kraken':['v_kraken.txt', r"Kraken version (\S+)"],
    'eigenstrat_snp_coverage':['v_eigenstrat_snp_coverage.txt',r"(\S+)"],
    'mapDamage2':['v_mapdamage.txt',r"(\S+)"],
    'bbduk':['v_bbduk.txt',r"(\S+)"]
}

results = OrderedDict()
results["nf-core/eager"] = '<span style="color:#999999;">N/A</span>'
results["Nextflow"] = '<span style="color:#999999;">N/A</span>'
results["FastQC"] = '<span style="color:#999999;">N/A</span>'
results["MultiQC"] = '<span style="color:#999999;">N/A</span>'
results['AdapterRemoval'] = '<span style="color:#999999;\">N/A</span>'
results['fastP'] = '<span style="color:#999999;\">N/A</span>'
results['BWA'] = '<span style="color:#999999;\">N/A</span>'
results['Bowtie2'] = '<span style="color:#999999;\">N/A</span>'
results['circulargenerator'] = '<span style="color:#999999;\">N/A</span>'
results['Samtools'] = '<span style="color:#999999;\">N/A</span>'
results['endorS.py'] = '<span style="color:#999999;\">N/A</span>'
results['DeDup'] = '<span style="color:#999999;\">N/A</span>'
results['Picard MarkDuplicates'] = '<span style="color:#999999;\">N/A</span>'
results['Qualimap'] = '<span style="color:#999999;\">N/A</span>'
results['Preseq'] = '<span style="color:#999999;\">N/A</span>'
results['GATK HaplotypeCaller'] = '<span style="color:#999999;\">N/A</span>'
results['GATK UnifiedGenotyper'] = '<span style="color:#999999;\">N/A</span>'
results['freebayes'] = '<span style="color:#999999;\">N/A</span>'
results['sequenceTools'] = '<span style="color:#999999;\">N/A</span>'
results['VCF2genome'] = '<span style="color:#999999;\">N/A</span>'
results['MTNucRatioCalculator'] = '<span style="color:#999999;\">N/A</span>'
results['bedtools'] = '<span style="color:#999999;\">N/A</span>'
results['DamageProfiler'] = '<span style="color:#999999;\">N/A</span>'
results['bamUtil'] = '<span style="color:#999999;\">N/A</span>'
results['pmdtools'] = '<span style="color:#999999;\">N/A</span>'
results['angsd'] = '<span style="color:#999999;\">N/A</span>'
results['sexdeterrmine'] = '<span style="color:#999999;\">N/A</span>'
results['multivcfanalyzer'] = '<span style="color:#999999;\">N/A</span>'
results['malt'] = '<span style="color:#999999;\">N/A</span>'
results['kraken'] = '<span style="color:#999999;\">N/A</span>'
results['maltextract'] = '<span style="color:#999999;\">N/A</span>'
results['eigenstrat_snp_coverage'] = '<span style="color:#999999;\">N/A</span>'
results['mapDamage2'] = '<span style="color:#999999;\">N/A</span>'
results['bbduk'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del results[k]

# Dump to YAML
print(
    """
id: 'software_versions'
section_name: 'nf-core/eager Software Versions'
section_href: 'https://github.com/nf-core/eager'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
"""
)
for k, v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k, v))
print("    </dl>")

# Write out regexes as csv file:
with open("software_versions.csv", "w") as f:
    for k, v in results.items():
        f.write("{}\t{}\n".format(k, v))
