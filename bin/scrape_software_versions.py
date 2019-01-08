#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'nf-core/eager': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'Picard MarkDuplicates': ['v_markduplicates.txt', r"([\d\.]+)-SNAPSHOT"],
    'Samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'Preseq': ['v_preseq.txt', r"Version: (\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
    'BWA': ['v_bwa.txt', r"Version: (\S+)"],
    'Qualimap': ['v_qualimap.txt', r"QualiMap v.(\S+)"],
    'GATK': ['v_gatk.txt', r"Version:([\d\.]+)"],
    'bamUtil' : ['v_bamutil.txt', r"Version: ([\d\.]+)"],
    'fastP': ['v_fastp.txt', r"([\d\.]+)"],
    'DamageProfiler' : ['v_damageprofiler.txt', r"version\": \"([\d\.]+)"],
}
results = OrderedDict()
results['nf-core/eager'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'
results['Picard MarkDuplicates'] = '<span style="color:#999999;\">N/A</span>'
results['Samtools'] = '<span style="color:#999999;\">N/A</span>'
results['Preseq'] = '<span style="color:#999999;\">N/A</span>'
results['BWA'] = '<span style="color:#999999;\">N/A</span>'
results['Qualimap'] = '<span style="color:#999999;\">N/A</span>'
results['GATK'] = '<span style="color:#999999;\">N/A</span>'
results['bamUtil'] = '<span style="color:#999999;\">N/A</span>'
results['fastP'] = '<span style="color:#999999;\">N/A</span>'
results['DamageProfiler'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = "v{}".format(match.group(1))

# Dump to YAML
print ('''
id: 'nf-core/eager-software-versions'
section_name: 'nf-core/eager Software Versions'
section_href: 'https://github.com/nf-core/eager'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd>{}</dd>".format(k,v))
print ("    </dl>")
