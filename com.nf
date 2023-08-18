// Enable DSL2 functionality
nextflow.enable.dsl = 2

// Workflow definition
workflow {
	// 1) input a file with fasta file their corresponding reference ID
	fasta_info = Channel.fromPath(params.fasta_info).splitCsv( sep: "\t" ).map { row -> [file(row[0]), row[1]] }

	//fasta_info.view()

	// 3) the reads to be mapped
	// params.reads

	// 2) input bowtie2 index of the concacenated fasta files
	// params.fasta_index

	// 4) params.threads

	// params.label

	//workflow------------------

	//INDEX() // memory-intensive so not included so far

	MAPPING(params.fasta_index, params.reads)
	//MAPPING.out.view()

	GET_CONTIGS_ID_AND_FASTA(fasta_info)
	//GET_CONTIGS_ID_AND_FASTA.out.view()

	//combine contig ID file into one: contig ID, reference ID
	COMBINE_CONFIG_ID(GET_CONTIGS_ID_AND_FASTA.out.collect(), params.label)
	//COMBINE_CONFIG_ID.out.view()

	SORT_AND_INDEX_BAM(MAPPING.out)
	//SORT_AND_INDEX_BAM.out.view()

	// contig ID, contig length, number of mapped reads, number of unmapped reads
	IDXSTAT(SORT_AND_INDEX_BAM.out)
	//IDXSTAT.out.view()

	JOIN_MAP_WITH_SOURCE_INFO_OF_CONTIGS(IDXSTAT.out, COMBINE_CONFIG_ID.out, params.label)
	//JOIN_MAP_WITH_SOURCE_INFO_OF_CONTIGS.out.view()

	SUM_MAPPED_READS_TO_FASTA( JOIN_MAP_WITH_SOURCE_INFO_OF_CONTIGS.out )
	//SUM_MAPPED_READS_TO_FASTA.out.view()

	EXTRACT_BAM_OF_REF_ID(GET_CONTIGS_ID_AND_FASTA.out, SORT_AND_INDEX_BAM.out.map{it[0]}, params.label)

	//some dependency issue
	//AUTHENTICATION( EXTRACT_BAM_OF_REF_ID.out, params.label )

}

process MAPPING {
	publishDir params.mapping, 
		mode: "copy"
	 
	input:
		val(index)
		path(reads)

	output:
		path("*.bam")
	
	script:	  
	"""
	bowtie2 --very-sensitive -p ${params.threads} -x $index -U $reads | \
  samtools view -@ ${params.threads} -Sb -q 1 - > \$(basename $reads).bam
	"""
}

process GET_CONTIGS_ID_AND_FASTA {
        publishDir params.contigs,
                mode: "copy"

        input:
		tuple path(fasta), val(ref_id)

        output:
                path("*.contigs")

        script:
        """
	#extract all contigs

	if [[ $fasta == *.gz ]]; then
		zcat "$fasta" | LC_ALL=C fgrep ">" | sed 's/>//g' > "${ref_id}.contigs"
	else
		cat "$fasta" | LC_ALL=C fgrep ">" | sed 's/>//g' > "${ref_id}.contigs"
	fi

	#add a second column as the taxon, output *.contigs
	awk -v taxon="$ref_id" '{ \$2 = \$2 "\t" taxon } 1' "${ref_id}.contigs" > temp_file && mv temp_file "${ref_id}.contigs"

	#correct format------------------
	sed -i 's/ /\t/g' ${ref_id}.contigs

        # Count the number of columns in the file
        num_columns=\$(awk '{print NF}' ${ref_id}.contigs | sort -nu | tail -n 1)

        # If the file has more than 2 columns, use cut to extract only the first and third columns
        if [ "\$num_columns" -gt 2 ]; then
                cut -f 1,3 ${ref_id}.contigs > tmp
                mv tmp ${ref_id}.contigs
        fi

	sed -i 's/\t\t/\t/g' ${ref_id}.contigs
        """
}

process COMBINE_CONFIG_ID{

	publishDir params.contigs,
                mode: "copy"
        
        input:  
                path(config_info)
		val(label)
        
        output: 
                path("*all.contigs")
        
        script:
        """
        config_info_bash=\$(echo $config_info | sed 's/[][]//g')

	cat \$config_info_bash >> ${label}_all.contigs
        """

}

process SORT_AND_INDEX_BAM {
        publishDir params.mapping,
                mode: "copy"

        input:
                path(bam)

        output:
		tuple path("*_sorted.bam*"), env(bam_name)

        script:
        """
        bam_name=\$(echo $bam | cut -f 1 -d".")

        samtools sort -o \${bam_name}_sorted.bam ${bam}

	#use csi index to remove the limitation of the size of bam.
        samtools index -c \${bam_name}_sorted.bam
        """
}

process IDXSTAT {
        publishDir params.results,
                mode: "copy"

        input:
		tuple path(sorted_bam), val(bam_name)

        output:
                path("*.idxstats")

        script:
        """
        samtools idxstats ${bam_name}_sorted.bam > ${bam_name}.idxstats
        """
}

process JOIN_MAP_WITH_SOURCE_INFO_OF_CONTIGS {
        publishDir params.results,
                mode: "copy"

        input:
                path(idxstats)
		path(config_from)
		val(label)

        output:
                path("*.mapped_config_from")

        script:
        """
	join_contig_with_speciesname.py $idxstats $config_from ${label}.mapped_config_from
        """
}

process SUM_MAPPED_READS_TO_FASTA {
        publishDir params.results,
                mode: "copy"

        input:
                path(mapped_config_from)

        output:
                path("*_reads_sum.csv")

        script:
        """
	filename=\$(echo "$mapped_config_from" | sed 's/.mapped_config_from//')
	echo \$filename
	contig_info_to_sum.py \$filename
        """
}

process EXTRACT_BAM_OF_REF_ID {
	// errorStrategy { task.exitStatus == '137' ? 'retry' : 'terminate' } 
	//memory { 6.GB * task.attempt }

	errorStrategy 'ignore'


        publishDir params.mapping,
                mode: "copy"

        input:
		path(contig_from)
		tuple path(sorted_bam), path(sorted_bam_index)
		val(label)
		
        output:
                path("*-ext.bam")

        script:
        """
	#extract all contig IDs
	cut -f 1 $contig_from > contig_IDs

	#get the reference ID
	ref_id=\$(cut -f 2 $contig_from | head -1)

	#output the header
	samtools view -H $sorted_bam >> ${label}-\${ref_id}-ext.sam

	samtools view $sorted_bam | LC_ALL=C fgrep -f contig_IDs >> ${label}-\${ref_id}-ext.sam

	#sam to bam
	samtools view -S -b ${label}-\${ref_id}-ext.sam > ${label}-\${ref_id}-ext.bam
        """
}


process AUTHENTICATION {

	errorStrategy 'ignore'

	cpus 10

        publishDir params.plot,
                mode: "copy"

        input:
		path(bam)
		val(label)
        output:
                path("*.pdf")

        script:
        """
	#extract reference ID from the bam name
	ref_id=\$(echo $bam | cut -d'-' -f 2)

	bam_name=\$(echo $bam | sed 's/.bam//')
	sorted_bam=\${bam_name}_sorted.bam
	samtools sort -o \$sorted_bam \${bam_name}.bam
	samtools index -c \$sorted_bam

	echo -e "\${ref_id}\\t\${sorted_bam}" > bam_file
	AMBER --bamfiles bam_file --out "${label}-\${ref_id}-damage"
        """
}
