/*
   nextflow \
    /ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/aertslab/scripts/src_dwmax/ngs/scrna-seq/drop-seq/nemesh/nextflow/main.nf \
        --outdir '.' \
        --reads '00.raw/*_S1_L001_R{1,2}_001.fastq.gz' \
        --selected_barcodes '../data/*.selected_cells.txt' \
        --selected_barcodes_tag 'paper' \
        -with-report report.html \
        -with-trace
 */

nextflow.preview.dsl=2

FASTP_DIR='/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/software/genius/fastp/0.20.0/bin'
FASTP_ADAPTERS="/vsc-hard-mounts/leuven-data/software/biomed/skylake_centos7/old_skylake_modules/fastp/0.20.0_6ff0ffa-foss-2018a/fastp.adapters"

/*
 * Create a channel for input read files
 */
Channel
    .fromFilePairs( params.reads, size: 2)
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .set { data }

// Check if genome exists in the config file
if (!params.selected_barcodes) {
    exit 1, "The argument selected_barcodes should be provided."
}
/*
 * Create a channel for input read files
 */
Channel
   .fromPath(params.selected_barcodes)
   .map {
      path -> tuple(path.baseName.split('\\.')[0], params.selected_barcodes_tag, path)
   }
   .set { selectedBarcodesByCustom }

data.subscribe { println it }
selectedBarcodesByCustom.subscribe { println it }

/**
 * Preprocess + FastQC
 */
process FASTP__FASTQ_CLEAN() {
    publishDir "${params.outdir}/01.clean", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        set val(sample), file(reads)
    output:
        tuple file('*_R{1,2}.clean.fastq.gz'), emit: fastq
        tuple file('*_fastp.{json,html}'), emit: report
    script:
        """
        ${FASTP_DIR}/fastp -w 6 \
            -i ${reads[0]} \
            -I ${reads[1]} \
            -o ${sample}_R1.clean.fastq.gz \
            -O ${sample}_R2.clean.fastq.gz \
            --length_required 20 \
            --adapter_fasta ${FASTP_ADAPTERS} \
            -j ${sample}_fastp.json \
            -h ${sample}_fastp.html
        """
}

process PICARD__FASTQ_TO_BAM {
    publishDir "${params.outdir}/01.clean", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        file(reads)
    output:
        tuple val(sample), file('*.unaligned.bam'), emit: bam
    script:
        sample = reads[0].toString() - ~/(_R1)?(\.clean)?(\.fq)?(\.fastq)?(\.gz)?$/
        """
        module load Java/1.8.0_192
        java -Djava.io.tmpdir=$DWMAX/tmp -jar \
            /ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/software/genius/picard/2.20.6/bin/picard.jar \
            FastqToSam \
                    FASTQ=${reads[0]} \
                    FASTQ2=${reads[1]} \
                    O=${sample}.unaligned.bam \
                    SAMPLE_NAME=${sample}
        """
}

params.baseRange = "1-12"
params.baseQuality = 10
params.barcodedRead = 1
params.discardRead = "False"
params.barcodeTagName = "XC"
params.numBasesBelowQuality = 1

process DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLBARCODE {
    publishDir "${params.outdir}/01.clean", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
        tuple val(sample), file('*.unaligned_tagged_Cell.bam'), emit: bam
        tuple file('*.unaligned_tagged_Cellular.bam_summary.txt'), emit: report
    script:
        """
		source $DWMAX/documents/aertslab/scripts/src_dwmax/bash-utils/utils.sh
		software load drop-seq_tools/1.12
		TagBamWithReadSequenceExtended \
			INPUT=${bam} \
			OUTPUT=${sample}.unaligned_tagged_Cell.bam \
			SUMMARY=${sample}.unaligned_tagged_Cellular.bam_summary.txt \
			BASE_RANGE=${params.baseRange} \
			BASE_QUALITY=${params.baseQuality} \
			BARCODED_READ=${params.barcodedRead} \
			DISCARD_READ=${params.discardRead} \
			TAG_NAME=${params.barcodeTagName} \
			NUM_BASES_BELOW_QUALITY=${params.numBasesBelowQuality}        
        """
}

params._baseRange = "13-20"
params._baseQuality = 10
params._barcodedRead = 1
params._discardRead = "True"
params._barcodeTagName = "XM"
params._numBasesBelowQuality = 1

process DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLMOLECULAR {
    publishDir "${params.outdir}/01.clean", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
        tuple val(sample), file('*.unaligned_tagged_CellMolecular.bam'), emit: bam
        tuple file('*.unaligned_tagged_Molecular.bam_summary.txt'), emit: report
    script:
        """
		source $DWMAX/documents/aertslab/scripts/src_dwmax/bash-utils/utils.sh
		software load drop-seq_tools/1.12
		TagBamWithReadSequenceExtended \
			INPUT=${bam} \
			OUTPUT=${sample}.unaligned_tagged_CellMolecular.bam \
			SUMMARY=${sample}.unaligned_tagged_Molecular.bam_summary.txt \
			BASE_RANGE=${params._baseRange} \
			BASE_QUALITY=${params._baseQuality} \
			BARCODED_READ=${params._barcodedRead} \
			DISCARD_READ=${params._discardRead} \
			TAG_NAME=${params._barcodeTagName} \
			NUM_BASES_BELOW_QUALITY=${params._numBasesBelowQuality}        
        """
}

params.tagReject = "XQ"

process DROP_SEQ_TOOLS__FILTER_UNALIGNED_TAGGED_BAM {
    publishDir "${params.outdir}/01.clean", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
        tuple val(sample), file('*.unaligned_tagged_filtered.bam'), emit: bam
    script:
        """
		source $DWMAX/documents/aertslab/scripts/src_dwmax/bash-utils/utils.sh
		software load drop-seq_tools/1.12
		FilterBAM \
			TAG_REJECT=${params.tagReject} \
			INPUT=${bam} \
			OUTPUT=${sample}.unaligned_tagged_filtered.bam
        """
}

params.adapterSequence = "AAGCAGTGGTATCAACGCAGAGTGAATGGG"
params.mismatches = 0
params.numBases = 5

process DROP_SEQ_TOOLS__TRIM_SMART_UNALIGNED_TAGGED_FILTERED_BAM {
    publishDir "${params.outdir}/01.clean", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
        tuple val(sample), file('*.unaligned_tagged_trimmed_smart.bam'), emit: bam
        tuple file('*.adapter_trimming_report.txt'), emit: report
    script:
        """
		source $DWMAX/documents/aertslab/scripts/src_dwmax/bash-utils/utils.sh
		software load drop-seq_tools/1.12
		TrimStartingSequence \
			INPUT=${bam} \
			OUTPUT=${sample}.unaligned_tagged_trimmed_smart.bam \
			OUTPUT_SUMMARY=${sample}.adapter_trimming_report.txt \
			SEQUENCE=${params.adapterSequence} \
			MISMATCHES=${params.mismatches} \
			NUM_BASES=${params.numBases}
        """
}

params._mismatches = 0
params._numBases = 6

process DROP_SEQ_TOOLS__TRIM_POLYA_UNALIGNED_TAGGED_TRIMMED_SMART {
    publishDir "${params.outdir}/01.clean", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
        tuple val(sample), file('*.unaligned_tagged_polyA_filtered.bam'), emit: bam
        tuple file('*.polyA_trimming_report.txt'), emit: report
    script:
        """
		source $DWMAX/documents/aertslab/scripts/src_dwmax/bash-utils/utils.sh
		software load drop-seq_tools/1.12
		PolyATrimmer \
			INPUT=${bam} \
			OUTPUT=${sample}.unaligned_tagged_polyA_filtered.bam \
			OUTPUT_SUMMARY=${sample}.polyA_trimming_report.txt \
			MISMATCHES=${params._mismatches} \
			NUM_BASES=${params._numBases}
        """
}

process PICARD__BAM_TO_FASTQ {
    publishDir "${params.outdir}/01.clean", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
        tuple val(sample), file('*.unaligned_tagged_polyA_filtered.fastq'), emit: fastq
    script:
        """
        module load Java/1.8.0_192
        java -Djava.io.tmpdir=$DWMAX/tmp -jar \
            /ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/software/genius/picard/2.20.6/bin/picard.jar \
                SamToFastq \
			        INPUT=${bam} \
			        FASTQ=${sample}.unaligned_tagged_polyA_filtered.fastq
        """    
}

process GZIP {

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(f)
    output:
        tuple val(sample), file("*.unaligned_tagged_polyA_filtered.fastq.gz"), emit: fastq_gz
    script:
        """
        gzip --force ${f}
        """
}

params.threads = 8
params.runMode = "genomeGenerate"
params.sjdbOverhang = (20+50)-1

process STAR__BUILD_INDEX {

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=1:00:00 -A ${params.qsubaccount}"

    input:
        file(annotation)
        file(genome)
    output:
        file("STAR_index")
    script:
        """
        mkdir STAR_index
        singularity run \
            -B /ddn1/vol1/staging/leuven/stg_00002/:/ddn1/vol1/staging/leuven/stg_00002/ \
            -B /staging/leuven/stg_00002/:/staging/leuven/stg_00002/ \
                /staging/leuven/res_00001/software/STAR/2.7.1a/STAR_2.7.1a.sif \
                --runThreadN ${params.threads} \
                --runMode genomeGenerate \
                --genomeDir STAR_index \
                --genomeFastaFiles ${genome} \
                --sjdbGTFfile ${annotation} \
                --sjdbOverhang ${params.sjdbOverhang} \
                --genomeSAindexNbases 13 # Suggested by STAR (default: 14), otherwise keeps on hanging
        """
}

process STAR__LOAD {
    input:
        file(f)
    output:
        file("STAR_LOADED")
    script:
        """
        singularity run \
            -B /ddn1/vol1/staging/leuven/stg_00002/:/ddn1/vol1/staging/leuven/stg_00002/ \
            -B /staging/leuven/stg_00002/:/staging/leuven/stg_00002/ \
                /staging/leuven/res_00001/software/STAR/2.7.1a/STAR_2.7.1a.sif \
                --genomeLoad LoadAndExit \
                --genomeDir $f
        rm -r _STARtmp Log.out Log.progress.out Aligned.out.sam
        touch STAR_LOADED
        """
}

params.STAR_threads = 8
params.STAR_limitBamSorRAM = 20000000000
params.STAR_outSamType = "BAM Unsorted"
// params.STAR_quantMode = "TranscriptomeSAM"
params.STAR_readFilesCommand = "zcat"

process STAR__ALIGN {
    publishDir "${params.outdir}/02.map", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(fastq)
        file(starIndex)
        file(starLoaded)
    output:
        tuple val(sample), file("*.STAR_Aligned.out.bam")
    script:
        """
        singularity run \
            -B /ddn1/vol1/staging/leuven/stg_00002/:/ddn1/vol1/staging/leuven/stg_00002/ \
            -B /staging/leuven/stg_00002/:/staging/leuven/stg_00002/ \
            /staging/leuven/res_00001/software/STAR/2.7.1a/STAR_2.7.1a.sif \
                --runThreadN ${params.STAR_threads} \
                --limitBAMsortRAM ${params.STAR_limitBamSorRAM} \
                --outSAMtype ${params.STAR_outSamType} \
                --readFilesCommand ${params.STAR_readFilesCommand} \
                --genomeDir ${starIndex} \
                --readFilesIn ${fastq} \
                --outFileNamePrefix ${sample}.STAR_
        """    
}

process PICARD__SORT_SAM {
    publishDir "${params.outdir}/02.map", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
        tuple val(sample), file("*.STAR_aligned_sorted.bam")
    script:
        """
        module load Java/1.8.0_192
        java -Djava.io.tmpdir=$DWMAX/tmp -jar \
            /ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/software/genius/picard/2.20.6/bin/picard.jar \
                SortSam \
                    I=${bam} \
                    O=${sample}.STAR_aligned_sorted.bam \
                    SO=coordinate
        """
}

process PICARD__CREATE_SEQUENCE_DICTIONARY {
    publishDir "${params.outdir}/00.refdata", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=1:00:00 -A ${params.qsubaccount}"

    input:
        file(genome)
    output:
        file "${genome.baseName}.dict"
    script:
        """
        module load Java/1.8.0_192
        java -Djava.io.tmpdir=$DWMAX/tmp -jar \
            /ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/software/genius/picard/2.20.6/bin/picard.jar \
                CreateSequenceDictionary \
                    R=${genome} \
                    O=${genome.baseName}.dict
        """
}

params.includeSecondaryAlignments = "false"
params.pairedRun = "false"

process PICARD__MERGE_BAM_ALIGNMENT {
    publishDir "${params.outdir}/02.map", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(unmappedBam)
        tuple val(sample), file(mappedBam)
        file(genome)
        file(dict)
    output:
        tuple val(sample), file("*.merged.bam")
    script:
        """
        module load Java/1.8.0_192
        java -Djava.io.tmpdir=$DWMAX/tmp -jar \
            /ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/software/genius/picard/2.20.6/bin/picard.jar \
                MergeBamAlignment \
                    REFERENCE_SEQUENCE=${genome} \
                    UNMAPPED_BAM=${unmappedBam} \
                    ALIGNED_BAM=${mappedBam} \
                    OUTPUT=${sample}.merged.bam \
                    INCLUDE_SECONDARY_ALIGNMENTS=${params.includeSecondaryAlignments} \
                    PAIRED_RUN=${params.pairedRun}
        """    
}

process FORMAT_GTF {
    publishDir "${params.outdir}/00.refdata", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=1:00:00 -A ${params.qsubaccount}"

    input:
        file(annotation)
    output:
        file "*.formatted.gtf"
    script:
        """
        sed -r 's/(.*); transcript_id (.*); (.*); gene_name (.*); \$/\\1; transcript_id \\2; \\3; gene_name \\4; transcript_name \\2;/' \
            ${annotation} \
            > ${annotation.baseName}.formatted.gtf
        """
}

process DROP_SEQ_TOOLS__CONVERT_TO_REFFLAT {
    publishDir "${params.outdir}/00.refdata", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=1:00:00 -A ${params.qsubaccount}"

    input:
        file(annotation)
        file(seqdict)
    output:
        file("${seqdict.baseName}.refFlat")
    script:
        """
		source $DWMAX/documents/aertslab/scripts/src_dwmax/bash-utils/utils.sh
		software load drop-seq_tools/1.12
		ConvertToRefFlat \
			ANNOTATIONS_FILE=${annotation} \
            SEQUENCE_DICTIONARY=${seqdict} \
            OUTPUT=${seqdict.baseName}.refFlat
        """
}

params.tag = "GE"

process DROP_SEQ_TOOLS__TAG_READ_WITH_GENE_EXON {
    publishDir "${params.outdir}/02.map", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
        file(annotation)
    output:
        tuple val(sample), file("*.merged_gene-exon-tagged.bam")
    script:
        """
		source $DWMAX/documents/aertslab/scripts/src_dwmax/bash-utils/utils.sh
		software load drop-seq_tools/1.12
		TagReadWithGeneExon \
			I=${bam} \
			O=${sample}.merged_gene-exon-tagged.bam \
			ANNOTATIONS_FILE=${annotation} \
			TAG=${params.tag}
        """    
} 

params.numBarcodes = 2000
params.primerSequence = "AAGCAGTGGTATCAACGCAGAGTAC"

process DROP_SEQ_TOOLS__DETECT_REPAIR_BARCODE_SYNTHESIS_ERRORS {
    publishDir "${params.outdir}/02.map", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
		tuple val(sample), file("*.final_cleaned.bam"), emit: bam
		tuple file("*.synthesis_stats.txt"), emit: stats
		// tuple file("*.synthesis_stats.summary.txt"), emit: statsSummary
    script:
        """
		source $DWMAX/documents/aertslab/scripts/src_dwmax/bash-utils/utils.sh
		software load drop-seq_tools/1.12
		DetectBeadSynthesisErrors \
			I=${bam} \
			O=${sample}.final_cleaned.bam \
			OUTPUT_STATS=${sample}.synthesis_stats.txt \
			SUMMARY=${sample}.synthesis_stats.summary.txt \
			NUM_BARCODES=${params.numBarcodes * 2} \
			PRIMER_SEQUENCE=${params.primerSequence} \
            TMP_DIR=$DWMAX/tmp
        """    
}

params.BAMTagHistogram_tag = "XC"

process DROP_SEQ_TOOLS__BAM_TAG_HISTOGRAM {
    publishDir "${params.outdir}/03.count", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=24:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam)
    output:
		tuple val(sample), file("*.cell_readcounts.txt.gz")
    script:
        """
		source $DWMAX/documents/aertslab/scripts/src_dwmax/bash-utils/utils.sh
		software load drop-seq_tools/1.12
		BAMTagHistogram \
			I=${bam} \
			O=${sample}.cell_readcounts.txt.gz \
			TAG=${params.BAMTagHistogram_tag}
        """
}

process DROPLET_UTILS__BARCODE_SELECTION {
    publishDir "03.count", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=1:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(readCounts)
    output:
		tuple val(sample), val("knee"), file("*.selected_cell_barcodes_by_knee.txt"), emit: selectedCellBarcodesByKnee
        tuple val(sample), val("inflection"), file("*.selected_cell_barcodes_by_inflection.txt"), emit: selectedCellBarcodesByInflection
        tuple file("*.barcode_rank_vs_total_umi_plot.png"), emit: plot
    script:
        """
        Rscript $DWMAX/documents/aertslab/scripts/src_dwmax/ngs/scrna-seq/drop-seq/nemesh/nextflow/droplet_utils_barcode_selection.R \
            ${readCounts} \
            ${sample}
        """    
}

process DROP_SEQ_TOOLS__DIGITAL_EXPRESSION {
    publishDir "03.count", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=1:00:00 -A ${params.qsubaccount}"

    input:
        tuple val(sample), file(bam), val(tag), file(selectedBarcodes)
    output:
        tuple file("*.${tag}.cells_dge.txt.gz"), emit: dgem
    shell:
        """
        source $DWMAX/documents/aertslab/scripts/src_dwmax/bash-utils/utils.sh
		software load drop-seq_tools/1.12
        DigitalExpression \
            I=${bam} \
            O=${sample}.${tag}.cells_dge.txt.gz \
            SUMMARY=${sample}.${tag}.cells_dge.summary.txt \
            CELL_BC_FILE=${selectedBarcodes} \
            TMP_DIR=$DWMAX/tmp
        """
}

params.genome = '/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/resources/refs/flybase/r6.16/dmel-all-chromosome-r6.16.fasta'
params.annotation = '/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/resources/refs/flybase/r6.16/dmel-all-r6.16.gtf'

workflow {
    main:
        FASTP__FASTQ_CLEAN( data )
        PICARD__FASTQ_TO_BAM( FASTP__FASTQ_CLEAN.out.fastq )
        DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLBARCODE( PICARD__FASTQ_TO_BAM.out.bam )
        DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLMOLECULAR( DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLBARCODE.out.bam )
        DROP_SEQ_TOOLS__FILTER_UNALIGNED_TAGGED_BAM( DROP_SEQ_TOOLS__TAG_UNALIGNED_BAM_WITH_CELLMOLECULAR.out.bam )
        DROP_SEQ_TOOLS__TRIM_SMART_UNALIGNED_TAGGED_FILTERED_BAM( DROP_SEQ_TOOLS__FILTER_UNALIGNED_TAGGED_BAM.out.bam )
        DROP_SEQ_TOOLS__TRIM_POLYA_UNALIGNED_TAGGED_TRIMMED_SMART( DROP_SEQ_TOOLS__TRIM_SMART_UNALIGNED_TAGGED_FILTERED_BAM.out.bam )
        PICARD__BAM_TO_FASTQ( DROP_SEQ_TOOLS__TRIM_POLYA_UNALIGNED_TAGGED_TRIMMED_SMART.out.bam )
        GZIP( PICARD__BAM_TO_FASTQ.out.fastq )
        STAR__BUILD_INDEX( file(params.annotation), file(params.genome) )
        STAR__LOAD( STAR__BUILD_INDEX.out )
        STAR__ALIGN( 
            GZIP.out.fastq_gz,
            STAR__BUILD_INDEX.out,
            STAR__LOAD.out
        )
        PICARD__SORT_SAM( STAR__ALIGN.out )
        PICARD__CREATE_SEQUENCE_DICTIONARY( file(params.genome) )
        PICARD__MERGE_BAM_ALIGNMENT( 
            DROP_SEQ_TOOLS__TRIM_POLYA_UNALIGNED_TAGGED_TRIMMED_SMART.out.bam,
            PICARD__SORT_SAM.out,
            file( params.genome ),
            PICARD__CREATE_SEQUENCE_DICTIONARY.out
        )
        FORMAT_GTF( file(params.annotation) )
        DROP_SEQ_TOOLS__CONVERT_TO_REFFLAT( 
            FORMAT_GTF.out, 
            PICARD__CREATE_SEQUENCE_DICTIONARY.out
        )
        DROP_SEQ_TOOLS__TAG_READ_WITH_GENE_EXON(
            PICARD__MERGE_BAM_ALIGNMENT.out,
            DROP_SEQ_TOOLS__CONVERT_TO_REFFLAT.out
        )
        DROP_SEQ_TOOLS__DETECT_REPAIR_BARCODE_SYNTHESIS_ERRORS( DROP_SEQ_TOOLS__TAG_READ_WITH_GENE_EXON.out )
        FINAL_BAM = DROP_SEQ_TOOLS__DETECT_REPAIR_BARCODE_SYNTHESIS_ERRORS.out.bam
        DROP_SEQ_TOOLS__BAM_TAG_HISTOGRAM( FINAL_BAM )
        DROPLET_UTILS__BARCODE_SELECTION( DROP_SEQ_TOOLS__BAM_TAG_HISTOGRAM.out )
        a = FINAL_BAM.combine(DROPLET_UTILS__BARCODE_SELECTION.out.selectedCellBarcodesByKnee, by: 0)
        b = FINAL_BAM.combine(DROPLET_UTILS__BARCODE_SELECTION.out.selectedCellBarcodesByInflection, by: 0)
        c = FINAL_BAM.combine(selectedBarcodesByCustom, by: 0)
        DROP_SEQ_TOOLS__DIGITAL_EXPRESSION(
            a.mix(b,c)
        )
}
