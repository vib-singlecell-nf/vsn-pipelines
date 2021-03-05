nextflow.enable.dsl=2

// binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/template/bin/" : ""

toolParams = params.tools.bwamaptools

process BWAMAPTOOLS__ADD_BARCODE_TAG {

    container toolParams.container
    label 'compute_resources__default','compute_resources__24hqueue'
    // todo: add storeDir instead of publishDir

    input:
        tuple val(sampleId),
              path(bam)

    output:
        tuple val(sampleId),
              path("${sampleId}.bwa.possorted.bam")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, toolParams.add_barcode_as_tag)
        processParams = sampleParams.local
        """
        samtools view -h ${bam} \
            | mawk '/^@/ {print;next} {
                N=split(\$1,n,"${processParams.delimiter_to_get_barcode_block}");
                NN=split(n[1],nbc,"${processParams.delimiter_to_split_barcodes}");
                ucorr_bc="";
                corr_bc="";
                if(NN==1){ # BioRad data with format {corrected_bc}_qname
                    corr_bc="\t${processParams.corrected_bc_tag}:Z:" nbc[1];
                }
                if(NN==2){ # standard format with {uncorrected_bc}-{corrected_bc}_qname
                    ucorr_bc="\t${processParams.uncorrected_bc_tag}:Z:" nbc[1];
                    if(length(nbc[2])>0){ # skip if corrected bc is empty: {uncorrected_bc}-_qname
                        corr_bc="\t${processParams.corrected_bc_tag}:Z:" nbc[2];
                    }
                }
                sub(/^[^${processParams.delimiter_to_get_barcode_block}]+${processParams.delimiter_to_get_barcode_block}/,"",\$0);
                print \$0 ucorr_bc corr_bc;
              }' \
            | samtools view -bS - -o ${sampleId}.bwa.possorted.bam
        """
}

