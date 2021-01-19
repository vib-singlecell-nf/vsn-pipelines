nextflow.enable.dsl=2

process EDIRECT__SRAID_TO_SAMPLENAME {
    
    container params.edirect.container
    label 'compute_resources__default'
    maxForks 1

    input:
        val(sraId)
    output:
        tuple val(sraId), stdout
    shell:
        """
        esearch -db sra -query ${sraId} \
           | efetch --format native \
           | sed -r 's/(.*)<TITLE>(.*)<\\/TITLE>(.*)/\\2/' \
           | grep "^[^<;]" \
           | tr -d '\\n' 
        """
}
