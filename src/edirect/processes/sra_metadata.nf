nextflow.preview.dsl=2

process EDIRECT__SRAID_TO_SAMPLENAME {
    
    maxForks 1
    container params.edirect.container

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
