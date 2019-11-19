nextflow.preview.dsl=2

if(!params.containsKey("test")) {
    binDir = "${workflow.projectDir}/src/dropletutils/bin/"
} else {
    binDir = ""
}

process SC__DROPLET_UTILS__EMPTY_DROPS {
    
    container params.sc.dropletutils.container

    input:
        tuple val(id), file(f)
    output:
		tuple val(id), file("${id}.SC__DROPLET_UTILS__EMPTY_DROPS.sce.rds")
    script:
        """
        Rscript ${binDir}empty_drops.R \
            --rds-file-path ${f} \
            --lower ${params.sc.dropletutils.empty_drops.lower} \
            --fdr-threshold ${params.sc.dropletutils.empty_drops.fdr_threshold} \
            --output "${id}.SC__DROPLET_UTILS__EMPTY_DROPS.sce.rds"
        """    
}
