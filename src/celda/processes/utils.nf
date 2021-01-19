nextflow.enable.dsl=2

import java.nio.file.Paths

moduleName = "celda"
binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/${moduleName}/bin" : Paths.get(workflow.scriptFile.getParent().getParent().toString(), "${moduleName}/bin")


process SC__CELDA__DECONTX_MERGE_OUTLIER_TABLES {

    publishDir "${params.global.outdir}/data/${moduleName}", mode: 'link'
    label 'compute_resources__default'

    input:
        file("*")

    output:
        path("all.CELDA__DECONTX.Contamination_Outlier_Table.tsv.gz")

    script:
        """
        cat * > tmp.CELDA__DECONTX.Contamination_Outlier_Table.tsv
        cat \
            <(head -n1 tmp.CELDA__DECONTX.Contamination_Outlier_Table.tsv) \
            <(cat tmp.CELDA__DECONTX.Contamination_Outlier_Table.tsv | sed "/^index/d") \
            | gzip -c > all.CELDA__DECONTX.Contamination_Outlier_Table.tsv.gz
        rm tmp.CELDA__DECONTX.Contamination_Outlier_Table.tsv
        """

}
