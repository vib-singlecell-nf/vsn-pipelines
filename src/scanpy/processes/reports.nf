nextflow.preview.dsl=2

/* general reporting function: 
takes a template ipynb and adata as input,
outputs ipynb named by the value in ${report_title}
*/
process SC__SCANPY__GENERATE_REPORT {

  container params.sc.scanpy.container
  clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  publishDir "${params.outdir}/notebooks/intermediate", mode: 'link', overwrite: true

  input:
    file ipynb
    tuple val(id), file(adata)
    val report_title
  output:
    file("${id}.${report_title}.ipynb")
  script:
    """
    papermill ${ipynb} \
        ${id}.${report_title}.ipynb \
        -p FILE $adata
    """
}

// QC report takes two inputs, so needs it own process
process SC__SCANPY__FILTER_QC_REPORT {

  container params.sc.scanpy.container
  clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  publishDir "${params.outdir}/notebooks/intermediate", mode: 'link', overwrite: true

  input:
    file(ipynb)
    tuple val(id), file(unfiltered), file(filtered)
    val report_title
  output:
    file("${id}.${report_title}.ipynb")
  script:
    """
    papermill ${ipynb} \
        ${id}.${report_title}.ipynb \
        -p FILE1 $unfiltered -p FILE2 $filtered
    """
}

process SC__SCANPY__REPORT_TO_HTML {

  container params.sc.scanpy.container
  clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  publishDir "${params.outdir}/notebooks/intermediate", mode: 'link', overwrite: true
  // copy final "merged_report" to notbooks root:
  publishDir "${params.outdir}/notebooks", pattern: '*merged_report*', mode: 'link', overwrite: true

  input:
    file ipynb
  output:
    file "*.html"
  script:
    """
    jupyter nbconvert ${ipynb} --to html
    """
}

process SC__SCANPY__MERGE_REPORTS {

  container params.sc.scanpy.container
  clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  publishDir "${params.outdir}/notebooks/intermediate", mode: 'link', overwrite: true
  // copy final "merged_report" to notbooks root:
  publishDir "${params.outdir}/notebooks", pattern: '*merged_report*', mode: 'link', overwrite: true

  input:
    file(ipynbs)
    val report_title

  output:
    file "${report_title}.ipynb"
  script:
    """
    nbmerge ${ipynbs} -o "${report_title}.ipynb"
    """
}

