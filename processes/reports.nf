nextflow.preview.dsl=2

include getBaseName from '../../utils/processes/files.nf'

/* general reporting function: 
takes a template ipynb and adata as input,
outputs ipynb named by the value in ${report_title}
*/
process SC__SCANPY__GENERATE_REPORT {

  container params.sc.scanpy.container
  publishDir "${params.outdir}/notebooks", mode: 'symlink'

  input:
    file ipynb
    file adata
    val report_title

  output:
    file("${getBaseName(adata)}.${report_title}.ipynb")
  script:
    """
    papermill ${workflow.projectDir}/src/scanpy/bin/reports/${ipynb} \
        ${getBaseName(adata)}.${report_title}.ipynb \
        -p FILE $adata
    """
}

// QC report takes two inputs, so needs it own process
process SC__SCANPY__FILTER_QC_REPORT {

  container params.sc.scanpy.container
  publishDir "${params.outdir}/notebooks", mode: 'symlink'

  input:
    file(ipynb)
    file(unfiltered)
    file(filtered)
    val report_title

  output:
    file("${getBaseName(unfiltered)}.${report_title}.ipynb")
  script:
    """
    papermill ${workflow.projectDir}/src/scanpy/bin/reports/${ipynb} \
        ${getBaseName(unfiltered)}.${report_title}.ipynb \
        -p FILE1 $unfiltered -p FILE2 $filtered
    """
}

process SC__SCANPY__REPORT_TO_HTML {

  container params.sc.scanpy.container
  publishDir "${params.outdir}/notebooks", mode: 'symlink'

  input:
    file ipynb
    val report_title

  output:
    file "*.html"
  script:
    """
    jupyter nbconvert ${ipynb} --to html
    """
}

process SC__SCANPY__MERGE_REPORTS {

  container params.sc.scanpy.container
  publishDir "${params.outdir}/notebooks", mode: 'symlink'

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

