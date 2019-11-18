nextflow.preview.dsl=2

/* general reporting function: 
takes a template ipynb and adata as input,
outputs ipynb named by the value in ${report_title}
*/
process GENERATE_REPORT {

  container params.sc.scenic.container
  clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  publishDir "${params.outdir}/notebooks", mode: 'link', overwrite: true

  input:
    file ipynb
    file loom
    val report_title

  output:
    file("${report_title}.ipynb")
  script:
    """
    papermill ${ipynb} \
        ${report_title}.ipynb \
        -p FILE $loom
    """
}

process REPORT_TO_HTML {

  container params.sc.scenic.container
  clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  publishDir "${params.outdir}/notebooks", mode: 'link', overwrite: true

  input:
    file ipynb

  output:
    file "*.html"
  script:
    """
    jupyter nbconvert ${ipynb} --to html
    """
}
