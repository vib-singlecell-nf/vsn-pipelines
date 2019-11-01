nextflow.preview.dsl=2

process PUBLISH {

    publishDir "out", mode: 'link', overwrite: true, saveAs: { filename -> params.out_filename }

    input:
    file f

    output:
    file f

    """
    """
}

process TOUCH {

    output:
    file("foo.txt")
    
    """
    touch foo.txt
    """

}

workflow {
    main:
        PUBLISH( TOUCH() )
}