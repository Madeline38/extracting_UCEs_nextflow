//====================================================//
// process for fastqc with raw data no decontaminated
//====================================================//

// the objective here is to do a multi fastqc to see if the reads have been well trimmed before



//===========================================================================================================


process FASTQC_RUN_AFTER {

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fb/fbd6555406e6243020a249373f0655364168c171d9f225b26fcb7dd2f053e212/data'
    // oras://community.wave.seqera.io/library/fastqc_bash_gzip_python:6b829ecf53e8f7d8

    input:
    tuple val(sample), path(clean_reads)

    output:
    path "fastqc_${sample}", emit: report

    script:
    """
    #module load bioinfo/FastQC/0.12.1

    mkdir -p fastqc_${sample}
    
    fastqc \\
         ${clean_reads} \\
         -o fastqc_${sample} \\
         --threads ${task.cpus}
    """
}

process MULTIQC_RUN_AFTER {

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/24/24b4858a6b8e8c7d28ab648c38ca2f3fc490d14897836e0c6b9f5600f8ae5e33/data'
    // oras://community.wave.seqera.io/library/multiqc_bash_gzip_python:789c63477a41a555

    input:
    path reports

    output:
    path "multiqc_report", emit: report

    script:
    """
    #module load bioinfo/MultiQC/1.35

    multiqc ${reports} -o multiqc_report
    """
}
