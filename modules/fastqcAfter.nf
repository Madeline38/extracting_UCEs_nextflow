//====================================================//
// process for fastqc with raw data no decontaminated
//====================================================//

// the objective here is to do a multi fastqc to see if the reads have been well trimmed before



//===========================================================================================================


process FASTQC_RUN_AFTER {

    tag "${sample}"

    input:
    tuple val(sample), path(clean_reads)

    output:
    path "fastqc_${sample}", emit: report

    script:
    """
    module load bioinfo/FastQC/0.12.1

    mkdir -p fastqc_${sample}
    
    fastqc \\
         ${clean_reads} \\
         -o fastqc_${sample} \\
         --threads ${task.cpus}
    """
}

process MULTIQC_RUN_AFTER {

    input:
    path reports

    output:
    path "multiqc_report", emit: report

    script:
    """
    module load bioinfo/MultiQC/1.35

    multiqc ${reports} -o multiqc_report
    """
}
