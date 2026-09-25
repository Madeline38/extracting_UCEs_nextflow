//=========================================================================//
// process for analyzing the quality of the assembly with quast and busco
//=========================================================================//

// the objective here is to understand better what we obtained with the assembly and to see if there is any warning to raise


//===========================================================================================================

process QUAST {

    input:
    tuple val(sample), path(contigs)

    output:
    tuple val(sample), path("quast_${sample}_${params.type_assembly}/report.txt"),     emit: report_txt
    tuple val(sample), path("quast_${sample}_${params.type_assembly}/report.tsv"),     emit: report_tsv
    tuple val(sample), path("quast_${sample}_${params.type_assembly}/report.html"),    emit: report_html
    tuple val(sample), path("quast_${sample}_${params.type_assembly}/basic_stats"),    emit: basic_stats
    tuple val(sample), path("quast_${sample}_${params.type_assembly}/icarus_viewers"), emit: icarus_viewers

    script:

    """
    module load bioinfo/QUAST/5.2.0

    quast.py ${contigs} --k-mer-stats -o quast_${sample}_${params.type_assembly} -t ${task.cpus}
    """
}


process BUSCO_DOWNLOAD {
    output:
    path "busco_downloads", emit: db

    script:
    """
    module load devel/Miniforge/Miniforge3
    module load bioinfo/BUSCO/6.0.0
    busco --download ${params.busco_db}
    """
}

process BUSCO {
    input:
    tuple val(sample), path(contigs)
    path busco_downloads

    output:
    tuple val(sample), path("busco_${sample}_${params.type_assembly}/logs"),                                                                                      emit: logs
    tuple val(sample), path("busco_${sample}_${params.type_assembly}/short_summary.specific.${params.busco_db}.busco_${sample}_${params.type_assembly}.txt"),     emit: summary_txt
    tuple val(sample), path("busco_${sample}_${params.type_assembly}/run_${params.busco_db}/full_table.tsv"),                                                     emit: full_table_tsv
    tuple val(sample), path("busco_${sample}_${params.type_assembly}/busco_figure.png"),                                                                          emit: plot
    
    script:
    """
    module load devel/Miniforge/Miniforge3
    module load bioinfo/BUSCO/6.0.0

    busco -i ${contigs} -o busco_${sample}_${params.type_assembly} -m genome -l ${params.busco_db} -c ${task.cpus} --offline --download_path busco_downloads
    busco --plot busco_${sample}_${params.type_assembly}/ --plot_percentages
    """
}