//====================================================================//
// process for decontamination of reads using kraken2 and kraken tools
//====================================================================//

// the objective here is to check the reads for contamination and remove any reads that are classified as human contamination
// later, it will be possible in the options to chose another type of decontamination too or at the opposite to select only the contaminated reads for example



//===========================================================================================================

process CONTAMINATION {

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_kraken_report.txt"), path("${sample}_kraken_output.txt"), emit: results

    script:
    """
    module load bioinfo/Kraken2/2.17.1 
   
    kraken2 --db ${params.kraken_db} \\
            --paired ${reads} \\
            --threads ${task.cpus} \\
            --report "${sample}_kraken_report.txt" \\
            --output "${sample}_kraken_output.txt"
    """
}


process REMOVE {

    input:
    tuple val(sample), path(reads), path(kraken_report), path(kraken_output)

    output:
    tuple val(sample), path("*.unclassified_*.f*q.gz"), emit: clean_reads

    script:

    """
    module load bioinfo/Kraken2/2.17.1
    module load bioinfo/KrakenTools/d4a2fbe
    module load devel/python/Python-3.12.4

    extract_kraken_reads.py \\
             -k ${kraken_output} \\
             -r ${kraken_report} \\
             -s1 ${reads[0]} \\
             -s2 ${reads[1]} \\
             -t 9606 \\
             --exclude \\
             --include-children \\
             -o "${sample}.unclassified_1.fq" \\
             -o2 "${sample}.unclassified_2.fq" \\
             --fastq-output

    gzip ${sample}.unclassified_1.fq
    gzip ${sample}.unclassified_2.fq
    """
}


process CHECK {

    input:
    tuple val(sample), path(clean_reads)

    output:
    path "${sample}_clean_kraken_report.txt", emit: report
    path "${sample}_clean_kraken_output.txt", emit: output

    script:
    """
    module load bioinfo/Kraken2/2.17.1
    
    kraken2 --db ${params.kraken_db} \\
            --paired ${clean_reads} \\
            --threads ${task.cpus} \\
            --report "${sample}_clean_kraken_report.txt" \\
            --output "${sample}_clean_kraken_output.txt"

    """
}
