//======================//
// process for fastqc
//======================//

// the objective here is to do a multi fastqc to see if the reads have been well trimmed before

process fastqc {

    input:
    val samples_and_dirs  // liste de tuples [sample, dir]

    output:
    path "spades-assemblies/contigs/*.contigs.fasta", emit: contigs
    path "assembly.conf", emit: conf

    script:
    def conf_lines = samples_and_dirs.join("\n")

    """
    echo "[samples]"  > assembly.conf
    echo "${conf_lines}" >> assembly.conf

    fastqc -o fastqc_report 
    
    multiqc fastqc_report


    """
}
