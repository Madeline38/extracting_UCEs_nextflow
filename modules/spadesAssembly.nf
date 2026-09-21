//======================//
// process pour assembly
//======================//

process spadesAssembly {

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

    phyluce_assembly_assemblo_spades \\
        --conf assembly.conf \\
        --output spades-assemblies \\
        --memory 500 \\
        --cores ${task.cpus}
    """
}
