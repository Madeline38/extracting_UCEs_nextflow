//=======================================//
// process pour match contigs to probes
//=======================================//

// the objective here is to realise an alignment with lastz tool between contigs obtained and probes


//===========================================================================================================

process ALIGNMENT_LASTZ {

    input:
    tuple val(sample), path (contigs)
    path probes

    output:
    tuple val(sample), path ("${sample}_fasta_vs_probes.sam"), emit: lastz_results

    script:
    
    """
    module load bioinfo/LASTZ/1.04.22
    
    lastz ${probes}[multiple,unmask] \\
                "${contigs}"[unmask] \\
                --gapped \\
                --format=softsam \\
                --ambiguous=iupac \\
                --output="${sample}_fasta_vs_probes.sam"
    """

}