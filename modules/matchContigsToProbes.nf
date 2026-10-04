//=======================================//
// process pour match contigs to probes
//=======================================//

// the objective here is to realise an alignment with lastz tool between contigs obtained and probes


//===========================================================================================================

process ALIGNMENT_LASTZ {

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/88/88752eb83bf732d6ecfa649fc86578142c50b6dc53b4e05c986b1f307430dd6d/data'
    // oras://community.wave.seqera.io/library/lastz_bash_gzip_python:f9d07a8edb2d3ffc

    input:
    tuple val(sample), path (contigs)
    path probes

    output:
    tuple val(sample), path ("${sample}_fasta_vs_probes.sam"), emit: lastz_results

    script:
    
    """
    #module load bioinfo/LASTZ/1.04.22
    
    lastz ${probes}[multiple,unmask] \\
                "${contigs}"[unmask] \\
                --gapped \\
                --format=softsam \\
                --ambiguous=iupac \\
                --output="${sample}_fasta_vs_probes.sam"
    """

}