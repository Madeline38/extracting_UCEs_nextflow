//=====================================================================================//
// Phylogenetic formattage : process to create a supermatrix for phylogenetic analysis
//=====================================================================================//

// warning : this part wasn't tested yet to see if trees can be generated well from the supermatrix

// tool : AMAS (v.1.0), python (v.3.15.2), bash (v.5.3)
// future : add nexus format, add partition file if needed

//===========================================================================================================

process CONCATENATE_MSA {

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fec0bd55a6ecb1423a40ff8e151be5c500ec561600f9ea22dce41d87d41d1451/data'
    //oras://community.wave.seqera.io/library/amas_bash_python:d1ff2911c095b420  

    input:
    path aligned_fasta

    output:
    path "*.${params.phylogenetic_format}", emit: concatenated

    script:
    """
    # on ne garde que le nom de l'échantillon dans les en-têtes
    for f in ${aligned_fasta}; do
        sed -i '/^>/ s/ .*//' "\$f"
    done

    AMAS.py concat \\
        -i ${aligned_fasta} \\
        -f fasta \\
        -d dna \\
        -u ${params.phylogenetic_format} \\
        -t supermatrix.${params.phylogenetic_format} \\
        -c ${task.cpus}
    """
}
