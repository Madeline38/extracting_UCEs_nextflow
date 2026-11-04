//====================================================================//
// process for formating the files to be used in phylogenetic analysis
//=====================================================================//

// the objective here is to format the multiple alignment to be used in phylogenetic analysis with .nexus files for example



//===========================================================================================================


process CONCATENATE_MSA {

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fec0bd55a6ecb1423a40ff8e151be5c500ec561600f9ea22dce41d87d41d1451/data'
    //oras://community.wave.seqera.io/library/amas_bash_python:d1ff2911c095b420   amas v1.0  python 3.15.2   bash 5.3

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
