//=========================================//
// process for cleaning multiple alignment 
//=========================================//

// the objective here is to clean multiple alignment by reducing gaps and removing sequences that are too divergent from the consensus sequence



//===========================================================================================================

process TRIMAL {

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b2ca0a9afe7cb7d5b16c51bcbad327c9505e83f8e609cf5cd8fed4c927847697/data'

    // oras://community.wave.seqera.io/library/trimal_bash_gzip_python:7dac6d07eb46eebd

    input:
    path clean_aln

    output:
    path "${clean_aln}_cleaned.fasta", emit: trimmed_aligned

    script:
    """    
    #module load bioinfo/trimAl/1.4.1

    trimal -in ${clean_aln} -out ${clean_aln}_cleaned.fasta
    """
}
