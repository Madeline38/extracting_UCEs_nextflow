//====================================================================//
// process for formating the files to be used in phylogenetic analysis
//=====================================================================//

// the objective here is to format the multiple alignment to be used in phylogenetic analysis with .nexus files for example



//===========================================================================================================

process NEXUS_FORMAT {

    maxForks 100

    // container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/00/00b5db083a0000c5c8b000b1bfb58d1c305ae708b43047131a54650e2ab20eaf/data'
    // oras://community.wave.seqera.io/library/bash_gzip_python:b4b02c6d19a1749b

    input:
    path aligned_fasta

    output:
    path "*.nexus", emit: nexus

    script:
    """
    module purge
    module load devel/python/Python-3.11.1
    convert_to_nexus.py ${aligned_fasta}
    """
}



process PHYLIP_FORMAT {

    maxForks 100

    // container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/00/00b5db083a0000c5c8b000b1bfb58d1c305ae708b43047131a54650e2ab20eaf/data'
    // oras://community.wave.seqera.io/library/bash_gzip_python:b4b02c6d19a1749b

    input:
    path aligned_fasta

    output:
    path "*.phy", emit: phylip

    script:
    """
    module purge
    module load devel/python/Python-3.11.1
    convert_to_phylip.py ${aligned_fasta}
    """
}
