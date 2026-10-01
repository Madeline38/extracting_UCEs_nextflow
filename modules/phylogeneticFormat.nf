//====================================================================//
// process for formating the files to be used in phylogenetic analysis
//=====================================================================//

// the objective here is to format the multiple alignment to be used in phylogenetic analysis with .nexus files for example



//===========================================================================================================

process NEXUS_FORMAT {

    maxForks 100

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