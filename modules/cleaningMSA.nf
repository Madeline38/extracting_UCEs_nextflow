//=========================================//
// process for cleaning multiple alignment 
//=========================================//

// the objective here is to clean multiple alignment by reducing gaps and removing sequences that are too divergent from the consensus sequence



//===========================================================================================================

process TRIMAL {

    input:
    path clean_aln

    output:
    path "${clean_aln}_cleaned.fasta", emit: trimmed_aligned

    script:
    """    
    module load bioinfo/trimAl/1.4.1

    trimal -in ${clean_aln} -out ${clean_aln}_cleaned.fasta
    """
}
