//=========================================//
// process for cleaning multiple alignment 
//=========================================//

// the objective here is to clean multiple alignment by reducing gaps and removing sequences that are too divergent from the consensus sequence



//===========================================================================================================

process CLIPKIT {
  
   input:
   path clean_aln


   output:
   path "${clean_aln}_clipkit.fasta", emit: clipkit_aligned


   script:
   """   
   #change this with a container image
   module load bioinfo/ClipKIT/2.14.0


   #by default the mode is kpic
   clipkit ${clean_aln} -m ${params.mode_clipkit} -o ${clean_aln}_clipkit.fasta
   """
}

process BMGE {


   input:
   path clean_aln


   output:
   path "${clean_aln}_bmge.fasta", emit: bmge_aligned


   script:
   """   
   module load bioinfo/bmge/1.12


   bmge -i ${clean_aln} -t NT -g 0.8 -b 3 -e 0.7 -of ${clean_aln}_bmge.fasta
   """
}


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
