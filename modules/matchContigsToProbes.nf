//=======================================//
// process pour match contigs to probes
//=======================================//

process matchContigsToProbes {

    input:
    path contigs, stageAs: 'contigs/*' //met directement les contigs.fasta dans contigs/
    path probes

    output:
    path "uce-search-results", emit: uce_results
    path "uce-search-results/probe.matches.sqlite", emit: sqlite

    script:
    
    """
    mkdir -p probes
    sed '/^>/! s/[RYSWKMBDHVryswkmbdhv]/N/g' ${probes} > probes/probes_hymeno_assembled_clean.fasta

    phyluce_assembly_match_contigs_to_probes \\
        --contigs contigs \\
        --probes probes/probes_hymeno_assembled_clean.fasta \\
        --output uce-search-results \\
        --regex "^(loci_\\d+)"
    """

}

