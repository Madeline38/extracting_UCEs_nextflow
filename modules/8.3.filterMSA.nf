//=========================================//
// process for filtering multiple alignment 
//=========================================//

// the objective here is to filter multiple alignment to remove loci that are not present in enough samples and to generate a summary of the filtering process



//===========================================================================================================

process FILTER_MSA {
    
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/91747ec66b65b4954100cafe626e9a452eecd1af8d802892b25912d472eb1463/data'
    // oras://community.wave.seqera.io/library/bash_gzip_python:dc01a8ce17bd9a8f

    input:
    path aligned_fastas
    val n_samples

    output:
    path "filtered_*/*_aligned.fasta", optional: true, emit: filtered
    path "locus_counts.tsv",                           emit: counts
    path "locus_membership.tsv",                       emit: membership

    script:
    def min_pct = params.min_percent.toString().toDouble()
    def min_n   = Math.max(2, Math.ceil(min_pct * n_samples / 100.0 - 1e-9) as int)

    """
    mkdir -p filtered_${min_pct}

    echo -e "locus\\tn_samples\\tstatus" > locus_counts.tsv
    echo -e "locus\\tsample\\tstatus"    > locus_membership.tsv

    for f in ${aligned_fastas}; do
        n=\$(grep -c '^>' "\$f")
        if [ "\$n" -lt ${min_n} ]; then
            status=removed
        else
            status=kept
            cp "\$f" filtered_${min_pct}/
        fi
        echo -e "\${f}\\t\${n}\\t\${status}" >> locus_counts.tsv

        grep '^>' "\$f" | awk -v l="\$f" -v s="\$status" 'BEGIN{OFS="\\t"} {sub(/^>/,""); sub(/^_R_/,""); print l, \$0, s}' >> locus_membership.tsv
    done
    """
}