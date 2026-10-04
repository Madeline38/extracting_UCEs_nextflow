//=======================================================//
// process for realizing multiple alignment of sequences
//========================================================//

// the objective here is to extract the sequences that aligned from the sam files and to do a multiple alignment of the sequences that have enough samples per locususing mafft 



//===========================================================================================================

process EXTRACT_FLANKED {

    maxForks 100

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/91747ec66b65b4954100cafe626e9a452eecd1af8d802892b25912d472eb1463/data'
    // oras://community.wave.seqera.io/library/bash_gzip_python:dc01a8ce17bd9a8f

    input:
    tuple val(sample), path(clean_sam)

    output:
    tuple val(sample), path("${sample}-*.fasta"), optional: true, emit: sequences

    script:
    """
    #module load devel/python/Python-3.12.4
    extract_flanked.py ${clean_sam} ${sample} ${params.flank}    
    
    """
}



process COMBINE_BY_LOCUS {

    maxForks 100

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/91747ec66b65b4954100cafe626e9a452eecd1af8d802892b25912d472eb1463/data'
    // oras://community.wave.seqera.io/library/bash_gzip_python:dc01a8ce17bd9a8f

    input:
    path all_fasta
    val n_samples

    output:
    path "locus_*.fasta", optional: true, emit: combined
    path "locus_counts.tsv",             emit: counts
    path "locus_membership.tsv",         emit: membership

    script:
    def min_n = Math.max(2, Math.ceil(params.min_percent * n_samples / 100.0 - 1e-9) as int)

    """
    for f in ${all_fasta}; do
        locus=\$(basename "\$f" .fasta)
        locus=\${locus#*-}
        cat "\$f" >> "locus_\${locus}.fasta"
    done

    echo -e "locus\\tn_samples\\tstatus" > locus_counts.tsv
    echo -e "locus\\tsample\\tstatus"    > locus_membership.tsv

    for f in locus_*.fasta; do
        [ -e "\$f" ] || continue
        n=\$(grep -c '^>' "\$f")
        if [ "\$n" -lt ${min_n} ]; then
            status=removed
        else
            status=kept
        fi
        echo -e "\${f}\\t\${n}\\t\${status}" >> locus_counts.tsv

        # qui était dans ce locus ? (à écrire AVANT le rm)
        grep '^>' "\$f" | awk -v l="\$f" -v s="\$status" 'BEGIN{OFS="\\t"} {sub(/^>/,""); print l, \$0, s}' >> locus_membership.tsv

        if [ "\$status" = removed ]; then rm "\$f"; fi
    done
    """
}

process MAFFT {

    maxForks 100

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c2c14e2ca6c1a12fe806dd8969f52484f2904d5cd457d2872c816da0597df7af/data'
    // oras://community.wave.seqera.io/library/mafft_bash_gzip_python:3e1548b2f1cdc594

    input:
    path locus_fasta

    output:
    path "*_aligned.fasta", emit: aligned

    script:
    def locus = locus_fasta.baseName.replaceFirst(/^locus_/, '')
    """
    #module load bioinfo/MAFFT/7.505

    mafft --auto --adjustdirection ${locus_fasta} > ${locus}_aligned.fasta
    """
}



process UCES_ANALYSIS {

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dc/dcba98ee6037ee71e2852483e0d2fd6275db255218b3425446b5788d30ef8362/data'
    // oras://community.wave.seqera.io/library/bash_gzip_pandoc_python:12ff9a2d0fc88bf2

    input:
    path locus_counts
    path locus_membership
    path aligned_fasta

    output:
    path "uce_summary.html", emit: summary_html
    path "uce_summary.md",   emit: summary_md

    script:
    """
    #module load tools/Pandoc/3.1.2

    uce_report.py ${locus_counts} ${locus_membership} uce_summary.md ${aligned_fasta}

    pandoc uce_summary.md -f markdown -t html -s \\
        --metadata title="Rapport final du pipeline d'extraction d'uces" \\
        -o uce_summary.html
    """
}

