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

    output:
    path "locus_*.fasta", optional: true, emit: combined

    script:
    """
    for f in ${all_fasta}; do
        locus=\$(basename "\$f" .fasta)
        locus=\${locus#*-}
        cat "\$f" >> "locus_\${locus}.fasta"
    done

    # on ne garde que les locus présent chez au moins deux taxons pour l'alignement multiple
    for f in locus_*.fasta; do
        [ -e "\$f" ] || continue
        n=\$(grep -c '^>' "\$f")
        if [ "\$n" -lt 2 ]; then rm "\$f"; fi
    done

    """
}

process MSA_MAFFT {

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

    mafft --auto ${locus_fasta} > ${locus}_aligned.fasta
    """
}

process SUMMARY_MSA {

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


process STATS_MSA {

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fec0bd55a6ecb1423a40ff8e151be5c500ec561600f9ea22dce41d87d41d1451/data'
    // oras://community.wave.seqera.io/library/amas_bash_python:d1ff2911c095b420

    input:
    path aligned_fasta

    output:
    path "summary_stats_MSA.txt", emit: summary_stats

    script:
    """
    AMAS.py summary \\
        -f fasta \\
        -d dna \\
        -i ${aligned_fasta} \\
        -o summary_stats_MSA.txt \\
        -c ${task.cpus}
    """
}

process PLOT_MSA_STATS {

    input:
    path summary_stats
    val n_samples

    output:
    path "msa_stats_*.png", optional: true, emit: plots

    script:
    """
    module load statistics/R/4.3.1

    plot_msa_stats.R ${summary_stats} ${n_samples}
    """
}
