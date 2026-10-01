//=======================================================//
// process for realizing multiple alignment of sequences
//========================================================//

// the objective here is to extract the sequences that aligned from the sam files and to do a multiple alignment of the sequences that have enough samples per locususing mafft 



//===========================================================================================================

process EXTRACT_FLANKED {

    maxForks 100

    input:
    tuple val(sample), path(clean_sam)

    output:
    tuple val(sample), path("${sample}-*.fasta"), optional: true, emit: sequences

    script:
    """
    module load devel/python/Python-3.12.4
    extract_flanked.py ${clean_sam} ${sample} ${params.flank}    
    
    """
}



process COMBINE_BY_LOCUS {

    maxForks 100

    input:
    path all_fasta
    val n_samples

    output:
    path "locus_*.fasta", optional: true, emit: combined
    path "locus_counts.tsv",             emit: counts

    script:
    def min_n = Math.max(2, Math.ceil(params.min_percent * n_samples / 100.0 - 1e-9) as int)

    """
    for f in ${all_fasta}; do
        locus=\$(basename "\$f" .fasta)
        locus=\${locus#*-}
        cat "\$f" >> "locus_\${locus}.fasta"
    done

    echo -e "locus\\tn_samples\\tstatus" > locus_counts.tsv
    for f in locus_*.fasta; do
        [ -e "\$f" ] || continue
        n=\$(grep -c '^>' "\$f")
        if [ "\$n" -lt ${min_n} ]; then
            echo -e "\${f}\\t\${n}\\tremoved" >> locus_counts.tsv
            rm "\$f"
        else
            echo -e "\${f}\\t\${n}\\tkept" >> locus_counts.tsv
        fi
    done
    """
}


process MAFFT {

    maxForks 100

    input:
    path locus_fasta

    output:
    path "*_aligned.fasta", emit: aligned

    script:
    def locus = locus_fasta.baseName.replaceFirst(/^locus_/, '')
    """
    module load bioinfo/MAFFT/7.505

    mafft --auto --adjustdirection ${locus_fasta} > ${locus}_aligned.fasta
    """
}


