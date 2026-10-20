//===============================//
// process for removing paralogs
//===============================//

// the objective here is to remove all paralogs from sam file we have



//===========================================================================================================

process REMOVE_PARALOGS {

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/91747ec66b65b4954100cafe626e9a452eecd1af8d802892b25912d472eb1463/data'
    // oras://community.wave.seqera.io/library/bash_gzip_python:dc01a8ce17bd9a8f

    input:
    tuple val(sample), path(sam_files)

    output:
    tuple val(sample), path("${sample}_clean.sam"),        emit: clean
    tuple val(sample), path("${sample}_exclus.sam"),       emit: excluded
    tuple val(sample), path("${sample}_loci_to_remove"),   emit: loci_remove
    tuple val(sample), path("${sample}_loci_to_keep"),     emit: loci_keep
    tuple val(sample), path("${sample}_summary.tsv"),      emit: stats

    script:

    """
    # all associations between contigs and UCEs
    grep -v '^@' ${sam_files} | cut -f1,3 > ${sample}_all


    # count of UCEs that are associated multi times with the same contig 
    sort ${sample}_all | uniq -c | awk '\$1>1' > ${sample}_same_all

    # count of contigs that are associated with more than one UCE
    cut -f1 ${sample}_all | sort |uniq -c | awk '\$1>1' > ${sample}_multi_uces

    # count of UCEs that are associated with more than one contig
    cut -f2 ${sample}_all | sort |uniq -c | awk '\$1>1' > ${sample}_multi_contigs


    # get the list of loci to remove
    awk '{print \$2}' ${sample}_multi_contigs | sort -u > ${sample}_rm_A

    awk '{print \$3}' ${sample}_same_all | sort -u > ${sample}_rm_C

    awk '{print \$2}' ${sample}_multi_uces > ${sample}_rm_B_contigs
    awk 'NR==FNR {bad[\$1]; next} (\$1 in bad) {print \$2}' ${sample}_rm_B_contigs ${sample}_all | sort -u > ${sample}_rm_B


    # combine all the loci to remove
    cat ${sample}_rm_A ${sample}_rm_B ${sample}_rm_C | sort -u > ${sample}_loci_to_remove


    # get the list of loci that we're going to keep
    cut -f2 ${sample}_all | sort -u > ${sample}_all_loci
    comm -23 ${sample}_all_loci ${sample}_loci_to_remove > ${sample}_loci_to_keep



    # create a new sam file with only the loci to keep

    awk -v loci_file="${sample}_loci_to_keep" '
    BEGIN {
        while ((getline line < loci_file) > 0) {
            gsub(/\\r/, "", line)
            gsub(/^[ \\t]+|[ \\t]+\$/, "", line)

            if (line != "")
                good_loci[line] = 1
        }
        close(loci_file)
    }

    {
        gsub(/\\r/, "", \$3)

        if (\$3 in good_loci)
            print
    }
    ' ${sam_files} > ${sample}_clean.sam


    # create a new sam file with only the loci to remove

    awk -v loci_file="${sample}_loci_to_remove" '
    BEGIN {
        while ((getline line < loci_file) > 0) {
            gsub(/\\r/, "", line)
            gsub(/^[ \\t]+|[ \\t]+\$/, "", line)

            if (line != "")
                bad_loci[line] = 1
        }
        close(loci_file)
    }

    {
        gsub(/\\r/, "", \$3)

        if (\$3 in bad_loci)
            print
    }
    ' ${sam_files} > ${sample}_exclus.sam



    # per-sample summary

    # counts per filter
    nb_all=\$(wc -l < ${sample}_all_loci)
    nb_A=\$(wc -l < ${sample}_rm_A)
    nb_B=\$(wc -l < ${sample}_rm_B)
    nb_remove=\$(wc -l < ${sample}_loci_to_remove)
    nb_keep=\$(wc -l < ${sample}_loci_to_keep)

    echo -e "sample\\tloci_total\\tloci_removed\\tloci_kept\\trm_A\\trm_B" > ${sample}_summary.tsv
    echo -e "${sample}\\t\${nb_all}\\t\${nb_remove}\\t\${nb_keep}\\t\${nb_A}\\t\${nb_B}" >> ${sample}_summary.tsv
    """
}





process SUMMARY_REMOVE_PARALOGS {

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dc/dcba98ee6037ee71e2852483e0d2fd6275db255218b3425446b5788d30ef8362/data'
    // oras://community.wave.seqera.io/library/bash_gzip_pandoc_python:12ff9a2d0fc88bf2

    input:
    path summaries

    output:
    path "remove_paralogs_summary.tsv",  emit: table
    path "remove_paralogs_report.md",    emit: report_md
    path "remove_paralogs_report.html",  emit: report_html

    script:

    """
    #module load tools/Pandoc/3.1.2

    echo -e "sample\\tloci_total\\tloci_removed\\tloci_kept\\trm_A\\trm_B\\tpct_removed\\tstatus" > remove_paralogs_summary.tsv

    for f in ${summaries}; do
        tail -n +2 "\$f"
    done | awk -F'\\t' '
    {
        pct = (\$2 > 0) ? (\$3 / \$2 * 100) : 0

        status = "OK"
        if (pct > 70)  status = "Beaucoup de loci exclus (plus de 70%)"
        if (\$4 < 100) status = (status == "OK") ? "Peu de loci restants (moins de 100)" : status " / Peu de loci restants (moins de 100)"

        printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%.1f\\t%s\\n", \$1,\$2,\$3,\$4,\$5,\$6,pct,status
    }' >> remove_paralogs_summary.tsv

    {
    echo "Cette étape retire les UCE pour lesquels l'assemblage ne permet pas de désigner un contig unique et fiable."
    echo ""
    echo "- **rm_A** : le locus est touché par plusieurs contigs (deux contigs différents, ou un même contig aligné à deux endroits sur ce locus). On ne peut pas savoir lequel est le bon."
    echo "- **rm_B** : le contig associé à ce locus s'aligne aussi sur d'autres UCE. C'est le signe d'une zone répétée du génome ou d'un contig qui regroupe plusieurs régions distinctes."
    echo ""
    echo "Un locus peut être compté dans rm_A et rm_B en même temps, ces deux colonnes ne s'additionnent donc pas forcément au nombre total de loci exclus."
    echo ""
    echo "## Résumé par échantillon"
    echo "| Échantillon | Loci totaux | Loci exclus | Loci conservés | rm_A | rm_B | % exclu | Statut |"
    echo "|---|---|---|---|---|---|---|---|"
    tail -n +2 remove_paralogs_summary.tsv | while IFS=\$'\\t' read -r sample total removed kept a b pct status; do
        echo "| \$sample | \$total | \$removed | \$kept | \$a | \$b | \${pct}% | \$status |"
    done
    echo ""

    n_warn=\$(tail -n +2 remove_paralogs_summary.tsv | awk -F'\\t' '\$8 != "OK"' | wc -l)
    if [ "\$n_warn" -gt 0 ]; then
        echo "## A verifier : \$n_warn echantillon(s)"
        echo ""
        tail -n +2 remove_paralogs_summary.tsv | awk -F'\\t' '\$8 != "OK" {print "- **" \$1 "** : " \$8}'
    else
        echo "## Tous les échantillons sont dans les normes attendues."
    fi
    } > remove_paralogs_report.md

    pandoc remove_paralogs_report.md -f markdown -t html -s \\
        --metadata title="Rapport de suppression des paralogues" \\
        -o remove_paralogs_report.html
    """
}