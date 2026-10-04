//====================================================================//
// process for decontamination of reads using kraken2 and kraken tools
//====================================================================//

// the objective here is to check the reads for contamination and remove any reads that are classified as human contamination
// later, it will be possible in the options to chose another type of decontamination too or at the opposite to select only the contaminated reads for example



//===========================================================================================================

process CONTAMINATION {

    maxForks 2

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e0ac1be3c45d53e5a1ed60a8cce72cd9832ce204bb44b0e908654866b712bdc9/data'
    // oras://community.wave.seqera.io/library/kraken2_bash_gzip_python:87f879546f2641d8

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_kraken_report.txt"), path("${sample}_kraken_output.txt"), emit: results
    path "${sample}_contam.tsv", emit: stats

    script:
    """
    #module load bioinfo/Kraken2/2.17.1 
   
    kraken2 --db ${params.kraken_db} \\
            --paired ${reads} \\
            --threads ${task.cpus} \\
            --report "${sample}_kraken_report.txt" \\
            --output "${sample}_kraken_output.txt"


    human_pct=\$(awk -F'\\t' '\$5 == 9606 {print \$1}' ${sample}_kraken_report.txt | tr -d ' ')
    human_pct=\${human_pct:-0}

    bacteria_pct=\$(awk -F'\\t' '\$5 == 2 {print \$1}' ${sample}_kraken_report.txt | tr -d ' ')
    bacteria_pct=\${bacteria_pct:-0}

    # potentiel problème si d'autres taxons ont un 2 dans leur ID, mais à voir

    printf "%s\\t%s\\t%s\\n" "${sample}" "\$human_pct" "\$bacteria_pct" > ${sample}_contam.tsv
    """
}


process REMOVE {

    maxForks 4

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/75/752b51b96c1e6b95264be32e9056d30f3d7f0c0e7ad756b269a8854704dc9af9/data'
    // oras://community.wave.seqera.io/library/kraken2_krakentools_bash_gzip_python:0663353035cc296d

    input:
    tuple val(sample), path(reads), path(kraken_report), path(kraken_output)

    output:
    tuple val(sample), path("*.unclassified_*.f*q.gz"), emit: clean_reads
    path "${sample}_removed_stats.tsv", emit: stats

    script:

    """
    #module load bioinfo/Kraken2/2.17.1
    #module load bioinfo/KrakenTools/d4a2fbe
    #module load devel/python/Python-3.12.4

    extract_kraken_reads.py \\
             -k ${kraken_output} \\
             -r ${kraken_report} \\
             -s1 ${reads[0]} \\
             -s2 ${reads[1]} \\
             -t 9606 \\
             --exclude \\
             --include-children \\
             -o "${sample}.unclassified_1.fq" \\
             -o2 "${sample}.unclassified_2.fq" \\
             --fastq-output > /dev/null


    # Partie statistiques
    reads_before=\$(( \$(zcat ${reads[0]} | wc -l) / 4 ))
    reads_after=\$(( \$(wc -l < ${sample}.unclassified_1.fq) / 4 ))
    reads_removed=\$(( reads_before - reads_after ))

    pct_removed=\$(awk -v a=\$reads_removed -v b=\$reads_before \\
        'BEGIN{ if (b>0) printf "%.2f", (a/b)*100; else print 0 }')

    printf "%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "${sample}" "\$reads_before" "\$reads_after" "\$reads_removed" "\$pct_removed" \\
        > ${sample}_removed_stats.tsv



    # Pour zipper les fichiers de sortie 
    gzip ${sample}.unclassified_1.fq
    gzip ${sample}.unclassified_2.fq

    rm -f "\$(readlink -f ${kraken_output})" #si on veut garder en cache la partie décontamination et la suite il ne faut pas supprimer le fichier kraken_output
    """
}


process CHECK {

    maxForks 2

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e0ac1be3c45d53e5a1ed60a8cce72cd9832ce204bb44b0e908654866b712bdc9/data'
    // oras://community.wave.seqera.io/library/kraken2_bash_gzip_python:87f879546f2641d8

    input:
    tuple val(sample), path(clean_reads)

    output:
    path "${sample}_clean_kraken_report.txt", emit: report
    path "${sample}_decontam.tsv", emit: stats

    script:
    """
    #module load bioinfo/Kraken2/2.17.1
    
    kraken2 --db ${params.kraken_db} \\
            --paired ${clean_reads} \\
            --threads ${task.cpus} \\
            --report "${sample}_clean_kraken_report.txt" \\
            --output "${sample}_clean_kraken_output.txt"

    
    rm -rf ${sample}_clean_kraken_output.txt

    human_pct=\$(awk -F'\\t' '\$5 == 9606 {print \$1}' ${sample}_clean_kraken_report.txt | tr -d ' ')
    human_pct=\${human_pct:-0}

    if awk -v h="\$human_pct" 'BEGIN {exit !(h > 0)}'; then
        echo "Contamination wasn't well removed" >&2
    fi

    bacteria_pct=\$(awk -F'\\t' '\$5 == 2 {print \$1}' ${sample}_clean_kraken_report.txt | tr -d ' ')
    bacteria_pct=\${bacteria_pct:-0}

    # potentiel problème si d'autres taxons ont un 2 dans leur ID, mais à voir

    printf "%s\\t%s\\t%s\\n" "${sample}" "\$human_pct" "\$bacteria_pct" > ${sample}_decontam.tsv
    """
}


process CONTAM_SUMMARY {

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dc/dcba98ee6037ee71e2852483e0d2fd6275db255218b3425446b5788d30ef8362/data'
    // oras://community.wave.seqera.io/library/bash_gzip_pandoc_python:12ff9a2d0fc88bf2

    input:
    path before_stats
    path removed_stats
    path after_stats

    output:
    path "contamination_summary.tsv", emit: table
    path "contamination_report.md",   emit: report
    path "contamination_report.html", emit: report_html

    script:
    """
    #module load tools/Pandoc/3.1.2

    echo -e "sample\\thuman_pct\\tbacteria_pct" > before.tsv
    cat ${before_stats} >> before.tsv

    echo -e "sample\\treads_before\\treads_after\\treads_removed\\tpct_removed" > removed.tsv
    cat ${removed_stats} >> removed.tsv

    echo -e "sample\\thuman_pct\\tbacteria_pct" > after.tsv
    cat ${after_stats} >> after.tsv



    join -t \$'\\t' \\
        <(tail -n +2 before.tsv | sort) \\
        <(tail -n +2 removed.tsv | sort) \\
        > joined.tsv

    echo -e "sample\\thuman_pct\\tbacteria_pct\\treads_before\\treads_after\\treads_removed\\tpct_removed" \\
        > contamination_summary.tsv
    cat joined.tsv >> contamination_summary.tsv

    {
    echo "Ce rapport présente les résultats de l'analyse de contamination des échantillons. Il inclut les pourcentages de contamination humaine et bactérienne avant la décontamination, ainsi que le nombre de reads avant et après le processus de décontamination. Si la décontamination a été mal réalisée cela sera aussi indiquée."
    echo ""
    echo "## Ressources utilisées"
    echo ""
    echo "Pour avoir la consommation réelle par processus il faut se référer à : \\`pipeline_info/\\` (colonnes \\`%cpu\\`, \\`rss\\`, \\`peak_rss\\` et \\`peak_vmem\\`)."
    echo ""
    echo "## Résumé de la décontamination"
    echo "| Sample | Humain % | Bactérie % | Reads avant | Reads après | Reads retirés | % retiré | Statut |"
    echo "|---|---|---|---|---|---|---|---|"
    tail -n +2 contamination_summary.tsv | while IFS=\$'\\t' read -r sample human bacteria before after removed pct; do
        status="OK"
        awk -v h="\$human" 'BEGIN{exit !(h>20)}'  && status="CRITIQUE : contamination humaine > 20%"
        awk -v h="\$human" 'BEGIN{exit !(h>10)}'  && [ "\$status" = "OK" ] && status="ATTENTION : contamination humaine > 10%"
        echo "| \$sample | \$human | \$bacteria | \$before | \$after | \$removed | \$pct | \$status |"
    done
    echo ""
    echo "## Processus de décontamination"
    if awk -F '\\t' 'NR > 1 && \$2 > 1 {found=1} END {exit !found}' after.tsv; then
        echo "ATTENTION : la décontamination n'a pas été efficace pour certains échantillons."
    else
        echo "La décontamination a été efficace pour tous les échantillons."
    fi
    } > contamination_report.md

    pandoc contamination_report.md -f markdown -t html -s \\
        --metadata title="Rapport de décontamination" \\
        -o contamination_report.html

    """
}
