//=========================//
// process for the assembly
//=========================//

// the objective here is to do an assembly of the reads using spades or megahit



//===========================================================================================================

process SPADES_ASSEMBLY {

    maxForks 2

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fed69f1d0cebfa54c4d0500b3415844dce34f36a7fe3e35223bc38b3d38496f5/data'
    // oras://community.wave.seqera.io/library/spades_bash_gzip_python:08863429beaea175

    input:
    tuple val(sample), path(clean_reads)

    output:
    tuple val(sample), path("spades/${sample}/scaffolds.fasta"), emit: contigs
    tuple val(sample), path("spades/${sample}/spades.log"),      emit: logs

    script:

    """
    #module load devel/python/Python-3.12.4
    #module load bioinfo/SPAdes/4.2.0

    CORRECTION_TIMEOUT=\$((3*3600))
    CHECK_INTERVAL=300

    echo "===================================================================="
    echo "Échantillon      : ${sample}"
    echo "Noeud            : \$(hostname)"
    echo "Date de début    : \$(date --iso-8601=seconds)"
    echo "===================================================================="

    kill_tree() {
        local pid=\$1
        local pgid
        pgid=\$(ps -o pgid= -p "\$pid" 2>/dev/null | tr -d ' ')
        if [ -n "\$pgid" ]; then
            echo "Envoi de SIGTERM au groupe de process \$pgid"
            kill -TERM -"\$pgid" 2>/dev/null
            sleep 10
            kill -KILL -"\$pgid" 2>/dev/null
        fi
    }

    watch_correction_timeout() {
        local pid=\$1
        local elapsed=0

        while kill -0 "\$pid" 2>/dev/null; do
            sleep \$CHECK_INTERVAL
            elapsed=\$((elapsed + CHECK_INTERVAL))

            if [ -f "spades/${sample}/spades.log" ] && grep -q "===== Assembling started" "spades/${sample}/spades.log"; then
                wait "\$pid"
                SPADES_EXIT=\$?
                TIMED_OUT=0
                return
            fi

            if [ \$elapsed -ge \$CORRECTION_TIMEOUT ]; then
                echo "Correction d'erreurs > 3h pour ${sample}, on bascule en --only-assembler"
                kill_tree "\$pid"
                sleep 5
                TIMED_OUT=1
                SPADES_EXIT=124
                return
            fi
        done

        wait "\$pid"
        SPADES_EXIT=\$?
        TIMED_OUT=0
    }

    echo "======== ${sample} lancé ! =========="
    setsid spades.py -t ${task.cpus} -m ${task.memory.toGiga()} --checkpoints last \\
        -k 21,33,55 \\
        -1 ${clean_reads[0]} -2 ${clean_reads[1]} \\
        -o spades/${sample} &
    SPADES_PID=\$!
    watch_correction_timeout \$SPADES_PID

    if [ \$TIMED_OUT -eq 1 ]; then
        echo "======== Relance de ${sample} en --only-assembler après timeout ========"
        mv "spades/${sample}" "spades/${sample}_correction_timeout_\$(date +%s)"
        spades.py --only-assembler -t ${task.cpus} -m ${task.memory.toGiga()} --checkpoints last \\
            -k 21,33,55 \\
            -1 ${clean_reads[0]} -2 ${clean_reads[1]} \\
            -o spades/${sample}
        STATUS=\$?
    else
        STATUS=\$SPADES_EXIT
    fi

    echo "Date de fin       : \$(date --iso-8601=seconds)"

    if [ \$STATUS -eq 0 ]; then
        rm -f "\$(readlink -f ${clean_reads[0]})"
        rm -f "\$(readlink -f ${clean_reads[1]})"
        rm -rf spades/${sample}/{corrected,tmp,misc,K21,K33,K55}
        echo "======== ${sample} : terminé avec succès ========"
    else
        echo "======== ${sample} : ÉCHEC (code retour \$STATUS) ========"
        exit \$STATUS
    fi
    """
}


process MEGAHIT_ASSEMBLY {

    maxForks 3

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/81a9d5cd7c71be35af0e0bf93ff9a5362de00d1ecbd31af148bd8f9b43a46efa/data'

    input:
    tuple val(sample), path(clean_reads)

    output:
    tuple val(sample), path("megahit/${sample}/final.contigs.fa"), emit: contigs
    tuple val(sample), path("megahit/${sample}/log"),              emit: logs
    tuple val(sample), path("megahit/${sample}/done"),             emit: done

    script:

    """
    #module load devel/python/Python-3.12.4
    #module load bioinfo/MEGAHIT/1.2.9

    echo "===================================================================="
    echo "Échantillon : ${sample}"
    echo "Noeud       : \$(hostname)"
    echo "Date début  : \$(date --iso-8601=seconds)"
    echo "===================================================================="


mkdir -p megahit

echo "======== ${sample} lancé ! =========="

megahit \\
    -t ${task.cpus} \\
    -m ${task.memory.toBytes()} \\
    -1 "${clean_reads[0]}" \\
    -2 "${clean_reads[1]}" \\
    -o "megahit/${sample}"

STATUS=\$?

if [ \$STATUS -eq 0 ] && [ -f "megahit/${sample}/done" ]; then
    echo "======== Nettoyage des contigs intermédiaires =========="
    rm -rf "megahit/${sample}/intermediate_contigs"
    rm -f "\$(readlink -f ${clean_reads[0]})"
    rm -f "\$(readlink -f ${clean_reads[1]})"
    echo "======== ${sample} fini ! =========="
else
    echo "======== ${sample} échec... =========="
    echo "Code retour MEGAHIT : \${STATUS}"
    exit \$STATUS
fi

echo "Date de fin : \$(date --iso-8601=seconds)"
echo "======== ${sample} : terminé avec succès ========"

"""
}