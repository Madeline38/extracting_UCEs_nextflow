#!/usr/bin/env nextflow

//========//
// modules
//========//

include { FASTQC_RUN ; MULTIQC_RUN                          } from './modules/fastqcBefore.nf'
include { CONTAMINATION ; REMOVE ; CHECK ; CONTAM_SUMMARY   } from './modules/removeContamination.nf'
include { FASTQC_RUN_AFTER ; MULTIQC_RUN_AFTER              } from './modules/fastqcAfter.nf'
include { SPADES_ASSEMBLY ; MEGAHIT_ASSEMBLY                } from './modules/assembly.nf'
include { QUAST ; BUSCO_DOWNLOAD ; BUSCO                    } from './modules/qualityAssembly.nf'
include { ALIGNMENT_LASTZ                                   } from './modules/matchContigsToProbes.nf'
include { REMOVE_PARALOGS ; SUMMARY_REMOVE_PARALOGS         } from './modules/removeParalogs.nf'
include { EXTRACT_FLANKED ; COMBINE_BY_LOCUS ; MAFFT        } from './modules/multipleAlignement.nf'



//=========//
// workflow
//=========//

workflow {
    
    main:

if ( !params.probes )
    error "Missing obligatory parameter : try again with --probes <name_probes.fasta>"
if ( !(params.probes ==~ /.*\.(fasta|fa)$/) )
    error "Your probe file is not a fasta file, please use the format <name_probes_file.fasta>"
probes_file = file(params.probes, checkIfExists: true)


// canaux par défaut (vides) pour les étapes pouvant être sautées
def fastqc_before_ch        = channel.empty()
def multiqc_before_ch       = channel.empty()
def contam_report_before_ch = channel.empty()
def contam_reads_before_ch  = channel.empty()
def contam_stats_ch         = channel.empty()
def clean_reads_pub_ch      = channel.empty()
def decontam_stats_ch       = channel.empty()
def contam_report_after_ch  = channel.empty()
def contam_reads_after_ch   = channel.empty()
def contam_stats_after_ch   = channel.empty()
def contam_summary_table_ch = channel.empty()
def contam_summary_full_ch  = channel.empty()
def contam_summary_html_ch  = channel.empty()
def fastqc_after_ch         = channel.empty()
def multiqc_after_ch        = channel.empty()
def logs_ch                 = channel.empty()
def contigs_ch              = channel.empty()


if ( !params.skip_assembly ) {

    if( !params.reads )
        error "Missing obligatory parameter : try again with --reads <directory_with_fastq.gz>"

    // create a channel for inputs from fastq.gz file
    reads_ch = channel.fromFilePairs(
        "${params.reads}/*{1,2}.f*q.gz",
        checkIfExists: true //check if there is pairs to avoid running empty channel
        )
        .view()

    // do check data quality with fastqc
    FASTQC_RUN(reads_ch)
    MULTIQC_RUN(FASTQC_RUN.out.report.collect())
    fastqc_before_ch  = FASTQC_RUN.out.report
    multiqc_before_ch = MULTIQC_RUN.out.report

    reads_clean_ch = reads_ch

        if ( !params.skip_decontamination ) {
            
            // decontaminate the reads from human contamination
            CONTAMINATION(reads_ch)

            REMOVE(reads_ch.join(CONTAMINATION.out.results))
            reads_clean_ch = REMOVE.out.clean_reads

            CHECK (reads_clean_ch)

            CONTAM_SUMMARY(
                CONTAMINATION.out.stats.collect(),
                REMOVE.out.stats.collect(),
                CHECK.out.stats.collect()
            )


            // do check data quality with fastqc again
            FASTQC_RUN_AFTER(reads_clean_ch)
            MULTIQC_RUN_AFTER(FASTQC_RUN_AFTER.out.report.collect())

            contam_report_before_ch = CONTAMINATION.out.results.map { sample, report, output -> tuple(sample, report) }
            contam_reads_before_ch  = CONTAMINATION.out.results.map { sample, report, output -> tuple(sample, output) }
            contam_stats_ch         = CONTAMINATION.out.stats
            clean_reads_pub_ch      = REMOVE.out.clean_reads
            decontam_stats_ch       = REMOVE.out.stats
            contam_report_after_ch  = CHECK.out.report
            contam_reads_after_ch   = CHECK.out.output
            contam_stats_after_ch   = CHECK.out.stats
            contam_summary_table_ch = CONTAM_SUMMARY.out.table
            contam_summary_full_ch  = CONTAM_SUMMARY.out.report
            contam_summary_html_ch  = CONTAM_SUMMARY.out.report_html
            fastqc_after_ch         = FASTQC_RUN_AFTER.out.report
            multiqc_after_ch        = MULTIQC_RUN_AFTER.out.report

        }


    // do assembly

    def valid_assembly = ['spades', 'megahit']
    if ( !valid_assembly.contains(params.type_assembly) ) {
        error "Unrecognised type_assembly: '${params.type_assembly}'. Options: ${valid_assembly.join(', ')}"
    }

    if (params.type_assembly == 'spades') {
        SPADES_ASSEMBLY(reads_clean_ch)
        contigs_ch = SPADES_ASSEMBLY.out.contigs
        logs_ch = SPADES_ASSEMBLY.out.logs
        }
          
    else {
        MEGAHIT_ASSEMBLY(reads_clean_ch)
        contigs_ch = MEGAHIT_ASSEMBLY.out.contigs
        logs_ch = MEGAHIT_ASSEMBLY.out.logs
          }
  }

else  {
    if ( !params.contigs )
        error "Missing obligatory parameter : try again with --contigs <directory_with_contigs>"
    contigs_ch = channel
        .fromPath("${params.contigs}/*.fasta", checkIfExists: true)
        .map { f -> tuple(f.simpleName, f) }
}

    // do analysis of the assembly
    QUAST(contigs_ch)
    BUSCO_DOWNLOAD()
    BUSCO(contigs_ch, BUSCO_DOWNLOAD.out.db.first())


    // do alignment contigs with probes
    ALIGNMENT_LASTZ (contigs_ch, probes_file)

    // do UCEs extraction
    REMOVE_PARALOGS(ALIGNMENT_LASTZ.out.lastz_results)
    SUMMARY_REMOVE_PARALOGS(REMOVE_PARALOGS.out.stats.map { sample, path -> path }.collect())

    // do multiple alignment
    EXTRACT_FLANKED(REMOVE_PARALOGS.out.clean)

    COMBINE_BY_LOCUS(
        EXTRACT_FLANKED.out.sequences
            .map { sample, files -> files }
            .flatten()
            .collect()
    )

    MAFFT(COMBINE_BY_LOCUS.out.combined.flatten())


    publish:
    fastqc_before                 = fastqc_before_ch
    multiqc_before                = multiqc_before_ch
    report_contamination_before   = contam_report_before_ch
    reads_c_or_u_before           = contam_reads_before_ch
    info_contamination            = contam_stats_ch
    clean_reads                   = clean_reads_pub_ch
    info_decontamination          = decontam_stats_ch
    report_contamination_after    = contam_report_after_ch
    reads_c_or_u_after            = contam_reads_after_ch
    info_contamination_after      = contam_stats_after_ch
    summary_decontamination       = contam_summary_table_ch
    summary_decontamination_full  = contam_summary_full_ch
    summary_decontamination_html  = contam_summary_html_ch
    fastqc_after                  = fastqc_after_ch
    multiqc_after                 = multiqc_after_ch
    assembly_contigs              = contigs_ch
    assembly_logs                 = logs_ch
    quast_report_txt              = QUAST.out.report_txt
    quast_report_tsv              = QUAST.out.report_tsv
    quast_report_html             = QUAST.out.report_html
    quast_basic_stats             = QUAST.out.basic_stats
    quast_icarus_viewers          = QUAST.out.icarus_viewers
    busco_logs                    = BUSCO.out.logs
    busco_summary_txt             = BUSCO.out.summary_txt
    busco_full_table              = BUSCO.out.full_table_tsv
    busco_plot                    = BUSCO.out.plot
    lastz_results                 = ALIGNMENT_LASTZ.out.lastz_results
    clean_sam                     = REMOVE_PARALOGS.out.clean
    excluded_sam                  = REMOVE_PARALOGS.out.excluded
    loci_remove                   = REMOVE_PARALOGS.out.loci_remove
    loci_keep                     = REMOVE_PARALOGS.out.loci_keep
    stats_paralogs                = REMOVE_PARALOGS.out.stats
    table_tsv                     = SUMMARY_REMOVE_PARALOGS.out.table
    report_md                     = SUMMARY_REMOVE_PARALOGS.out.report_md
    report_html                   = SUMMARY_REMOVE_PARALOGS.out.report_html
    flanked_sequences             = EXTRACT_FLANKED.out.sequences
    locus_combined                = COMBINE_BY_LOCUS.out.combined.flatten()
    alignments                    = MAFFT.out.aligned
}



//=========//
// output
//=========//

output {
    fastqc_before                {path 'Etape1_analyses/fastqc'}
    multiqc_before               {path 'Etape1_analyses/multiqc'}
//    report_contamination_before  {path 'Etape2_decontamination/contamination_before'}
//    reads_c_or_u_before          {path 'Etape2_decontamination/contamination_before'}
    info_contamination           {path 'Etape2_decontamination/contamination_before'}
    clean_reads                  {path 'Etape2_decontamination/new_reads/clean_reads'}
    info_decontamination         {path 'Etape2_decontamination/new_reads/stats'}
    report_contamination_after   {path 'Etape2_decontamination/contamination_after'}
    reads_c_or_u_after           {path 'Etape2_decontamination/contamination_after'}
    info_contamination_after     {path 'Etape2_decontamination/contamination_after/stats'}
    summary_decontamination      {path 'Etape2_decontamination'}
    summary_decontamination_full {path 'Etape2_decontamination'}
    summary_decontamination_html {path 'Etape2_decontamination'}
    fastqc_after                 {path 'Etape3_analyses_post_decontamination/fastqc'}
    multiqc_after                {path 'Etape3_analyses_post_decontamination/multiqc'}
    assembly_contigs             {path 'Etape4_assemblage/contigs'}
    assembly_logs                {path 'Etape4_assemblage/logs'}
    quast_report_txt             {path 'Etape5_analyses_assemblage/QUAST'}
    quast_report_tsv             {path 'Etape5_analyses_assemblage/QUAST'}
    quast_report_html            {path 'Etape5_analyses_assemblage/QUAST'}
    quast_basic_stats            {path 'Etape5_analyses_assemblage/QUAST'}
    quast_icarus_viewers         {path 'Etape5_analyses_assemblage/QUAST'}
    busco_logs                   {path 'Etape5_analyses_assemblage/BUSCO'}
    busco_summary_txt            {path 'Etape5_analyses_assemblage/BUSCO'}
    busco_full_table             {path 'Etape5_analyses_assemblage/BUSCO'}
    busco_plot                   {path 'Etape5_analyses_assemblage/BUSCO'}
    lastz_results                {path 'Etape6_alignment_probes_contigs'}
    clean_sam                    {path 'Etape7_remove_paralogs/files/sam_clean'}
    excluded_sam                 {path 'Etape7_remove_paralogs/files/sam_removed'}
    loci_remove                  {path 'Etape7_remove_paralogs/files'}
    loci_keep                    {path 'Etape7_remove_paralogs/files'}
    stats_paralogs               {path 'Etape7_remove_paralogs/files'}
    table_tsv                    {path 'Etape7_remove_paralogs'}
    report_md                    {path 'Etape7_remove_paralogs'}
    report_html                  {path 'Etape7_remove_paralogs'}
    flanked_sequences            {path 'Etape8_alignements/flanked'}
    locus_combined               {path 'Etape8_alignements/loci'}
    alignments                   {path 'Etape8_alignements/mafft'}

}
