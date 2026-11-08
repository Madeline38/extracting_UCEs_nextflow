#!/usr/bin/env nextflow

//========//
// modules
//========//

include { FASTQC ; MULTIQC                                 		                                    } from './modules/1.qualityControl.nf'
include { CHECK ; REMOVE ; SUMMARY_DECONTAMINATION          	        	                        } from './modules/2.removeContamination.nf'
include { FASTQC_AFTER ; MULTIQC_AFTER                     		                                    } from './modules/3.qualityControlAfter.nf'
include { ASSEMBLY_SPADES ; ASSEMBLY_MEGAHIT                        		                        } from './modules/4.assembly.nf'
include { QUAST ; BUSCO_DOWNLOAD ; BUSCO                          	                                } from './modules/5.qualityAssembly.nf'
include { ALIGNMENT_LASTZ                                           		                        } from './modules/6.matchContigsToProbes.nf'
include { REMOVE_PARALOGS ; SUMMARY_REMOVE_PARALOGS                 		                        } from './modules/7.removeParalogs.nf'
include { EXTRACT_FLANKED ; COMBINE_BY_LOCUS ; MSA_MAFFT ; SUMMARY_MSA ; STATS_MSA ; PLOT_MSA_STATS } from './modules/8.1.multipleSequenceAlignment.nf'
include { CLIPKIT ; BMGE ; TRIMAL                                                                	} from './modules/8.2.trimmingMSA.nf'
include { FILTER_MSA                                                                             	} from './modules/8.3.filterMSA.nf'
include { CONCATENATE_MSA                             		                                        } from './modules/9.phylogeneticFormat.nf'



//=========//
// workflow
//=========//

workflow {
    
    main:

// check if parameters are correct
if ( !params.probes )
    error "Missing obligatory parameter : try again with --probes <name_probes.fasta>"

if ( !(params.probes ==~ /.*\.(fasta|fa)$/) )
    error "Your probe file is not a fasta file, please use the format <name_probes_file.fasta>"
probes_file = file(params.probes, checkIfExists: true)


// convert min_percent to float and check if it is between 0 and 100
def min_percent = params.min_percent.toString().toFloat()

if ( min_percent < 0 || min_percent > 100 )
    error "--min_percent doit être compris entre 0 et 100"


// empty channel by default for skipping processes
def fastqc_before_ch        = channel.empty()
def multiqc_before_ch       = channel.empty()
def contam_report_before_ch = channel.empty()
def contam_reads_before_ch  = channel.empty()
def contam_stats_ch         = channel.empty()
def clean_reads_pub_ch      = channel.empty()
def decontam_stats_ch       = channel.empty()
def contam_summary_table_ch = channel.empty()
def contam_summary_full_ch  = channel.empty()
def contam_summary_html_ch  = channel.empty()
def fastqc_after_ch         = channel.empty()
def multiqc_after_ch        = channel.empty()
def logs_ch                 = channel.empty()
def contigs_ch              = channel.empty()
def mafft_ch                = channel.empty()
def mafft_filtered_ch       = channel.empty()
def mafft_trimmed_ch        = channel.empty()
def concatenated_ch         = channel.empty()


// skip_assembly is false (by default)
if ( !params.skip_assembly ) {

    // check if reads parameter is provided
    if( !params.reads )
        error "Missing obligatory parameter : try again with --reads <directory_with_fastq.gz>"

    // create a channel for inputs from fastq.gz files (reads)
    reads_ch = channel.fromFilePairs(
        "${params.reads}/*{1,2}.f*q.gz", // reads files have to be named something like sample_1.fastq.gz and sample_2.fastq.gz / don't use "-" in the name
        checkIfExists: true              // check if there is pairs to avoid running empty channel
        )
        .view()

    // 1. check data quality with fastqc and multiqc and collect the reports
    FASTQC(reads_ch)
    MULTIQC(FASTQC.out.report.collect())
    fastqc_before_ch  = FASTQC.out.report
    multiqc_before_ch = MULTIQC.out.report

    reads_clean_ch = reads_ch   // if skip_decontamination is true, reads_clean_ch will be the same as reads_ch

        // 2. decontamination step if skip_decontamination is false (by default)
        if ( !params.skip_decontamination ) {
            
            // check where are the contamination 
            CHECK(reads_ch)

            // decontaminate the reads from human contamination
            REMOVE(reads_ch.join(CHECK.out.results))
            reads_clean_ch = REMOVE.out.clean_reads

            // summary of the decontamination process
            SUMMARY_DECONTAMINATION(
                CHECK.out.stats.collect(),
                REMOVE.out.stats.collect(),
            )


            // 3. check data quality with fastqc and multiqc with the cleaned reads and collect the reports (only if decontamination step is true)
            FASTQC_AFTER(reads_clean_ch)
            MULTIQC_AFTER(FASTQC_AFTER.out.report.collect())

            // publish the channels for output
            contam_report_before_ch = CHECK.out.results.map { sample, report, output -> tuple(sample, report) }
            contam_reads_before_ch  = CHECK.out.results.map { sample, report, output -> tuple(sample, output) }
            contam_stats_ch         = CHECK.out.stats
            clean_reads_pub_ch      = REMOVE.out.clean_reads
            decontam_stats_ch       = REMOVE.out.stats
            contam_summary_table_ch = SUMMARY_DECONTAMINATION.out.table
            contam_summary_full_ch  = SUMMARY_DECONTAMINATION.out.report
            contam_summary_html_ch  = SUMMARY_DECONTAMINATION.out.report_html
            fastqc_after_ch         = FASTQC_AFTER.out.report
            multiqc_after_ch        = MULTIQC_AFTER.out.report

        }


    // 4. assembly step

    // check if the parameter type_assembly is valid
    def valid_assembly = ['spades', 'megahit'] 
    if ( !valid_assembly.contains(params.type_assembly) ) {
        error "Unrecognised type_assembly: '${params.type_assembly}'. Options: ${valid_assembly.join(', ')}"
    }

    // spades
    if (params.type_assembly == 'spades') {
        ASSEMBLY_SPADES(reads_clean_ch)
        contigs_ch = ASSEMBLY_SPADES.out.contigs
        logs_ch = ASSEMBLY_SPADES.out.logs
        }

    // megahit
    else {
        ASSEMBLY_MEGAHIT(reads_clean_ch)
        contigs_ch = ASSEMBLY_MEGAHIT.out.contigs
        logs_ch = ASSEMBLY_MEGAHIT.out.logs
          }
  }


// if skip_assembly is true
else  {

    // check if the contigs parameter is provided
    if ( !params.contigs )
        error "Missing obligatory parameter : try again with --contigs <directory_with_contigs>"

    // create a channel for inputs from contigs files
    contigs_ch = channel
        .fromPath("${params.contigs}/*.fasta", // contigs files have to be named something like sample.fasta
        checkIfExists: true 
        ) 
        .map { f -> tuple(f.simpleName, f) } // create a tuple with the sample name and the path to the contigs file
}

// 5. quality assembly step
QUAST(contigs_ch) // analyse metrics of the assembly (N50, L50...)
BUSCO_DOWNLOAD() // download the busco database (genes of hymenoptera_odb12 database by default)
BUSCO(contigs_ch, BUSCO_DOWNLOAD.out.db.first()) // analyse the completeness of the assembly with busco


// 6. do alignment contigs with probes (reference loci)
ALIGNMENT_LASTZ (contigs_ch, probes_file)

// 7. remove paralogs
REMOVE_PARALOGS(ALIGNMENT_LASTZ.out.lastz_results)
n_samples_ch = contigs_ch.count()  // count the number of samples to use it in the FILTER_MSA process

SUMMARY_REMOVE_PARALOGS(REMOVE_PARALOGS.out.stats.map { sample, path -> path }.collect())

// 8. do multiple sequence alignment with the cleaned sam files and the loci to keep
EXTRACT_FLANKED(REMOVE_PARALOGS.out.clean)


COMBINE_BY_LOCUS(
    EXTRACT_FLANKED.out.sequences
        .map { sample, files -> files }
        .flatten()
        .collect(),
)

MSA_MAFFT(COMBINE_BY_LOCUS.out.combined.flatten())
mafft_ch = MSA_MAFFT.out.aligned

// differents possibility of trimming the MSA(clipkit, trimal, bmge) or no trimming at all
msa_ch = mafft_ch   //canal de travail
    
if ( params.trimming == 'clipkit' ) {
    CLIPKIT(mafft_ch)
    mafft_trimmed_ch = CLIPKIT.out.clipkit_aligned
}

else if ( params.trimming == 'bmge' ) {
    BMGE(mafft_ch)
    mafft_trimmed_ch = BMGE.out.bmge_aligned
}

else if ( params.trimming == 'trimal' ) {
    TRIMAL(mafft_ch)
    mafft_trimmed_ch = TRIMAL.out.trimmed_aligned
}

// si un trimming a eu lieu, on filtre les alignements trimmés
if ( params.trimming != 'none' ) {
    msa_ch = mafft_trimmed_ch
}

STATS_MSA(msa_ch.collect())
PLOT_MSA_STATS(STATS_MSA.out.summary_stats, n_samples_ch)

// par défault le filtre est à 0% donc on garde tous les locus, mais on peut changer le paramètre min_percent pour filtrer les locus avec moins de n% d'échantillons   
FILTER_MSA(msa_ch.collect(), n_samples_ch)
msa_ch = FILTER_MSA.out.filtered


SUMMARY_MSA(
    FILTER_MSA.out.counts,
    FILTER_MSA.out.membership,
    FILTER_MSA.out.filtered.collect()
    )


    // do reformat

    def valid_phylogenetic_format = ['nexus', 'phylip']
    if ( !valid_phylogenetic_format.contains(params.phylogenetic_format) ) {
        error "Unrecognised phylogenetic_format: '${params.phylogenetic_format}'. Options: ${valid_phylogenetic_format.join(', ')}"
    }

CONCATENATE_MSA(FILTER_MSA.out.filtered.collect())
concatenated_ch = CONCATENATE_MSA.out.concatenated
    


    publish:
    fastqc_before                 = fastqc_before_ch
    multiqc_before                = multiqc_before_ch
    report_contamination_before   = contam_report_before_ch
    info_contamination            = contam_stats_ch
    info_decontamination          = decontam_stats_ch
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
    locus_combined                = COMBINE_BY_LOCUS.out.combined
    alignments_all                = mafft_ch          
    alignments_trimmed            = mafft_trimmed_ch
    summary_stats_msa             = STATS_MSA.out.summary_stats
    plot_stats_msa                = PLOT_MSA_STATS.out.plots
    alignments_filtered           = mafft_filtered_ch
    uce_summary_html              = SUMMARY_MSA.out.summary_html
    uce_summary_md                = SUMMARY_MSA.out.summary_md
    locus_counts                  = FILTER_MSA.out.counts        
    locus_membership              = FILTER_MSA.out.membership  
    final_file                   = concatenated_ch
}



//=========//
// output
//=========//

output {
    fastqc_before                {path 'Etape_1_analyses/fastqc'}
    multiqc_before               {path 'Etape_1_analyses/multiqc'}
    report_contamination_before  {path 'Etape_2_decontamination/contamination'}
    info_contamination           {path 'Etape_2_decontamination/contamination/stats'}
    info_decontamination         {path 'Etape_2_decontamination/new_reads/stats'}
    summary_decontamination      {path 'Etape_2_decontamination'}
    summary_decontamination_full {path 'Etape_2_decontamination'}
    summary_decontamination_html {path 'Etape_2_decontamination'}
    fastqc_after                 {path 'Etape_3_analyses_post_decontamination/fastqc'}
    multiqc_after                {path 'Etape_3_analyses_post_decontamination/multiqc'}
    assembly_contigs             {path 'Etape_4_assemblage/contigs'}
    assembly_logs                {path 'Etape_4_assemblage/logs'}
    quast_report_txt             {path 'Etape_5_analyses_assemblage/QUAST'}
    quast_report_tsv             {path 'Etape_5_analyses_assemblage/QUAST'}
    quast_report_html            {path 'Etape_5_analyses_assemblage/QUAST'}
    quast_basic_stats            {path 'Etape_5_analyses_assemblage/QUAST'}
    quast_icarus_viewers         {path 'Etape_5_analyses_assemblage/QUAST'}
    busco_summary_txt            {path 'Etape_5_analyses_assemblage/BUSCO'}
    busco_full_table             {path 'Etape_5_analyses_assemblage/BUSCO'}
    busco_plot                   {path 'Etape_5_analyses_assemblage/BUSCO'}
    lastz_results                {path 'Etape_6_alignment_probes_contigs'}
    clean_sam                    {path 'Etape_7_remove_paralogs/files/sam_clean'}
    excluded_sam                 {path 'Etape_7_remove_paralogs/files/sam_removed'}
    loci_remove                  {path 'Etape_7_remove_paralogs/files'}
    loci_keep                    {path 'Etape_7_remove_paralogs/files'}
    stats_paralogs               {path 'Etape_7_remove_paralogs/files'}
    table_tsv                    {path 'Etape_7_remove_paralogs'}
    report_md                    {path 'Etape_7_remove_paralogs'}
    report_html                  {path 'Etape_7_remove_paralogs'}
    flanked_sequences            {path 'Etape_8_alignments/flanked'}
    locus_combined               {path 'Etape_8_alignments/loci'}
    alignments_all               {path 'Etape_8_alignments/mafft/all'}
    alignments_trimmed           {path 'Etape_8_alignments/mafft/trimmed'}
    summary_stats_msa            {path 'Etape_8_alignments/stats'}
    plot_stats_msa               {path 'Etape_8_alignments/stats'}
    alignments_filtered          {path 'Etape_8_alignments/mafft/filtered'}
    locus_membership             {path 'Etape_8_alignments'}
    uce_summary_html             {path 'Etape_9_final_files'}
    uce_summary_md               {path 'Etape_9_final_files'}
    locus_counts                 {path 'Etape_8_alignments'}
    final_file                   {path 'Etape_9_final_files/'}
}
