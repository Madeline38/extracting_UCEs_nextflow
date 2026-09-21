#!/usr/bin/env nextflow

//========//
// modules
//========//

include { fastqcBefore          } from './modules/fastqcBefore.nf'
include { removeContamination   } from './modules/removeContamination.nf'
include { fastqcAfter           } from './modules/fastqcAfter.nf'
include { spadesAssembly        } from './modules/spadesAssembly.nf'
include { megahitAssembly       } from './modules/megahitAssembly.nf'
include { qualityAssembly       } from './modules/qualityAssembly.nf'
include { matchContigsToProbes  } from './modules/matchContigsToProbes.nf'
include { removeParalogs        } from './modules/removeParalogs.nf'



//=========//
// workflow
//=========//

workflow {

    main:
    // create a channel for inputs from fastq.gz file
    reads_ch = channel
        .fromFilePairs("${params.input}/*-READ{1,2}.fastq.gz")
        .view()

    // do check data quality with fastqc
    fastqcBefore(reads_ch)

    // do decontamination with kraken
    removeContamination(fastqc.out)

    // do re check data quality with fastqc
    fastqcAfter(removeContamination.out)

    // do assembly
    if (params.skip_assembly == 'false') {
        if (params.type_assembly == 'spades') {
          spadesAssembly(clean.out)
          contigs_ch = spadesAssembly.out.contigs
          }
          
        else {
          megahitAssembly(clean.out)
          contigs_ch = megahitAssembly.out.contigs
          }
  }

    else  {
         contigs_ch = channel
           .fromPath("${params.contigs}/*.contigs.fasta")
           .collect()
        assembly_conf_ch = channel.empty()
    }


    // do match contigs to probes
    matchContigsToProbes(contigs_ch, 
                         file(params.probes))

    // remove all the paralogs
    removeParalogs (matchContigsToProbes.out)



    publish:
    fastqc_before                 = fastqcBefore.out.report
    report_decontamination_before = removeDecontamination.out.before
    data_decontaminated           = removeDecontamination.out.data
    report_decontamination_after  = removeDecontamination.out.after
    fastqc_after                  = fastqcAfter.out.report
    contigs                       = assembly_ch
    quast                         = qualityAssembly.out.quast
    busco                         = qualityAssembly.out.busco
    alignment                     = matchContigsToProbes.out.lastz
    extraction_uces               = removeParalogs.out.data
    removed_paralogs              = removeParalogs.out.removed    
    
    



//=========//
// output
//=========//

output {
}
