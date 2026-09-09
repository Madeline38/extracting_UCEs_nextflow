#!/usr/bin/env nextflow

//========//
// modules
//========//

include { fastqc                } from './modules/fastqc.nf'
include { removeContamination   } from './modules/removeContamination.nf'
include { spadesAssembly        } from './modules/spadesAssembly.nf'
include { megahitAssembly       } from './modules/megahitAssembly.nf'
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
        .map { sample, reads -> "${sample}:${reads[0].parent}" }        
        .collect()

    // do check data quality with fastqc
    fastqc(reads_ch)

    // do decontamination with kraken
    removeContamination(fastqc.out)

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



//=========//
// output
//=========//

output {
}
