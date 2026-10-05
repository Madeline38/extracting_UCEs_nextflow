***

<center>
🚧 This pipeline is still under construction 🚧
  
Feel free to try to use anything in here, but I make no promises that it will work as it is supposed to.
<center>

***


# Extracting UCEs - Nextflow

National French Projet granted by *PEPR Dynamiques de la Biodiversité Terrestre (Dynabiod)*

## Objective

This pipeline extracts **Ultra-Conserved Elements (UCEs)** from DNA sequencing reads in order to build a phylogenetic tree.

## Pipeline steps (`modules/')

| Step | Process name | Description |
|---|---|---|
| 1 | `fastqcBefore` | Fastqc and multiqc analysis before decontamination |
| 2 | `removeContamination` | Decontamination of human reads |
| 3 | `fastqcAfter` | Fastqc and multiqc analysis after decontamination |
| 4 | `asssembly` | Assembly of reads into contigs with [Megahit](https://github.com/voutcn/MEGAHIT) or [SPAdes](https://github.com/ablab/spades) |
| 5 | `qualityAssembly` | QUAST and BUSCO for analysis of the assembly |
| 6 | `matchContigsToProbes` | Matching contigs to UCE probes  |
| 7 | `removeParalogs` | Removing loci that are supposed to be paralogs |
| 8 | `multipleAlignment` | Extraction of loci and sequences (with flanking regions). And selecting loci shared by a minimum number of taxa and realizing a multiple alignment with them by loci with [MAFFT](https://github.com/GSLBiotech/mafft) |
| 9 | `phylogeneticFormat` | Put alignment in .nexus or .phylip format for phylogenetic analysis |


Outputs of each step are published to `results/Etape_1/` through `results/Etape_9/` (see [Outputs](#outputs)).

## Requirements

- A SLURM cluster with at least the `nextflow/25.10.2` module available:
  ```bash
  module avail nextflow
  module load nextflow/25.10.2
  ```
- The `singularity` module is accessible.

## Input data

- **Reads**: a directory containing paired `.fastq.gz` files, for example:
  ```
  clean-fastq/
  └── ANIC00017_0101.1.fastq.gz
  └── ANIC00017_0101.2.fastq.gz
  ```
- **Probes**: a fasta file containing the full set of UCE probes (e.g. `probes_hymeno_assembled_wrap.fasta`).



> **Note**: reads must already be adapter- and quality-trimmed 
> (e.g. with trimmomatic or illumiprocessor) before running this pipeline.

## Installation and usage

```bash
git clone <url-of-this-repo>
cd extracting_UCEs_nextflow
nextflow run main.nf
```

To resume an interrupted run, or re-run after a downstream change without recomputing everything:

```bash
nextflow run main.nf -resume
```

## Parameters (`nextflow.config`)

| Parameter | Default value | Description |
|---|---|---|
| `reads` | `'raw_data'` | Directory containing the reads |
| `skip_decontamination` | `false or true'` | By default it realize decontamination but if you don’t want to do it you can activate --skip_decontamination true |
| `kraken_db` | ` '/bank/kraken/k2_pluspfp_20251015'` | By default it is the one that exist on genoutoul cluster but you can change it if you want one of another cluster |
| `skip_assembly` | `false or true'` | If you already have your contigs or scaffolds ready you can skip all what is before with --skip_assembly true ; by default it is false |
| `contigs` | `'contigs_dir'` | If you activate skip_decontamination you have to put the path to your contigs |
| `type_assembly` | `'spades' or ‘megahit’` | You can choose between spades and megahit for the assembling part ; by default it is megahit |
| `busco_db` | `'hymenoptera_odb12'` | Database used to analyze assembly quality with presence or absence of genes ; by default it is an hymenoptera database |
| `probes` | `'probes_hymeno_assembled_wrap.fasta'` | Probe fasta file |
| `flank` | `160` | Size of the flanking region you want to obtain ; by default it is 160pb |
| `min_percent` | `75` | Minimum percentage of taxa required per locus ; by default it is 75% |
| `trimming` | `'none'` | This is trimming possibility after MSA ; by default it isn’t realize but you can ask for --trimming mafft |
| `phylogenetic_format` | `'nexus'` or `'phylip'` | You can choose the final format you want ; by default it is nexus |
| `outdir` | `'results'` | Name of the output directory |



Edit these values directly in `nextflow.config` before launching the pipeline, or override them from the command line, e.g.:

```bash
nextflow run main.nf --reads raw_data --skip_decontamination --type_assembly ‘spades’ --min_percent 50 --phylogenetic_format ‘phylip’ --outdir results_shotgun
```

## Outputs

```
results/
├── Etape_1/   # assembled contigs + assembly config
├── Etape_2/   # contig-to-probe matching results
├── Etape_3/   # taxon set, configs, logs
├── Etape_4/   # incomplete fasta per sample
├── Etape_5/   # trimmed alignments (MAFFT)
├── Etape_6/   # cleaned alignments (Gblocks)
├── Etape_7/   # post-cleaning alignments
├── Etape_8/   # loci filtered by minimum taxa
├── Etape_9/   # results in nexus or phylip format
└── pipeline_info/   # report.html, timeline.html, trace.txt
```

## Resources (e.g for 21 samples of shotgun data, each read about 5G assembly with megahit)
- Approximate runtime: about 17h 41m 40s
- Disk space: about 500 G 

## Learning more

To get familiar with Nextflow and understand each process in this pipeline: [training.nextflow.io](https://training.nextflow.io/latest/hello_nextflow/)

## Citation

If you use this tool for your analyses, please don't forget to cite this GitHub repository.

## Author

Madeline GUENE

<div>
  <img src="https://www.cnrs.fr/sites/default/files/logo/logo.svg" alt="Logo CNRS" width="75"/>
    &nbsp;&nbsp;&nbsp;&nbsp;
  <img src="https://cat.opidor.fr/images/7/7b/Screenshot_2022-07-05_at_11-35-16_M%C3%A9dias_-_LECA_-_Laboratoire_d%E2%80%99%C3%A9cologie_alpine.png" alt="Logo LECA" width="125"/>
</div>





