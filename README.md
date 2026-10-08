***

<center>
🚧 This pipeline is still under construction 🚧
  
Feel free to try to use anything in here, but I make no promises that it will work as it is supposed to.

Without any modifications it will only work on genotoul cluster
</center>

***


# Extracting UCEs - Nextflow

French national project funded by *PEPR Dynamiques de la Biodiversité Terrestre (Dynabiod)*

## Objective

{name_pipeline} is a bioinformatics pipeline that can be used to extract **Ultra-Conserved Elements (UCEs)** from DNA sequencing reads in order to build a phylogenetic tree.

## Pipeline steps (`modules/')

| Step | Process name | Tools | Description |
|---|---|---|---|
| 1 | `fastqcBefore` | [Fastqc](https://github.com/s-andrews/fastqc) & [multiqc](https://github.com/multiqc/multiqc) | Assess read quality before decontamination (read-input mode only). |
| 2 | `removeContamination` | Kraken2 | Identify and remove human-classified reads; skipped with `--skip_decontamination true`. |
| 3 | `fastqcAfter` | [Fastqc](https://github.com/s-andrews/fastqc) & [multiqc](https://github.com/multiqc/multiqc) | Assess read quality after decontamination (when decontamination is enabled). |
| 4 | `assembly` | [MEGAHIT](https://github.com/voutcn/MEGAHIT) or [SPAdes](https://github.com/ablab/spades) | Assemble paired reads into contigs, unless existing contigs are supplied with `--skip_assembly true`. |
| 5 | `qualityAssembly` | [QUAST](https://github.com/ablab/quast) & [BUSCO](https://busco.ezlab.org/) | Evaluate assembly statistics and completeness of the contigs. |
| 6 | `matchContigsToProbes` | LASTZ | Align contigs against the supplied UCE probe FASTA. |
| 7 | `removeParalogs` | Paralogue filtering script | Exclude loci with ambiguous or multiple contig-to-locus matches. |
| 8 | `multipleAlignment` | Flank extraction, locus filtering, [MAFFT](https://github.com/GSLBiotech/mafft) | Extract flanked sequences, retain loci meeting `min_percent`, align them by locus, and produce the UCE summary report. TrimAl can optionally trim alignments. |
| 9 | `phylogeneticFormat` | Biopython package | Convert alignments to Nexus or Phylip. |

Outputs are published under the `outdir` directory (default: `results/`; see [Outputs](#outputs)).


## Requirements

- A SLURM cluster configured for the `workq` queue, as specified in `nextflow.config` which normally is related to Genotoul cluster
- Nextflow 25.10.2 (or a later version) available as a module, for example:
  ```bash
  module avail nextflow
  module load bioinfo/Nextflow/26.04.6
  ```
- Singularity/Apptainer available for running the containers :
  ```bash
  module avail apptainer
  module load containers/Apptainer/1.4.1
  ```
- Java available :
  ```bash
  module avail java
  module load devel/java/24.0.2
  ```
- If decontamination is enabled, an accessible Kraken2 database. The default path is genotoul cluster-specific and can be overridden with `--kraken_db`.

## Input data

- **Reads**: a directory containing paired `.fastq.gz` files. The pair names must follow the `1`/`2` convention matched by the workflow, for example:
  ```
  clean-fastq/
  ├── ANIC00017_0101.1.fastq.gz
  ├── ANIC00017_0101.2.fastq.gz
  ├── GDEL00541_0101.1.fastq.gz
  └── GDEL00541_0101.2.fastq.gz
  ```
- **Probes**: a fasta file containing the full set of UCE probes (e.g. `probes_hymeno_assembled_wrap.fasta`).
- **Existing contigs (optional)**: when using `--skip_assembly true`, provide a directory containing one or more `.fasta` files with the contigs for each sample.



> **Note:** this pipeline does not trim input reads. Adapter and quality trimming must be performed beforehand (for example, with Trimmomatic or Illumiprocessor).

## Installation and usage

```bash
git clone https://github.com/Madeline38/extracting_UCEs_nextflow.git
cd extracting_UCEs_nextflow
nextflow run main.nf --reads [raw_data] --probes [probes_file] [OPTION]
```


## Parameters (`nextflow.config`)

| Parameter | Default value | Description |
|---|---|---|
| `reads` | `none` | Directory containing paired, compressed FASTQ reads. Required unless `--skip_assembly true` is used. |
| `skip_decontamination` | `false` | Set to `true` to skip human-read decontamination when starting from reads. |
| `kraken_db` | `'/bank/kraken/k2_pluspfp_20251015'` | Path to the Kraken2 database used for decontamination. This default is environment-specific. |
| `keep_kraken_output` | `false` | Keep the large intermediate files produced by Kraken2 during the decontamination step. Enabling this preserves more cached data, which speeds up reruns with -resume, but requires significant disk space. By default, these files are deleted. |
| `skip_assembly` | `false` | Set to `true` to skip read processing and assembly and use existing contigs instead. |
| `contigs` | `none` | Directory containing input contig FASTA files. Required when `--skip_assembly true` is used; the workflow looks for `*.fasta` files. |
| `type_assembly` | `megahit` | Assembler to use: `megahit` or `spades`. |
| `busco_db` | `hymenoptera_odb12` | BUSCO lineage database used to assess assembly completeness. |
| `probes` | `none` | FASTA file of UCE probes. Required in both execution modes. |
| `flank` | `160` | Number of bases to extract on each side of a matched UCE. |
| `min_percent` | `75` | Minimum percentage of samples in which a locus must be present to retain it (0–100). At least two samples are required for a locus to be retained. |
| `trimming` | `none` | Optional alignment trimming. Use `trimal` to enable TrimAl; default is no trimming. |
| `phylogenetic_format` | `phylip` | Output alignment format: `nexus` or `phylip`. |
| `outdir` | `results` | Output directory for published files and pipeline reports. |



Parameters can be overridden on the command line. For paired reads (the default mode):

```bash
nextflow run main.nf \
  --reads clean-fastq \
  --probes probes_hymeno_assembled_wrap.fasta \
  --outdir results
```

To start from existing contigs instead of reads and assembly, provide a directory containing `.fasta` files:

```bash
nextflow run main.nf \
  --skip_assembly true \
  --contigs contigs_dir \
  --probes probes_hymeno_assembled_wrap.fasta \
  --min_percent 50 \
  --phylogenetic_format phylip \
  --outdir results
```

To resume an interrupted run or reuse valid cached tasks after changing downstream steps, repeat the original command and add `-resume`, for example:

```bash
nextflow run main.nf \
  --reads clean-fastq \
  --probes probes_hymeno_assembled_wrap.fasta \
  --outdir results \
  -resume
```

## Outputs

```
results/
├── Etape_1_analyses/                      # FastQC and MultiQC before decontamination
├── Etape_2_decontamination/               # Kraken2 reports, cleaned reads, and summary
├── Etape_3_analyses_post_decontamination/ # FastQC and MultiQC after decontamination
├── Etape_4_assemblage/                    # contigs and assembly logs
├── Etape_5_analyses_assemblage/           # QUAST and BUSCO results
├── Etape_6_alignment_probes_contigs/      # LASTZ probe-matching results
├── Etape_7_remove_paralogs/               # filtered SAM files and paralogs reports
├── Etape_8_alignements/                   # flanked sequences, loci, and MAFFT alignments
├── Etape_9_final_files/                   # UCE summary and Nexus/Phylip output
└── pipeline_info/                         # Nextflow report, timeline, and trace
```

Some directories are only populated when their corresponding steps run (for example, decontamination outputs when `--skip_decontamination true` is not set).

## Notes from previous runs

- The number of reads removed during decontamination depends on the Kraken2 database and the samples.
- SPAdes took substantially longer than MEGAHIT. Runtime and assembly quality depend on the dataset and available resources; benchmark both assemblers for your data if needed.
- Alignment trimming is optional (`--trimming trimal`). Compare trimmed and untrimmed alignments to decide whether trimming is appropriate for your data.
- Nextflow can reuse cached tasks with `-resume` when task inputs and outputs are available. The decontamination process removes its Kraken2 output file after read extraction, so this intermediate file is not retained. Keeping it would require modifying the decontamination process and would consume additional disk space.

## Resources

As an example, a previous run with 21 shotgun samples, decontamination and MEGAHIT assembly took about 17 hours 42 minutes and used approximately 500 GB of disk space. These figures are dataset- and cluster-specific; actual resource requirements may differ substantially.

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





