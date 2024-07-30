# RNA-seq analayis pipeline script

This pipeline leverages widely used bioinformatics tools to process raw RNA-seq fastq files using the following tools:
* FastQC: Perform quality control
* MultiQC: Perform quality control
* Cutadapt: Trim poor quality and overrepresented sequence reads
* STAR: Align the fastq files to the reference genome
* FeatureCounts: Create mRNA read count matrix using sample prefixes for groups
* EdgeR: Differential expression analysis

Structure the directory as follows:

```
rna-seq-pipeline/
│
├── data/
│   ├── CTL1_R1.fastq.gz
│   ├── CTL1_R2.fastq.gz
│   ├── DIS1_R1.fastq.gz
│   ├── DIS1_R2.fastq.gz
│   └── ... (additional samples)
│
├── genome/
│   ├── STAR_index/  (STAR genome index files)
│   └── annotation.gtf  (GTF annotation file)
│
├── results/  (output directory, will be created during pipeline execution)
│
├── main.nf  (Nextflow script)
│
├── nextflow.config  (Nextflow configuration file)
│
└── Dockerfile  (Dockerfile for building the container - if using Docker)
```
