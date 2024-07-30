# RNA-seq analayis pipeline script

Structure the directory as follows:

```
rna-seq-pipeline/
│
├── data/
│   ├── sample1_R1.fastq.gz
│   ├── sample1_R2.fastq.gz
│   ├── sample2_R1.fastq.gz
│   ├── sample2_R2.fastq.gz
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
