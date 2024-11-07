# FASTQ Quality Control Pipeline

This Nextflow pipeline performs quality control on FASTQ files using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and aggregates the results with [MultiQC](https://multiqc.info/). The pipeline supports both gzipped (`.fastq.gz`) and uncompressed (`.fastq`) FASTQ files.

## Table of Contents
- [Requirements](#requirements)
- [Usage](#usage)
  - [Basic Command](#basic-command)
  - [Parameters](#parameters)
  - [Examples](#examples)
- [Output](#output)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## Requirements

- **Nextflow**: Install Nextflow by following the [installation instructions](https://www.nextflow.io/docs/latest/getstarted.html).
- **FastQC**: Ensure FastQC is installed and accessible in your PATH.
- **MultiQC**: Ensure MultiQC is installed and accessible in your PATH.

## Usage

### Basic Command

To run the pipeline, use the following command:
```bash
nextflow run fastq_qc.nf --fastq_dir '<path_to_fastq_files>' --output_dir '<output_directory>' --zipped <true|false>
```

### Parameters

- `--fastq_dir` (required): Path to the directory containing FASTQ files.
  - Example: `--fastq_dir '/path/to/fastq_files/'`
  
- `--output_dir` (optional): Path to the output directory where results will be stored. The output will contain two subdirectories: `FastQC` for individual FastQC reports and `MultiQC` for the consolidated MultiQC report.
  - Default: `./qc_results/`
  - Example: `--output_dir '/desired/output_directory/'`
  
- `--zipped` (optional): Specify whether the input FASTQ files are gzipped (`true` for `.fastq.gz`) or uncompressed (`false` for `.fastq`).
  - Default: `true`
  - Example: `--zipped false`
 
### Examples

#### Example 1: Running with gzipped FASTQ files and default output directory
```bash
nextflow run fastq_qc.nf --fastq_dir '/path/to/fastq_files/' --zipped true
```

#### Example 2: Running with uncompressed FASTQ files and custom output directory
```bash
nextflow run fastq_qc.nf --fastq_dir '/path/to/fastq_files/' --output_dir '/custom/output_directory/' --zipped false
```

### Output

The pipeline creates the following directory structure in the specified output directory:

<output_directory>/
├── FastQC/
│   ├── [FastQC output files: .zip, .html for each FASTQ file]
└── MultiQC/
    └── multiqc_report.html

- **FastQC**: Contains individual FastQC reports for each FASTQ file.
- **MultiQC**: Contains a single MultiQC report summarizing all FastQC results.

### License 

This pipeline is licensed under the MIT license.

### Acknowledgments

This pipeline was created using Nextflow for streamlining quality control of FASTQ files with FastQC and MultiQC. 

This `README.md` file provides a clear overview of the pipeline usage, parameter settings, example commands, expected output structure, and troubleshooting tips. Let me know if you'd like to add more details!

