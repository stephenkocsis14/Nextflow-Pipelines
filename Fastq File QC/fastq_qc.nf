// fastq_qc.nf

// Define parameters with defaults
params.fastq_dir = './path_to_fastq_files/'    // Default FASTQ file directory
params.output_dir = './qc_results/'                 // Default output directory
params.zipped = true                                // Default to TRUE, expecting .fastq.gz files

// Set file extension based on zipped parameter
def fastq_extension = params.zipped ? "*.fastq.gz" : "*.fastq"

// Create channels for FASTQ files based on the extension
Channel
    .fromPath("${params.fastq_dir}/${fastq_extension}")
    .set { fastq_files }

// Process to run FastQC
process FastQC {
    input:
    path fastq_file from fastq_files

    output:
    path "*.zip" into fastqc_zips
    path "*.html" into fastqc_reports

    script:
    """
    mkdir -p ${params.output_dir}/FastQC
    fastqc $fastq_file --outdir ${params.output_dir}/FastQC
    """
}

// Process to run MultiQC
process MultiQC {
    input:
    path fastqc_zips.collect()

    output:
    path "*.html" into multiqc_report

    script:
    """
    mkdir -p ${params.output_dir}/MultiQC
    multiqc ${params.output_dir}/FastQC -o ${params.output_dir}/MultiQC
    """
}

// Workflow definition
workflow {
    fastqc_zips = FastQC()
    multiqc_report = MultiQC(fastqc_zips)
}
