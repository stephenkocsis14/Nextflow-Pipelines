#!/usr/bin/env nextflow

params.reads = "$baseDir/data/*_{1,2}.fastq.gz"
params.outdir = "$baseDir/results"

process CreateOutputDir {
    script:
    """
    mkdir -p ${params.outdir}
    """
}

process FastQC {
    beforeScript 'module load fastqc'

    input:
    file fastq from params.reads

    output:
    file "${fastq.baseName}_fastqc.zip" into fastqc_results
    file "${fastq.baseName}_overrepresented.txt" into overrepresented_files

    script:
    """
    fastqc ${fastq} --outdir ${params.outdir}
    unzip ${params.outdir}/${fastq.baseName}_fastqc.zip -d ${params.outdir}
    grep '>>Overrepresented sequences' -A 1000 ${params.outdir}/${fastq.baseName}_fastqc/fastqc_data.txt | grep -v '>>END_MODULE' | grep -v '>>Overrepresented sequences' > ${params.outdir}/${fastq.baseName}_overrepresented.txt
    """
}

process CollectOverrepresentedSequences {
    input:
    file overrepresented_files from fastqc_results.collect()

    output:
    file "combined_overrepresented_sequences.txt" into combined_sequences

    script:
    """
    cat ${overrepresented_files} > ${params.outdir}/combined_overrepresented_sequences.txt
    """
}

process MultiQC {
    beforeScript 'module load multiqc'

    input:
    file zipfiles from fastqc_results.collect()

    output:
    file "multiqc_report.html" into multiqc_report

    script:
    """
    multiqc ${params.outdir} -o ${params.outdir}
    """
}

process Cutadapt {
    beforeScript 'module load cutadapt'

    input:
    file fastq from params.reads
    file overrepresented from combined_sequences

    output:
    file "trimmed_${fastq.baseName}.fastq.gz" into trimmed_reads

    script:
    """
    adapters=$(grep -v '#' ${overrepresented} | awk '{print $1}' | paste -sd ',' -)
    cutadapt -a ${adapters} -A ${adapters} -o ${params.outdir}/trimmed_${fastq.baseName}.fastq.gz ${fastq}
    """
}

process STAR {
    beforeScript 'module load star'

    input:
    file trimmed_fastq from trimmed_reads.collect()

    output:
    file "${trimmed_fastq.baseName}.bam" into aligned_reads

    script:
    """
    STAR --genomeDir genome/STAR_index --readFilesIn ${trimmed_fastq.join(' ')} --runThreadN 4 --outFileNamePrefix ${params.outdir}/${trimmed_fastq.baseName}_ --outSAMtype BAM SortedByCoordinate
    """
}

process FeatureCounts {
    beforeScript 'module load subread'

    input:
    file bam from aligned_reads.collect()

    output:
    file "counts_${bam.baseName}.txt" into counts_table

    script:
    """
    featureCounts -a genome/annotation.gtf -o ${params.outdir}/counts_${bam.baseName}.txt ${bam.join(' ')}
    """
}

process EdgeR {
    beforeScript 'module load R'

    input:
    file counts from counts_table.collect()

    output:
    file "edgeR_results_*.csv" into edgeR_results

    script:
    """
    Rscript -e "
    library(edgeR)
    file.names <- Sys.glob('${params.outdir}/counts_*.txt')
    count.tables <- lapply(file.names, read.table, header=TRUE, row.names=1)
    sample.names <- gsub('counts_|.txt', '', file.names)
    counts <- do.call(cbind, count.tables)
    colnames(counts) <- sample.names
    unique_groups <- unique(gsub('_[^_]+$', '', sample.names))
    for (group in unique_groups) {
        if (sum(grepl(group, sample.names)) > 1) {
            comparison <- gsub('_[^_]+$', '', sample.names) == group
            y <- DGEList(counts=counts[, comparison])
            y <- calcNormFactors(y)
            design <- model.matrix(~ 0 + group)
            colnames(design) <- levels(group)
            y <- estimateDisp(y, design)
            fit <- glmFit(y, design)
            for (i in 2:length(levels(group))) {
                contrast <- makeContrasts(levels(group)[i] - levels(group)[1], levels=design)
                lrt <- glmLRT(fit, contrast=contrast)
                write.csv(topTags(lrt, n=Inf), paste0('${params.outdir}/edgeR_results_', levels(group)[i], '_vs_', levels(group)[1], '.csv'))
            }
        }
    }
    "
    """
}

workflow {
    create_output_dir = CreateOutputDir()
    fastqc_results = FastQC()
    combined_sequences = CollectOverrepresentedSequences(fastqc_results)
    multiqc_report = MultiQC(fastqc_results)
    trimmed_reads = Cutadapt(combined_sequences)
    aligned_reads = STAR(trimmed_reads)
    counts_table = FeatureCounts(aligned_reads)
    edgeR_results = EdgeR(counts_table)
}
