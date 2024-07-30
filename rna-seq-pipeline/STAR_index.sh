# Generate the STAR genome index given any reference fasta and gtf
STAR --runMode genomeGenerate --genomeDir genome/STAR_index --genomeFastaFiles /path/to/genome.fasta --sjdbGTFfile /path/to/annotation.gtf --runThreadN 4
