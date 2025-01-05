#!/usr/bin/env nextflow

// Définir les paramètres pour les outils
params.trimmomatic = "path to trimmomatic"
params.adapters = "path to adapters file"
params.reference = "path to reference sequence file"
params.bwa = "paht to bwa"
params.ivar = "path to ivar"
params.min_depth = 10
params.min_quality = 20
params.consensus_threshold = 0.0

// Définir les fichiers d'entrée
params.fastq = "path to read fastq file "

workflow {
    // Étape 1: Quality control and trimming 
    Channel
        .fromFilePairs(params.fastq, flat: true)
        .set { fastq_files }

    fastq_files | processQualityControl | processAlignment | processConsensus | processFilterConsensus
}

process processQualityControl {
    tag "$sample"
    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("*.trimmed.fq.gz")

    script:
    """
    java -jar ${params.trimmomatic} PE \
        -threads 4 \
        ${reads[0]} ${reads[1]} \
        ${sample}_paired_R1.trimmed.fastq.gz ${sample}_unpaired_R1.trimmed.fastq.gz \
        ${sample}_paired_R2.trimmed.fastq.gz ${sample}_unpaired_R2.trimmed.fastq.gz \
        ILLUMINACLIP:${params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    """
}

process processAlignment {
    tag "$sample"
    input:
    tuple val(sample), path(trimmed_reads)

    output:
    tuple val(sample), path("${sample}.bam")

    script:
    """
    ${params.bwa} mem -t 4 ${params.reference} ${trimmed_reads[0]} ${trimmed_reads[1]} | \
        samtools view -Sb - | \
        samtools sort -o ${sample}.bam
    samtools index ${sample}.bam
    """
}

process processConsensus {
    tag "$sample"
    input:
    tuple val(sample), path(bam_file)

    output:
    tuple val(sample), path("${sample}.consensus.fa")

    script:
    """
    samtools mpileup -aa -A -d 0 -Q ${params.min_quality} -B -f ${params.reference} ${bam_file} | \
        ${params.ivar} consensus -t ${params.consensus_threshold} -m ${params.min_depth} -p ${sample}.consensus
    """
}

process processFilterConsensus {
    tag "$sample"
    input:
    tuple val(sample), path(consensus_file)

    output:
    path("${sample}.final.fa")

    script:
    """
    if grep -q "N" ${consensus_file}; then
        missing_bases=\$(grep -o "N" ${consensus_file} | wc -l)
        total_bases=\$(grep -v ">" ${consensus_file} | wc -c)
        missing_ratio=\$(echo "\${missing_bases} / \${total_bases}" | bc -l)
        if (( \$(echo "\${missing_ratio} >= 0.1" | bc -l) )); then
            echo "Skipping ${sample} due to >10% missing bases" > ${sample}.log
        else
            cp ${consensus_file} ${sample}.final.fa
        fi
    fi
    """
}

