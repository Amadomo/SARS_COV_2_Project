#!/usr/bin/env nextflow

// Define input files
params.fastq = "path to fastq files"
params.index = "Path to reference index files"
params.reference = "Path to reference sequence"
params.adapters = "Path to adapter sequence file"
params.min_quality = // define the min quality 
params.consensus_threshold = // define the consensus threshold 
params.min_depth = // define min depth 
params.consensus = "path where should consensus sequence directories"

workflow {
    fastq_ch = Channel.fromPath(params.fastq)
        .map { file ->
            def sample = file.baseName
            [sample, file]
        }

    // Pass the formatted tuples to the quality control process
    trimmed_ch = processQualityControl(fastq_ch)
    trimmed_ch.view { "Trimmed: ${it}" }

    // Continue with the rest of the pipeline
    ref_index_ch = Channel.fromPath(params.index)
    reference_ch = Channel.fromPath(params.reference)

    aligned_ch = trimmed_ch
    .combine(ref_index_ch)
    .combine(reference_ch)
    .map { sample, trimmed_reads, ref_index, reference ->
        [sample, trimmed_reads, ref_index, reference]
    }

    aligned_ch = processAlignment(aligned_ch)
    aligned_ch.view { "Aligned: ${it}" }
    
    consensus_ch = processConsensus(aligned_ch)
    filtered_ch = processFilterConsensus(consensus_ch)

    alignment2_ch = processAlignment2(filtered_ch)
    phylogeny_ch = processPhylogeny(alignment2_ch)
    //temporal_analysis_ch = processTemporalAnalysis(phylogeny_ch)
}

process processQualityControl {

    publishDir("${params.control}", mode: 'copy')

    tag "$sample"
    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}.trimmed.fq.gz"), emit: trimmed_reads

    script:
    """
    echo "Performing QC for $sample"
    trimmomatic SE -threads 4 -phred33 ${reads} \
        ${sample}.trimmed.fq.gz \
        ILLUMINACLIP:${params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    """
}

process processAlignment {

    publishDir("${params.aligned}", mode: 'copy')

    input:
    tuple val(sample), path(trimmed_reads), path(ref_index), path(reference)

    output:
    tuple val(sample), path("${sample}.sorted.bam"), emit: bam

    script:
    """
    echo "Aligning $sample"
    echo "Input file: ${trimmed_reads}"
    bwa mem -t 4 ${ref_index}/${reference} ${trimmed_reads} | \
        samtools view -Sb -o ${sample}.bam -
    samtools sort -o ${sample}.sorted.bam ${sample}.bam
    samtools index ${sample}.sorted.bam
    """
}

process processConsensus {

    publishDir("${params.consensus_nf}", mode: 'copy')

    input:
    tuple val(sample), path(aligned_reads)

    output:
    tuple val(sample), path("${sample}.consensus.fa"), emit: consensus

    script:
    """
    echo "Generating consensus for $sample"
    samtools mpileup -aa -A -d 0 -Q ${params.min_quality} -B -f ${params.reference} ${aligned_reads} | \
    ivar consensus -t ${params.consensus_threshold} -m ${params.min_depth} -p ${sample}.consensus
    """
}

process processFilterConsensus {

    publishDir("${params.consensus}", mode: 'copy')

    input:
    tuple val(sample), path(consensus_sequences)

    output:
    tuple val(sample), path("${sample}.filtered.fa"), emit: filtered

    script:
    """
    echo "Filtering consensus for $sample"
    if grep -q "N" ${consensus_sequences}; then
        missing_bases=\$(grep -o "N" ${consensus_sequences} | wc -l)
        total_bases=\$(grep -v ">" ${consensus_sequences} | wc -c)
        missing_ratio=\$(echo "\${missing_bases} / \${total_bases}" | bc -l)
        if (( \$(echo "\${missing_ratio} >= 0.1" | bc -l) )); then
            echo "Skipping ${sample} due to >10% missing bases" > ${sample}.log
            echo ">${sample}_skipped" > ${sample}.filtered.fa  # Create a placeholder file
        else
            cp ${consensus_sequences} ${sample}.filtered.fa
        fi
    else
        cp ${consensus_sequences} ${sample}.filtered.fa
    fi
    """
}

process processAlignment2 {
    tag "$sample"
    input:
    tuple val(sample), path(filtered_sequences)

    output:
    tuple val(sample), path("${sample}.aligned2.fasta"), emit: aligned2

    script:
    """
    mafft --auto ${filtered_sequences} > ${sample}.aligned2.fasta
    """
}

process processPhylogeny {
    tag "$sample"
    input:
    tuple val(sample), path(aligned2_fasta)

    output:
    tuple val(sample), path("${sample}_ML_tree.nwk"), path("${sample}_model.txt")

    script:
    """
    iqtree -s ${aligned2_fasta} -m MF -nt AUTO -pre ${sample}_ML_tree
    """
}
