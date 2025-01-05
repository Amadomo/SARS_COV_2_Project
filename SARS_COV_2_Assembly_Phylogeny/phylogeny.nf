#!/usr/bin/env nextflow

// Définir les paramètres
params.mafft = "mafft"
params.iqtree = "iqtree2"
params.treetime = "treetime"
params.reference_tree = "/path/to/reference_tree.nwk"

// Fichiers d'entrée
params.fasta = "data/*.fasta"

workflow {
    // Étape 1 : Alignment with MAFFT
    Channel
        .fromFilePairs(params.fasta, flat: true)
        .set { fasta_files }

    fasta_files | processAlignment | processPhylogeny | processTemporalAnalysis
}

process processAlignment {
    tag "$sample"
    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}.aligned.fasta")

    script:
    """
    ${params.mafft} --auto ${fasta} > ${sample}.aligned.fasta
    """
}

process processPhylogeny {
    tag "$sample"
    input:
    tuple val(sample), path(aligned_fasta)

    output:
    tuple val(sample), path("${sample}_ML_tree.nwk"), path("${sample}_model.txt")

    script:
    """
    ${params.iqtree} -s ${aligned_fasta} -m MF -nt AUTO -bb 1000 -alrt 1000 -pre ${sample}_ML_tree
    """
}

process processTemporalAnalysis {
    tag "$sample"
    input:
    tuple val(sample), path(ml_tree), path(model_file)

    output:
    path("${sample}_temporal_tree.nwk")

    script:
    """
    ${params.treetime} clock --tree ${ml_tree} --alignment ${sample}.aligned.fasta --outdir temporal_analysis --model ${model_file} --gtr --coalescent --reference ${params.reference_tree}
    """
}

