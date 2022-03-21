#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/**
process tree_model {

    input:
    file(validated_msa)

    output:
    tuple val(model), file("${params.prefix}*.") into tree_ch

    """
    iqtree -s ${validated_msa} --prefix ${params.prefix} -m MF -T AUTO
    #bash script to retrieve model from log file
    """
}
*/

//-m MFP finds best model and completes analysis using it.
//can accept phy and fa msa format
process construct_tree {

    tag "contruct maximum likelihood (ML) tree with IQ-TREE \
using ${params.evo_model} model."
    publishDir "${params.out_dir}/${params.prefix}/", mode: 'copy'

    input:
    file validated_msa

    output:
    file("*.treefile")

    script:
    """
    iqtree -s ${validated_msa} \
    --prefix ${params.prefix} \
    -m ${params.evo_model} \
    ${params.iqtree_args}
    """

}


workflow TREE {

  take:
  msa1

  main:
  construct_tree(msa1)

  construct_tree.out.subscribe onComplete: { log.info """\nCheck ${params.out_dir}/${params.prefix}.treefile
in TempEst for temporal signal & mistakenly named sequences.\n
Proceed when ready by adding < --run_part2 true -resume > to original command.\n""" }

  emit:
  construct_tree.out

}
