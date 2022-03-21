#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process msa {

  tag "aligning sequences from ${multi_fasta} with MAFFT"
  publishDir "${params.out_dir}/${params.prefix}/", mode: 'copy'

  input:
  file(multi_fasta)

  output:
  path("*${params.prefix}*_msa_raw.fasta")

  """
  mafft --auto ${multi_fasta} > ${params.prefix}_msa_raw.fasta
  """
}

process msa_validation {

  tag "deduplicating & checking readability of ${msa}"
  publishDir "${params.out_dir}/${params.prefix}/", mode: 'copy'

    input:
    file(msa)

    output:
    file("*${params.prefix}*.phy")

    """
    raxml-ng \
    --check  \
    --msa ${msa} \
    --model ${params.evo_model} \
    --prefix ${params.prefix}
    """
}


workflow MSA {

  //Define channel for input multifasta file
  multi_ch = Channel.fromPath(params.multi_fa, checkIfExists:true)

  msa(multi_ch) | msa_validation

  msa_out = msa_validation.out

  emit:
  msa_out

}
