#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//BEAUti does NOT accept PHY format
process phy_to_fa {

  tag "convert PHYLIP MSA to fasta for BEAUti2"
  publishDir "${params.out_dir}/${params.prefix}/", mode: 'copy'

  input:
  file(phy_msa)

  output:
  file("*${params.prefix}*_msa.fasta")

  shell: //use 'script:' when conditional script block or regular calling, incl.
  template "phy_to_fa.sh"

}

process xml_gen {

  tag "generate xml file for BEAST2"
  publishDir "${params.out_dir}/${params.prefix}/", mode: 'copy'

  input:
  tuple file(fa_msa), file(template)

  output:
  file("*${params.prefix}*.xml")

  when:
  params.run_part2 == true

  script:
  """
  beast2-xml.py \
  --templateFile ${template}\
  ${params.xml_args} \
  --logFileBasename ${params.prefix}\
  --fastaFile ${fa_msa} > ${params.prefix}.xml
  """
}

process beast {

  tag "Bayesian evolutionary analysis sampling trees"
  publishDir "${params.out_dir}/${params.prefix}/", mode: 'copy'

  input:
  file xml

  output:
  file "*${params.prefix}*.trees"

  """
  beast ${params.beast_args} ${xml}
  """
}

process mcc_tree {

  tag "maximum credibility clade tree from BEAST2 output"
  publishDir "${params.out_dir}/${params.prefix}/", mode: 'copy'

  input:
  file(trees)

  output:
  file("*${params.prefix}*mcc.tree")

  """
  treeannotator \
  ${params.treeann_args} \
  ${trees} \
  ${params.prefix}_mcc.tree
  """
}

workflow BEAST {

  take:
  msa2 //needs to be renamed or ? not consumed

  main:
  //TempEst
  if (msa2
    .filter("*.phy")) {
      phy_to_fa(msa2)
      msa_ch = phy_to_fa.out
    }
  else {
    msa_ch = msa2
  }

  template_ch = Channel.fromPath(params.xml_temp, checkIfExists:true)

  xml_ch = msa_ch.combine(template_ch)

  xml_gen(xml_ch)

  beast(xml_gen.out)

  mcc_tree(beast.out)

  //FigTree

  emit: //emit channel from workflow to be accessed as name_of_workflow.out
  mcc_tree.out
  //view mol. clock in figtree
}
