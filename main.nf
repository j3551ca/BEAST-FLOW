#!/usr/bin/env nextflow

//enable domain-specific-language 2
nextflow.enable.dsl=2

/**
---------------------------------------------------------------------------------
function definition
---------------------------------------------------------------------------------
*/

def helpMe() {
  log.info """

Overview:
Nextflow pipeline for estimation of evolutionary rate using Bayesian statistics.

Usage:
nextflow run main.nf [OPTIONS]

Mandatory arguments:
 --multi_fa                     Multi-fasta file containing consensus seqs of interest
 --prefix                       Prefix used for result files - avoids overwriting

Optional arguments:
--out_dir                      Output directory to place final BEAST-FLOW results [./output]
--evo_model                    Evolutionary model used to generate phylogenetic tree [GTR+G]
--run_part2                    Continue with rest of analysis after manual check of temporal signal [false]
--xml_args                     Input parameters for generating xml for BEAST [--clockModel strict --chainLength 1000000 ...]
--xml_temp                     Template used to generate xml, differs depending on clock used [strict.xml]
--beast_args                   Parameters for BEAST run [-overwrite]
--version                      Current BEAST-FLOW version number
--help                         This usage statement
        """
}

def version() {
  log.info """
  BEAST-FLOW version: ${workflow.manifest.version}
  """
}

//displays help upon request
if (params.help) {
  helpMe()
  exit 0 //stop running
}

//version upon request
if (params.version) {
  version()
  exit 0
}

/**
---------------------------------------------------------------------------------
program introduction
---------------------------------------------------------------------------------
*/

// this prints program header with mandatory input
log.info """


██████╗░███████╗░█████╗░░██████╗████████╗░░░░░░███████╗██╗░░░░░░█████╗░░██╗░░░░░░░██╗
██╔══██╗██╔════╝██╔══██╗██╔════╝╚══██╔══╝░░░░░░██╔════╝██║░░░░░██╔══██╗░██║░░██╗░░██║
██████╦╝█████╗░░███████║╚█████╗░░░░██║░░░█████╗█████╗░░██║░░░░░██║░░██║░╚██╗████╗██╔╝
██╔══██╗██╔══╝░░██╔══██║░╚═══██╗░░░██║░░░╚════╝██╔══╝░░██║░░░░░██║░░██║░░████╔═████║░
██████╦╝███████╗██║░░██║██████╔╝░░░██║░░░░░░░░░██║░░░░░███████╗╚█████╔╝░░╚██╔╝░╚██╔╝░
╚═════╝░╚══════╝╚═╝░░╚═╝╚═════╝░░░░╚═╝░░░░░░░░░╚═╝░░░░░╚══════╝░╚════╝░░░░╚═╝░░░╚═╝░░
\n=====================================================================================
multi-fasta: ${params.multi_fa}
prefix: ${params.prefix}
results folder: ${params.out_dir}
"""

/**
---------------------------------------------------------------------------------
processes
---------------------------------------------------------------------------------
*/

process msa {

  tag "aligning sequences from ${multi_fasta} with MAFFT"
  publishDir "${params.out_dir}/${params.prefix}/"

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
  publishDir "${params.out_dir}/${params.prefix}/"

    input:
    file(msa)

    output:
    file("*${params.prefix}*.phy")

    """
    raxml-ng --check  --msa ${msa} --model ${params.evo_model} --prefix ${params.prefix}
    """
}
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

//-m MFP finds best model and completes analysis using it
//can accept phy and fa msa format
process construct_tree {

    tag "contruct maximum likelihood (ML) tree with IQ-TREE \
using ${params.evo_model} model."
    publishDir "${params.out_dir}/${params.prefix}/"

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

//BEAUti does NOT accept PHY format
process phy_to_fa {

  tag "convert PHYLIP MSA to fasta for BEAUti2"
  publishDir "${params.out_dir}/${params.prefix}/"

  input:
  file(phy_msa)

  output:
  file("*${params.prefix}*_msa.fasta")

  shell: //use 'script:' when conditional script block or regular calling, incl.
  template "phy_to_fa.sh"

}

process xml_gen {

  tag "generate xml file for BEAST2"
  publishDir "${params.out_dir}/${params.prefix}/"

  input:
  tuple file(fa_msa), file(template)

  output:
  file("*${params.prefix}*.xml")

  when:
  params.run_part2 == true

  script://TODO: date parsing with regex not working properly
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
  publishDir "${params.out_dir}/${params.prefix}/"

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
  publishDir "${params.out_dir}/${params.prefix}/"

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

//view mol. clock in figtree
/**
---------------------------------------------------------------------------------
import processes from modules
---------------------------------------------------------------------------------
*/
/**like python "import pandas as pd"
direct main script to modules
*/
/**
include {
  msa;
  msa_validation
} from "${projectDir}/modules/msa.nf"
*/

/**
include { MSA } from "./modules/msa.nf"
include { TREE } from "${projectDir}/workflows/tree.nf"
include { BEAST } from "${projectDir}/workflows/beast.nf"

workflow {

MSA | TREE | BEAST

}
*/

/**
---------------------------------------------------------------------------------
main workflow commences here
---------------------------------------------------------------------------------
*/

workflow {

  //Define channel for input multifasta file
  multi_ch = Channel.fromPath(params.multi_fa, checkIfExists:true)

  msa(multi_ch) | msa_validation

  construct_tree(msa_validation.out)

  construct_tree.out.subscribe onComplete: { log.info """\nCheck ${params.out_dir}/${params.prefix}.treefile
in TempEst for temporal signal & mistakenly named sequences.\n
Proceed when ready by adding < --run_part2 true -resume > to original command.\n""" }

  if (msa_validation.out //channel
    .filter("*.phy")) {
      phy_to_fa(msa_validation.out)
      msa_ch = phy_to_fa.out
    }
  else {
    msa_ch = msa_validation.out
  }

  //tree_model(msa_ch) TODO

  //TempEst

  template_ch = Channel.fromPath(params.xml_temp, checkIfExists:true)

  xml_ch = msa_ch.combine(template_ch)

  xml_gen(xml_ch)

  beast(xml_gen.out)

  mcc_tree(beast.out)

  //FigTree

  emit: //emit channel from workflow to be accessed as name_of_workflow.out
  mcc_tree.out

}
