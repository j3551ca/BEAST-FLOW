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
--xml_args                     Input parameters for generating xml for BEAST. Use <beast2-xml.py --help>
                               for full list of options [--chainLength 1000000 --sequenceIdAgeRegex...]
--xml_temp                     Template used to generate xml, differs depending on clock used [strict.xml]
--beast_args                   Parameters for BEAST run. Use <beast -help> for all options [-overwrite]
--treeann_args                 TreeAnnotator arguments. Use <treeannotator -help> for options [-burnin 10]
--iqtree_args                  IQ-TREE arguments. Use <iqtree --help> for more options [-T 4 -B 1000]
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
import processes from modules
---------------------------------------------------------------------------------
*/

include { MSA } from "${projectDir}/modules/msa.nf"
include { TREE } from "${projectDir}/modules/tree.nf"
include { BEAST } from "${projectDir}/modules/beast.nf"

workflow {

MSA | TREE
BEAST(MSA.out)

}
