[Image](/pics/beast-flow_logo.png)

## Introduction

A Nextflow pipeline for molecular clock analysis using Bayesian Evolutionary Analysis Sampling Trees (BEAST v2.6.6).

BEAST-FLOW automates the process of estimating molecular evolutionary rate where possible, 
while integrating breaks for manual revision of data-suitability in TempEst and Tracer before proceeding. The pipeline mandatorily accepts a multi-fasta file of various sample consensus sequences and a prefix string as its input. The pipeline uses MAFFT, IQ-TREE, BEAST2 XML, BEAST2, and TreeAnnotator. 

BEAST-FLOW is written in the glue language and workflow engine of Nextflow. The parallelization, portability, and modularity of dataflow programming facilitate high-throughput data analysis, reproducibility, and customization. Modifications to the default parameters can be easily adjusted to fit the purpose of a project.  

BEAST-FLOW was first published as part of a paper classifying SARS-CoV-2 reinfections. 

## Table of Contents
- [Introduction](#introduction)
- [Workflow](#workflow)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)
- [Input](#input)
- [Output](#output)
- [References](#references)

## Workflow

![image](/pics/beast-flow_dag.png)

## Dependencies 

This bioinformatic pipeline requires Nextflow:

```
conda install -c bioconda nextflow
```

or download and add the nextflow executable to a location in your user $PATH variable:

```
curl -fsSL get.nextflow.io | bash
mv nextflow ~/bin/
```

Nextflow requires Java v8.0+, so check it is installed:

```
java -version
```

All other dependencies, including MAFFT, BEAST2-XML, BEAST2, & IQ-TREE, can be found in the ‘beastflow_env.yml’ file and are activated upon running the program.

## Installation

To copy the program into a directory of your choice, from desired directory run:

```
git clone https://github.com/j3551ca/BEAST-FLOW.git
cd BEAST-FLOW
nextflow run main.nf -profile conda
```

or run directly using:

```
nextflow run j3551ca/BEAST-FLOW -profile conda
```

## Usage

Change into working directory:
```
cd /home/user/directory/containing/BEAST-FLOW
```
Run BEAST-FLOW pipeline:
```
nextflow run main.nf -profile conda --multi_fa dengue_multi.fasta --prefix  dengue_run1 [OPTIONS]
```
For BEAST-FLOW help message:
```
nextflow run main.nf --help
```

## Input

1.	Multi-fasta file with at least 5 sequences \(MAFFT will throw error otherwise\): 

\>Seq1_2021-05-31\
ATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTTTACAGGTTCGCGACGTGCTCGTACGTGGCTTTGGAGACTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATCAA\
\>Seq2_2021-09-16\
ATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTTTACAGGTTCGCGACGTGACCTTATGAAGATTTTCAAGAAAACTGGAACACTAAACATAGCAGTGGTGTTACCCGTGAACTCATGCGTGAGCT\
\>Seq3_2022-01-04\
ATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTTTACAGGTTCGCGACGTGCTCGTACGTGGCTTTGGAGACTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATCTTAAAGATGGCACTTGTGG\
.\
.\
.

2.	Prefix string. This label will be used to name output files generated in BEAST-FLOW and avoids overwriting:

dengue_run1

## Output

Depending on the options specified in the command line, the following result files should be placed in the output directory under a folder containing the same name as the prefix specified in the initial command (see step 2 under [Usage](#usage):

- A multiple sequence alignment \(MSA\): \*_msa.fasta
- A maximum-likelihood \(ML\) phylogenetic tree in Newick format: \*.treefile
- A trace log file from BEAST2: \*.log
- The posterior distribution of trees from BEAST2: \*.trees
- A maximum clade credibility \(MCC\) tree: \*_mcc.tree
 

## References


