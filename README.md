![image](/pics/beast-flow_logo_outline.png)

## Introduction

A Nextflow pipeline for molecular clock analysis using Bayesian Evolutionary Analysis Sampling Trees (BEAST v2.6.6).

BEAST-FLOW automates the process of estimating molecular evolutionary rate where possible, 
while integrating breaks for manual revision of data-suitability in TempEst and Tracer before proceeding. The pipeline mandatorily accepts a multi-fasta file of various sample consensus sequences and a prefix string as its input. The pipeline uses MAFFT, IQ-TREE, BEAST2 XML, BEAST2, and TreeAnnotator. 

BEAST-FLOW is written in the glue language and workflow engine of Nextflow. The parallelization, portability, and modularity of dataflow programming facilitate high-throughput data analysis, reproducibility, and customization. Modifications to the default parameters can be easily adjusted to fit the purpose of a project.  

BEAST-FLOW was first published as part of a paper (manuscript in preparation) classifying SARS-CoV-2 reinfections. 

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

1. Bouckaert, R., Heled, J., Kühnert, D., Vaughan, T., Wu, C-H., Xie, D., Suchard, MA., Rambaut, A., & Drummond, A. J. (2014). BEAST 2: A Software Platform for Bayesian Evolutionary Analysis. PLoS Computational Biology, 10(4), e1003537. doi:10.1371/journal.pcbi.1003537
2. Drummond, A. J. and A. Rambaut. 2007. BEAST: Bayesian evolutionary analysis by sampling trees. BMC Evolutionary Biology 7:214.
3. Jones, T. (2018). BEAST2 XML. https://github.com/acorg/beast2-xml
4. Katoh, K., Misawa, K., Kuma, K., & Miyata, T. (2002). MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic Acids Research, 30(14), 3059–3066. https://doi.org/10.1093/nar/gkf436
5. Kozlov AM, Darriba D, Flouri T, Morel B, Stamatakis A. RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. Bioinformatics. 2019 Nov 1;35(21):4453-4455. doi: 10.1093/bioinformatics/btz305. PMID: 31070718; PMCID: PMC6821337. 
6. Nguyen, L.-T., Schmidt, H. A., Haeseler, A. von, & Minh, B. Q. (2014). IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies. Molecular Biology and Evolution, 32(1), 268–274. https://doi.org/10.1093/molbev/msu300
