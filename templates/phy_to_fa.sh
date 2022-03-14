#!/usr/bin/env bash

#name=!{params.prefix}
# convert from phylip to fasta aln for beauti
sed '/_/s/^/>/' !{phy_msa} | sed '1d' | fold -w 60 -s > !{params.prefix}_msa.fasta

#NOTE: parameter expansion of input/output vars defined in nextflow
#process is !{} instead of usual ${}

#-i '' '1d'
