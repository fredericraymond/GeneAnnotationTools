---
title: "Readme"
author: "Rogia Kpanou"
date: "10 August, 2017"
output: 
    html_document:
       toc: true
       theme: cerulean
       highlight: pygments
---


#### CluBB
Identify groups of contiguous genes belonging to the same metabolic pathways.
The pathways used come from Brenda.

#### USAGE
A sample command line run might look like: 
```r
$python CluBB.py  inputfolder
```
To see a detailed description of all command line options, do:
```r
$python CluBB.py -h
```
#### REQUIEREMENTS
+ Python 2.7
+ numpy
+ xwlt
+ pickle
+ docopt
+ xlrd

#### FILES DESCRIPTION
The directory contains three files:

1. CluBB.py : It' s the main file.This file performs the following actions:
    + Step 1: Runs prodigal on input folder files if -p option is specified
    + Step 2: Parses prodigal results
    + Step 3: Runs cd-hit-2d to identify the input file proteins(gene) that are similar to Brenda enzymes with a percentage identity of 90
    + Step 4: Genes found are grouped by pathways and by contig
    + Step 5: For each pathways and contig, only contiguous genes are considered to form cluster
    + Step 6: Calculates cluster density
2. Sanalysis.py : Contains functions used to perform previous actions
3. Brenda.py : This script is used to collect information about brenda enzymes and the pathways to which they belong

#### OUTPUT DESCRIPTION
Two types of files are obtained at the output:

1. .clstr file : Contains the list of all found genes groups, the name of the pathways to which they belong, the density of the pathways and that of the cluster
2. .xlsx : Contains statistics on pathways and brenda enzymes found in genomes.

NB: Pour des quelconques tests preliminaires, vous pouvez utiliser le fichier de sequences brenda: `/lustre4/CHUdeQuebec/nne-790-ad/projects/header_modif+NDn_BRENDA_sequences.fasta`