# Protein Protein Interaction - Alternative Splicing

## Introduction

This work pretends to create a program to perform some different analysis once we have the results of an RNA-Seq analysis done in different tissues and samples.

The main goal for doing this is to be able to determine differences in the expression of the same gene with different isoforms expressed in different tissues.

Initially, we have the results of an RNA-Seq analysis of different tisues: the expression of each isoform of each sample inside one tissue file (one *.csv.gz* or *.csv* file per tissue). This is an example of the file needed:

| transcript_id | gene_name | gene_type | transcript_type | gene_id | Samp1 | Samp2 | Samp3 | ... |
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| ENST00000373020.4 | TSPAN6 | protein_coding | protein_coding | ENSG00000000003.10 | 0.71 | 0.94 | 0.77 | ... |
| ENST00000494424.1 | TSPAN6 | protein_coding | processed_transcript | ENSG00000000003.10 | 0.08 | 0.00 | 0.08 | ... |
| ENST00000494424.1 | TSPAN6 | protein_coding | processed_transcript | ENSG00000000003.10 | 0.21 | 0.06 | 0.15 | ... |
| ENST00000373031.4 | TNMD | protein_coding | protein_coding | ENSG00000000005.5 | -1 | -1 | -1 | ... |
| ENST00000485971.1 | TNMD | protein_coding | processed_transcript | ENSG00000000005.5 | -1 | -1 | -1 | ... |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |
| ENST00000359637.2 | CFH | protein_coding | protein_coding | ENSG00000000971.11 | 0.09 | 0.09 | 0.08 | ... |
| ENST00000367429.4 | CFH | protein_coding | protein_coding | ENSG00000000971.11 | 0.63 | 0.57 | 0.53 | ... |
| ENST00000439155.2 | CFH | protein_coding | protein_coding | ENSG00000000971.11 | 0.21 | 0.34 | 0.25 | ... |
| ENST00000466229.1 | CFH | protein_coding | retained_intron | ENSG00000000971.11 | 0.04 | 0.00 | 0.03 | ... |
| ENST00000470918.1 | CFH | protein_coding | retained_intron | ENSG00000000971.11 | 0.03 | 0.00 | 0.08 | ... |
| ENST00000496761.1 | CFH | protein_coding | processed_transcript | ENSG00000000971.11 | 0.00 | 0.00 | 0.03 | ... |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |

You may see that there are some gene with -1 as values. When this happens, it means that that gene is not expressed in that sample (but it could be expressed in another sample).

Once we have this, we can perform an analysis of each gene, each tissues or all of theme. So, in order to do this, the project has one main Python file:

> Functions.py

and three more files that call Functions.py to perform each analysis:

> gene_analysis.py
> one_tissue_analysis.py
> all_tissues_analysis.py