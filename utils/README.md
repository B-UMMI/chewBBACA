# TestGenomeQuality.py

Evaluate wgMLST call quality per genome

Usage:

	% chewBBACA.py TestGenomeQuality -i out.txt -n 12 -t 250
	
`-i` raw output file from an allele calling

`-n` maximum number of iterations, each iteration removes a set of genomes over the threshold and recalculates all variables

`-t` maximum threshold, will start at 5 increasing in a step of 5 until t

The output consists on an interactive plot and a `removedGenomes.txt`. The file contains the information of which genomes were removed per threshold.

Example of the plot output can be seen [here](http://im.fm.ul.pt/chewBBACA/GenomeQual/GenomeQualityPlot_all_genomes.html). This examples use an 
original set of 710 genomes and a scheme of 1264 loci, using a parameter `-n 12`,`-s 5` and `-t 300`.

Based on the plot and specific knowledge of the species under analysis the user has to decide which threshold he wants to use and what, if any, genomes should be excluded the analysis due to possible draft genome assembly/annotation (i.e. CDS definition) problems.

# init_schema_4_bbaca.py

chewBBACA allele call requires a "short" version of each gene file, which regular public schemas d'ont provide. This script creates those files in order to prepare external schemas to be run with chewBBACCA.

This program prepares a schema for a chewBBACA allele call, removing alleles that are not a cds, creating a short version of each fast with only the 1st cds allele and removing loci that have no cds alleles.

	% init_schema_4_bbaca.py -i schema/
	
`-i` path to folder

`--cpu` number of cpu to use. Default 1.

`-v` (optional)more verbose output

# AutoAlleleCuration.py

Remove alleles from a specific list of genomes

	% AutoAlleleCuration.py -i listgenes.txt -g listgenomes.txt
	
`-i` list of path of gene files to process

`-g` list of genomes names from which the alleles must be removed


# CountNumberMissingData.py

Program that prints the number of missing loci from a raw allele call output file, per genome

	% CountNumberMissingData.py -i raw_genes.tsv
	
`-i` raw allele call output file


# RemoveGenes.py

This program removes genes from a tab separated allele profile file. Writes a file (oh the path the program is being run) called new.tsv .

	% RemoveGenes.py -i raw_output.tsv -g listGenes2Remove.txt -o newFile
	
`-i` main matrix file from which to remove

`-g` list of genes to remove

`-o` output file name

`--inverse` (optional) instead of removing the list of genes, keep the ones on the list and remove the others

# RemoveGenomes.py

This program removes genomes from a tab separated allele profile file

	% RemoveGenomes.py -i raw_output.tsv -l listgenomes.txt -o new_file.tsv
	
`-i` main matrix file from which to remove

`-l` list of genomes to remove

`-o` output file name

