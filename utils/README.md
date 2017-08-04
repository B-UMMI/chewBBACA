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

# RemoveGenomes.py

This program removes genomes from a tab separated allele profile file

	% RemoveGenomes.py -i raw_output.tsv -l listgenomes.txt -o new_file.tsv
	
`-i` main matrix file from which to remove

`-l` list of genomes to remove

`-o` output file name

# TestGenomeQuality.py

Evaluate wgMLST call quality per genome

Usefull to determine a core genome and remove genomes that may have technical issues.

	% TestGenomeQuality.py -i out.txt -n 12 -t 250
	
`-i` raw output file from an allele calling

`-n` maximum number of iterations, each iteration removes a set of genomes over the threshold and recalculates all variables

`-t` maximum threshold, will start at 5 increasing in a step of 5 until t

The output consists in a set of plots per iteration and a removedGenomes.txt file where its informed of which genomes are removed per threshold when it reaches a stable point (no more genomes are removed)


# init_schema_4_bbaca.py

chewBBACA allele call requires a "short" version of each gene file, which regular public schemas d'ont provide. This script creates those files in order to prepare external schemas to be run with chewBBACCA.

This program prepares a schema for a chewBBACA allele call, creating a short version of each fast with only the 1st allele.

	% init_schema_4_bbaca.py -i listgenes.txt
	
`-i` list of path of gene files to process
