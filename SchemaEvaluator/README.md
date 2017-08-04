# Validate Schema

Analyze your alleles for a set of parameters, taking special consideration on the allele CDS and the allele sizes per gene/locus

Python dependencies:
* numpy
* [Biopython 1.66 ](http://biopython.org/wiki/Main_Page)

Other dependencies:
* [mafft](http://mafft.cbrc.jp/alignment/software/linux.html)
* [clustalw2](http://www.clustal.org/download/2.1/)

The following visualization tools are later used on the html (not to install):
* [MSAViewer](http://msa.biojs.net/)
* [Phylocanvas](http://phylocanvas.org/)

# How to run

Main script that calls alleSizeStats.py and CheckCDS.py creating an html file for an easier results reading.

The optional parameters refer to the analysis of the alleles size. This analysis will calculate the mode size per gene and using that value -/+ a threshold (0.05 default) it will consider the allele as being a good allele if it is still considered within the threshold. It is given to the user the choice of the threshold and the choice to consider that a gene is considered "Conserved" if only one of the alleles is outside the threshold (default) or all the alleles must be within the threshold to be considered a "Conserved allele.

	% chewBBACA.py SchemaValidation -i ./genes/ -ta 11 -l ./rms/ratemyschema.html --cpu 3
		
`-i` directory where the genes .fasta files are located or alternatively a .txt file containing the full path for each gene .fasta file per line

`-ta` which translation table to use ( 11 in case of bacteria)

`--title` (optional) title to appear on the final html.

`-l` Location/name of the final html output

`--cpu` number of cpu to use, will be used for mafft and clustal

`--log` (optional) number of alleles per locus plot will be ploted in log10 scale

`-p` (optional) True if all alleles must be within threshold (default=False)

`-t` (optional) threshold used to calculate the range at which the allele is considered good (default=0.05)

`--light` (optional) skip mafft and clustal run, faster but less information on individual loci pages

### Two extra tab separated output files

An abridged example `locus_stats.tsv` file:

```
Locus	Mode_value	number_alleles
b0002.fasta	2463	1493
b0019.fasta	1167	906
b0007.fasta	1431	1198
b0009.fasta	588	524
b0006.fasta	777	711
b0020.fasta	906	672
b0010.fasta	567	529
b0008.fasta	954	752

```

An abridged example `non_cds_alleles.tsv` file:

```
Locus	Frameshift	No Start	More than 1 Stop	 Other
b0002.fasta	9	-	-	-
b0009.fasta	-	81,91	-	-
b0008.fasta	36	-	-	-

```
