# chewBBACA: Quick Usage

**chewBBACA** stands for "BSR-Based Allele Calling Algorithm". The "chew" part could be thought of as "Comprehensive and  Highly Efficient Workflow" 
but at this point still it needs a bit of work to make that claim so we just add "chew" to add extra coolness to the software name. BSR stands for 
BLAST Score Ratio as proposed by [Rasko DA et al.](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-2) 

chewBBACA is a comprehensive pipeline including a set of functions for the creation and validation of whole genome and core genome MultiLocus Sequence 
Typing (wg/cgMLST) schemas, providing an allele calling algorithm based on Blast Score Ratio that can be run in multiprocessor 
settings and a set of functions to visualize and validate allele variation in the loci.

----------
## Check the [wiki pages](https://github.com/mickaelsilva/chewBBACA/wiki) ...
...for a much more thorough chewBBACA walkthrough. 
Below you can find a list of commands for a quick usage of the software.

## Use [BBACA gitter](https://gitter.im/BBACA/Lobby)...

... if you have any pressing question. Chat can be faster and better than email for troubleshooting purposes 

**Important Notes before starting:**

 - **chewBBACA** define an allele as a complete Coding DNA Sequence, with start and stop codon according 
 to the [NCBI genetic code table 11](http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) identified using [Prodigal 2.6.0 ](https://github.com/hyattpd/prodigal/releases/). It will 
 automatically exclude any allele for which the DNA sequence does not contain start or stop codons and for which the length is not multiple of three. 
 - All the referenced lists of files *must contain full path* for the files.
 - Make sure that your fasta files are UNIX format. If they were created in Linux or MacOS systems they should be in the correct format, but if they were created in Windows systems, you should do a a quick conversion using for example [dos2unix](http://linuxcommand.org/man_pages/dos2unix1.html).

## An extensive [tutorial repository](https://github.com/mickaelsilva/chewBBACA_tutorial) ...
...is available as example on how to run an analysis pipeline using chewBBACA. 

## A ready to use [docker image](https://hub.docker.com/r/mickaelsilva/chewbbaca/) ...
...automatically built from the latest version of chewBBACA in Ubuntu 16.04. 

----------

## 0. Setting up the analysis

**Installing chewBBACA**

Git clone the whole repository

```
git clone https://github.com/mickaelsilva/chewBBACA.git
```

You need to install the following dependencies. Prodigal and BLAST must be added to the PATH variables.

You may use this pip command to install the python dependencies automatically:

```
pip install -r requirements.txt
```


Python dependencies:
* [Biopython 1.70 ](http://biopython.org/wiki/Main_Page)


Other dependencies:
* BLAST 2.5.0+ ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.5.0/ or above
* [Prodigal 2.6.0 ](https://github.com/hyattpd/prodigal/releases/) or above

----------

## 1. wgMLST schema creation

Create your own wgMLST schema based on a set of genomes fasta files. The command is the following:

`chewBBACA.py CreateSchema -i ./genomes/ -o OutputFolderName --cpu 4`

**Parameters**

`-i` Folder containing the genomes from which schema will be created. Alternatively a file 
 containing the path to the list of genomes. One file path (must be full path) 
 to any fasta/multifasta file containing all the complete or draft genomes you want to call alleles for.

`-o` prefix for the output folder for the schema

`--cpu` Number of cpus to use

`-t` (Optional but recommended, contact for new species) taxon name (e.g. Streptococcus agalactiae). It will call the taxon-specific file to be used for training prodigal

`--bsr` (Optional) Minimum BSR for defining locus similarity. Default at 0.6. 

**Outputs:** 

One fasta file per gene in the `-o`directory that is created. 
The fasta file names are the given according the FASTA annotation for each coding sequence. 

----------

## 2.  Allele call using the wgMLST schema 


Then run is the following:

`chewBBACA.py Allelecall -i ./genomes/ -g genes/ -o OutPrefix --cpu 3 `

**Parameters** 

`-i` Folder containing the query genomes. Alternatively a file
 containing the list with the full path of the location of the query genomes.
`-g` Folder containing the reference genes of the schema. Alternatively a file
 containing the list with the full path of the location of the reference genes.  

`-o` prefix for the output directory. ID for the allele call run.

`--cpu` Number of cpus to use 

`-t` (Optional but recommended, contact for new species) taxon name (e.g. Streptococcus agalactiae). It will call the taxon-specific file to be used for training prodigal

`-b` (optional)Blastp full path. In case of slurm system BLAST version being outdated it may 
be hard to use a different one, use this option using the full path of the updated blastp executable



**Outputs files**:
```
./< outPrefix >_< datestamp>/< outPrefix >/results_statistics.txt
./< outPrefix >_< datestamp>/< outPrefix >/results_contigsInfo.txt
./< outPrefix >_< datestamp>/< outPrefix >/results_Alleles.txt 
./< outPrefix >_< datestamp>/< outPrefix >logging_info.txt 
./< outPrefix >_< datestamp>/< outPrefix >RepeatedLoci.txt
```


----------

## 3. Evaluate wgMLST call quality per genome


Usage:


`chewBBACA.py TestGenomeQuality -i alleles.tsv -n 12 -t 200 -s 5 -o OutFolder`  
	
`-i` raw output file from an allele calling (i.e. results_Alleles.txt)

`-n` maximum number of iterations. Each iteration removes a set of genomes over the defined threshold (-t) and recalculates all loci presence percentages.

`-t` maximum threshold, will start at 5. This threshold represents the maximum number of missing loci allowed, for each genome independently, before removing it (genome).

`-s` step to add to each threshold (suggested 5)

`-o` Folder for the analysis files

The output consists in a plot with all thresholds and a removedGenomes.txt file where its 
informed of which genomes are removed per threshold when it reaches a stable point (no more genomes are removed).

Example of an output can be seen [here](http://i.imgur.com/jlTV2vg.png) . This example uses an 
original set of 714 genomes and a scheme consisting of 3266 loci, using a parameter `-n 12`,`-s 5` and `-t 300`.

----------
## 4. Defining the cgMLST schema

 **Creating a clean allelic profile for PHYLOViZ** 
 
Clean a raw output file from an allele calling to a phyloviz readable file.


Basic usage:

`chewBBACA.py ExtractCgMLST -i rawDataToClean.tsv -o output_folders`
	
`-i` raw output file from an allele calling

`-o` output folder (created by the script if not existant yet)

`-r` (optional) list of genes to remove, one per line (e.g. the list of gene detected by ParalogPrunning.py)

`-g` (optional) list of genomes to remove, one per line (e.g. list of genomes to be removed selected based on testGenomeQuality results) 

`-p` (optional) minimum percentage of loci presence (e.g 0.95 to get a matrix with the loci that are present in at least 95% of the genomes)

----------
## 5. All in one option 

### Run all chewBBACA analyses in a single command: 

**if you have no schema**

`cd` to the folder where you want to do the analysis, create folders for:
folder 1 : genomes fasta files as base for schema creation
folder 2 : genomes fasta files to call the alleles (genomes from folder 1 will already be used)

`/home/user/chewBBACA/fullBBACA.py --cs genomes/cg/ 
--cpu 6 -o myAnalysis -t Streptococcus agalactiae --genomes ./genomes/other/`

**if you have a schema**

`cd` to the folder where you want to do the analysis, create folders for:
folder 1 : genomes fasta files to call the alleles
folder 2 : target genes fasta files 


`/home/user/chewBBACA/fullBBACA.py --genes schema_seed/ 
--cpu 6 -o myAnalysis -t Streptococcus agalactiae --genomes ./genomes/other/`

`--cs` path to folder with genomes to create schema

`--genes` path to folder with target genes

`--genomes` path to folder with target genomes

`--cpu` Number of cpus to use (dont use all your cpu if using a laptop)

`-o` folder name with the analysis files

`-t` (Optional) taxon to use for prodigal training input


----------
## FAQ

### Q: Which is the fastest way to run?  
A: Check step 1

### Q: Step 2 is taking hours, will it ever end?  
A: Depending on the variability of the strains used to create the schema and the number 
of CPUs you have selected, the computing time used will vary. The more variable the strains, the more BLAST 
comparisons will be made, meaning more time will be needed for finishing the analysis.

### Q: Step 3 just crashed at 99% after 2 days running, do I need to start over?  
A: chewBBACA should allow you to continue where you stopped, just re-run the same command and you should be prompted to continue the allele call or use the flag --fc.

### Q: I ran all the steps and my cgMLST loci size is smaller than traditional MLST, does this even work?  
A: You probably forgot to eliminate from the analysis genomes responsible for a considerable loss of loci. 
Try to run again step 4, remove some of those genomes and check if the cgMLST loci number rises.

### Q: Which species already have a training file?  
A: At the moment:
 - *Campylobacter jejuni*
 - *Acinetobacter baumannii*
 - *Haemophilus influenzae*
 - *Streptococcus agalactiae*
 - *Yersinia enterocolitica*
 - *Enterococcus faecium*
 - *Staphylococcus haemolyticus*
 - *Salmonella enterica enteritidis*
 - *Staphylococcus aureus*

----------

