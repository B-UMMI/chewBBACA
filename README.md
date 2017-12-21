# chewBBACA_NS: Quick Usage

chewBBACA teams with nomenclature server

example: http://137.205.69.51/app/v1/NS/species/2/schemas/1


----------

## 0. Setting up the analysis

**Installing chewBBACA**

### chewBBACA_ns works with python 3

Git clone the whole repository (advised to clone to a new folder as it uses a different branch)

```
git clone -b chewie_NS https://github.com/B-UMMI/chewBBACA.git
```

You need to install the following dependencies. Prodigal and BLAST must be added to the PATH variables.

You may use this pip command to install the python dependencies automatically:

```
pip3 install -r requirements.txt
```


Python dependencies:
* [Biopython 1.70 ](http://biopython.org/wiki/Main_Page)


Other dependencies:
* BLAST 2.5.0+ ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.5.0/ or above
* [Prodigal 2.6.0 ](https://github.com/hyattpd/prodigal/releases/) or above

----------

## 1. get a schema from NS

Download a schema from NS. The command is the following:

`chewBBACA.py DownloadSchema -i http://137.205.69.51/app/v1/NS/species/2/schemas/1 -p my_schema/ --cpu 6`

**Parameters**

`-i` url to schema

`-o` folder where to download the schema fastas

`--cpu` Number of cpus to use


**Outputs:** 

One fasta file per gene in the `-o`directory that is created. 

----------

## 2.  Allele call using the wgMLST schema 


Then run is the following:

`chewBBACA.py AlleleCall -i ./genomes/ -g genes/ -o OutPrefix --cpu 3 `

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

## 3. Check a local profile against the NS

Usage:


`chewBBACA.py Send2NS -s my_schema/ -p results/results_20171220T154443/results_alleles.tsv`  
	
`-s` folder where the schema is

`-p` profile result from the allele call (results_alleles.tsv)


The output will be a new modified file with the new modified profile (newProfile.tsv). Alleles id with an * are only present locally.

----------
## 4. Sync local schema with NS schema

Get new alleles present on the NS and update your local schema


Basic usage:

`chewBBACA.py SyncSchema -p my_schema/ --cpu 6`
	
`-p` schema folder

`--cpu` number of cpu to use

----------

## 4. Download all profiles from a given species for a given locus

Get all profiles on a tsv.

Basic usage:

`chewBBACA.py DownloadProfiles --sp http://137.205.69.51/app/v1/NS/species/2 --sc http://137.205.69.51/app/v1/NS/species/2/schemas/1 --cpu 6`
	
`--sp` species uri

`--sc` schema uri

`--cpu` number of cpu to use

The output will be a new file with all the profiles (profiles.tsv)


----------
