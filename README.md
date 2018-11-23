# chewBBACA_NS: Quick Usage

chewBBACA_NS is a special chewBBACA version (see https://github.com/B-UMMI/chewBBACA) that features special functions to communicate with a nomenclature server instance.

It's utilization comprehends two kind of users, authenticated and not authenticated. Authenticated users are allowed to send and store their data (profiles and alleles) on the server, while
non authenticated users are able to fetch all data but not to store. 

Users with sensible data can be matched to non authenticated users.


----------  
## 0. Setting up software

**Installing chewBBACA**


Install using pip:

```
pip3 install chewbbaca_ns
```


Python dependencies:
* numpy>=1.14.0
* scipy>=0.13.3
* biopython>=1.70
* plotly>=1.12.9
* SPARQLWrapper>=1.8.0
* pandas>=0.22.0
* requests==2.2.1

**Docker image**

```
docker pull mickaelsilva/chewbbaca_nserver
```

----------

## 1. Download a schema from NS

Download a schema from NS with the DownloadSchema. Schemas are associated to species, each species has its own schemas.
To check the schemas of species 1 you go to http://chewbbaca.online/app/v1/NS/species/1/schemas. The URI that defines the schema should be used
for the -i input.

**Note:**
Downloading a schema with thousands of loci and thousands of alleles from the server is fast, however
due to the need of some calculations for chewBBACA software usage (BSR scores), it may take a while to complete the process.
In order to mitigate this issue we have compressed versions of the schemas available that, even if outdated, are easily
updatable locally through the SyncSchema function (see step 4). The compressed versions can be downloaded from a
browser e.g. http://chewbbaca.online/app/v1/NS/species/1/schemas/2/compressed

The command is the following:

`chewBBACA.py DownloadSchema -i http://chewbbaca.online/app/v1/NS/species/1/schemas/2 -p my_schema/ --cpu 6`

**Parameters**

`-i` url to schema

`-p` folder where to download the schema fastas

`--cpu` Number of cpus to use

`--BSR` (optional) default: 0.6 - BSR value to consider the same allele, see (https://github.com/B-UMMI/chewBBACA/wiki/2.-Allele-Calling for further info)


**Outputs:** 

One fasta file per gene in the `-o`directory that is created. A short folder, the content of the folder is
used on the chewBBACA allele call.

----------

## 2.  Allele call using the schema 

The AlleleCall function is to be used on a schema downloaded from the nomenclature server. Local
alleles will be present both on the fastas from the schema and the profile results with an *.

Then run is the following:

`chewBBACA.py AlleleCall -i ./genomes/ -g genes/ -o OutPrefix --cpu 3 `

**Parameters** 

`-i` Folder containing the query genomes. Alternatively a file
 containing the list with the full path of the location of the query genomes.

`-g` Folder containing the reference genes of the schema. Alternatively a file
 containing the list with the full path of the location of the reference genes.  

`-o` prefix for the output directory. ID for the allele call run.

`--cpu` Number of cpus to use 

`--ptf` (Optional but recommended, contact for new species) provide the prodigal training file (ptf) path. It will call the taxon-specific file to be used for training prodigal

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

## 3. Push a local profile and it's respective alleles to the NS

This function is to be used to send alleles and profiles to the Nomenclature server. Only the alleles
with a * present in the profile will be uploaded to the server. The new alleles and profiles will only
be submitted if the authentication token is provided, else they are queried against the server alleles database and will return
the allele server number if it has already been attributed


Usage:


`chewBBACA.py Send2NS -s my_schema/ -p results/results_20171220T154443/results_alleles.tsv`  
	
`-s` path to folder where the schema is

`-p` profile result from the allele call (results_alleles.tsv)

`-t` (optional) private authentication token

`--cpu` Number of cpus to use 

`-m` (optional) metadata file


Metadata file should be a tsv file with the following headers:
```FILE ACCESSION COUNTRY ST STRAIN```

Country should be one of the present on the DBpedia list:
http://dbpedia.org/class/yago/WikicatMemberStatesOfTheUnitedNations


**Outputs files**:
The output will be a new modified file with the new modified profile (newProfile.tsv). Alleles id with an * are only present locally.

----------
## 4. Sync local schema with NS schema

This function is to be used to get alleles new alleles from the Nomenclature server. This will allow to
update your local schema.

Basic usage:

`chewBBACA.py SyncSchema -p my_schema/ --cpu 6`
	
`-p` path to schema folder

`--cpu` number of cpu to use

----------

## 4. Download profiles from a given species for a given schema

This function is to be used to get profiles from the Nomenclature server. Given a species id and a
schema id from that species, get a tsv with all the profiles. Due to the way the nomenclature server is built,
getting a profile is costly to the server, for that reason, downloading all the profiles may take a while. It
is advised to specify a list of genomes you are mostly interested in (-r option) or to keep a file where you update
new profiles (-p option)

Basic usage:

`chewBBACA.py DownloadProfiles --sp http://chewbbaca.online/app/v1/NS/species/2 --sc 1 --cpu 6`
	
`--sp` species uri

`--sc` schema id

`--cpu` number of cpu to use

`-r` (optional) list of specific genomes to download profiles

`-p` (optional) profile with already downloaded profiles for that schema.

The output will be a new file with all the profiles (profiles.tsv)


----------

## 5. Defining the cgMLST profile

This function is to be used to clean a raw output 
file from an allele calling to a PhyloViz readable file. Alleles that are only present locally (*alleles)
will be replaced by a high allele number (e.g. 99999)

Basic usage:

`chewBBACA.py ExtractCgMLST -i rawDataToClean.tsv -o output_folders`
	
`-i` raw output file from an allele calling

`-o` output folder (created by the script if not existent yet)

`-r` (optional) list of genes to remove, one per line (e.g. the list of gene detected by ParalogPrunning.py)

`-g` (optional) list of genomes to remove, one per line (e.g. list of genomes to be removed selected based on testGenomeQuality results) 

`-p` (optional) minimum percentage of loci presence (e.g 0.95 to get a matrix with the loci that are present in at least 95% of the genomes)
