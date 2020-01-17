# chewBBACA_NS: Quick Usage

chewBBACA_NS is a special chewBBACA version (see https://github.com/B-UMMI/chewBBACA) that features special functions to communicate with a nomenclature server instance.

It's utilization comprehends two kind of users, authenticated and not authenticated. Authenticated users are allowed to send and store their data (profiles and alleles) on the server, while
non authenticated users are able to fetch all data but not to store. 

<<<<<<< HEAD
Users with sensible data can be matched to non authenticated users.
=======
chewBBACA performs the schema creation and allele calls on complete or draft genomes resulting from de novo assemblers.

The code contained in this repository concerning release v2.0.5. has been published in Microbial Genomics under the title:  
**chewBBACA: A complete suite for gene-by-gene schema creation and strain identification**  - [Link to paper](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000166) 

When using chewBBACA please use the following citation:


Silva M, Machado M, Silva D, Rossi M, Moran-Gilad J, Santos S, Ramirez M, Carriço J. 15/03/2018. M Gen 4(3): [doi:10.1099/mgen.0.000166](doi:10.1099/mgen.0.000166)

# IMPORTANT
As from 25/01/2018 this repo is now fully on **python 3** (tested on >=3.4). The previous python 2 version can be found at https://github.com/B-UMMI/chewBBACA/tree/chewbbaca_py2  
As from 02/02/2018 **chewBBACA is a python package**, install with pip instead of git clone.  
The prodigal training files to be used are now provided with chewBBACA's source code, easing their usage.

# Latest updates

## New in 2.1.0 (05/11/2019)
* New `PrepExternalSchema` implementation: Algorithmic optimizations to improve speed and maintain memory efficiency. New output files with summary information about schema adaptation and excluded sequences and options to control the Blast Score Ratio, minimum sequence length and genetic code values passed to the process.

## chewBBACA released as a galaxy module! 
Many Thanks to Stefano Morabito and Arnold Knijn (https://github.com/aknijn) for EURL VTEC in ISS, Rome ! 
https://toolshed.g2.bx.psu.edu/repository?repository_id=88fd7663075eeae9&changeset_revision=093352878303

## New in 2.0.17 (10/2/2019)
* New alleles also have a timestamp added to the allele name.

## New in 2.0.16 (06/1/2018)
* Corrected bug from 2.0.15 when no prodigal training file provided.

## New in 2.0.15 (06/1/2018)
* Added prodigal training files to the package. They are now at CHEWBBACA/prodigal_training_files.
 
## New in 2.0.13 (18/09/2018)
* when using the function `PrepExternalSchema`, older behavior would remove any locus with a single translation error while the latest change(2.0.12) would not change the original source fasta, this would make the schema unusable. It is now enforced that the alleles that do not translate are removed from the fasta, be sure to backup your data before using this function.
 
## New in 2.0.11 (05/06/2018)
* corrected bug when -h on allele call
* new option for the schema creation. A schema can be created based on a single fasta file, jumping the prodigal gene prediction running. Use `--CDS` and provide a sinfle fasta file on the `-i` input.
 
## New in 2.0.10 (21/05/2018)
* cgMLST profile extraction function (ExtractCgMLST) more efficient (thanks Dillon Barker)
* new option for the allele call, size threshold previously hardcoded at 0.2 can now be changed using the `--st` option. Size threshold is important for the definition of ASM and ALM (alleles smaller/larger than mode).
 
## New in 2.0.9 (04/04/2018)
* blast results during allele call are not saved as a file, instead are piped directly for processing
* new option for the allele call, if genome fasta input is already a fasta of CDS use the `--CDS` option  
 
## New in 2.0.7 (23/02/2018)
* corrected bug that prevented usage of latest blast version (>=2.7.0)
* version flag can now be used `--version`
* instead of calling the main script `chewBBACA.py` you can now use `chewie` (if installed trought pip).
 
## New in 2.0.5 (14/02/2018)
* AlleleCall : -i option accepts a single fasta file now

---------
## Check the [wiki pages](https://github.com/B-UMMI/chewBBACA/wiki) ...
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
>>>>>>> master_cp


<<<<<<< HEAD
----------  
## 0. Setting up software

**Installing chewBBACA**


Install using pip:

```
pip3 install chewbbaca_nserver
```

=======
## A ready to use [docker image](https://hub.docker.com/r/mickaelsilva/chewbbaca_py3/) ...
...automatically built from the latest version of chewBBACA in Ubuntu 16.04. 

## Available [training files](https://github.com/mickaelsilva/prodigal_training_files) ...

...use for a better result, species specific. **Also inside the package now!** e.g. `--ptf  Acinetobacter_baumannii.trn` will now also work!

----------

## 0. Setting up the analysis

**Installing chewBBACA**

Install using pip

```
pip3 install chewbbaca
```

You need to install the following dependencies. Prodigal and BLAST must be added to the PATH variables.

>>>>>>> master_cp

Python dependencies:
* numpy>=1.14.0
* scipy>=0.13.3
* biopython>=1.70
* plotly>=1.12.9
* SPARQLWrapper>=1.8.0
* pandas>=0.22.0
<<<<<<< HEAD
* requests==2.2.1

**Docker image**

```
docker pull mickaelsilva/chewbbaca_nserver
```
=======

Main dependencies:
* BLAST 2.5.0+ ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.5.0/ or above
* [Prodigal 2.6.0 ](https://github.com/hyattpd/prodigal/releases/) or above
>>>>>>> master_cp

Other dependencies (for schema evaluation only):
* [ClustalW2](http://www.clustal.org/download/current/)
* [mafft](https://mafft.cbrc.jp/alignment/software/)


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

`chewie_ns DownloadSchema -i http://chewbbaca.online/app/v1/NS/species/1/schemas/2 -p my_schema/ --cpu 6`

**Parameters**

`-i` url to schema

`-p` folder where to download the schema fastas

`--cpu` Number of cpus to use

<<<<<<< HEAD
`--BSR` (optional) default: 0.6 - BSR value to consider the same allele, see (https://github.com/B-UMMI/chewBBACA/wiki/2.-Allele-Calling for further info)

=======
`--bsr` (Optional) Minimum BSR for defining locus similarity. Default at 0.6. 
>>>>>>> master_cp

`--ptf` (Optional but recommended, contact for new species) path to file of prodigal training file to use.

**Outputs:** 

<<<<<<< HEAD
One fasta file per gene in the `-o`directory that is created. A short folder, the content of the folder is
used on the chewBBACA allele call.
=======
One fasta file per gene in the `-o` directory that is created. 
The fasta file names are the given according the FASTA annotation for each coding sequence. 
>>>>>>> master_cp

**Optional:** 

Information about each locus is almost non existant at this point, the only information directly given by the schema creation is where are located each identified protein on the 
genome (proteinID_Genome.tsv file). A function was added to fetch information on each locus based on the [uniprot SPARQL endpoint](http://sparql.uniprot.org/sparql).

`chewBBACA.py UniprotFinder -i schema_seed/ -t proteinID_Genome.tsv --cpu 4`

**Parameters**

`-i` Folder containing the reference genes of the schema.

`-t` proteinID_Genome.tsv output from the schema creation

`--cpu` Number of cpus to use

**Outputs:** 

A tsv file with the information of each fasta (new_protids.tsv), location on the genome, a name for which the protein sequence was submitted on uniprot and a link to that identified protein. 

----------

## 2.  Allele call using the schema 

The AlleleCall function is to be used on a schema downloaded from the nomenclature server. Local
alleles will be present both on the fastas from the schema and the profile results with an *.

Then run is the following:

`chewie_ns AlleleCall -i ./genomes/ -g genes/ -o OutPrefix --cpu 3 `

**Parameters** 

`-i` Folder containing the query genomes. Alternatively a file
 containing the list with the full path of the location of the query genomes.

`-g` Folder containing the reference genes of the schema. Alternatively a file
 containing the list with the full path of the location of the reference genes.  

`-o` prefix for the output directory. ID for the allele call run.

`--cpu` Number of cpus to use 

<<<<<<< HEAD
`--ptf` (Optional but recommended, contact for new species) provide the prodigal training file (ptf) path. It will call the taxon-specific file to be used for training prodigal

=======
>>>>>>> master_cp
`-b` (optional)Blastp full path. In case of slurm system BLAST version being outdated it may 
be hard to use a different one, use this option using the full path of the updated blastp executable

`--ptf` (Optional but recommended, contact for new species) path to file of prodigal training file to use.


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


`chewie_ns Send2NS -s my_schema/ -p results/results_20171220T154443/results_alleles.tsv`
	
`-s` path to folder where the schema is

`-p` profile result from the allele call (results_alleles.tsv)

`-t` (optional) private authentication token

`--cpu` Number of cpus to use 

`-m` (optional) metadata file


**NOTICE** metadata may be sent later with another function, see 7.

Metadata file should be a tsv file with the following headers:
```
FILE	ST	ACCESSION	COUNTRY	STRAIN	collection_date	host	host_disease	lat	long	isol_source
```

Metadata fields with no info should be completed with NA.

**FILE** column must have the same name as the genome given in `-p`.

**Country** should be one of the present on the DBpedia list:
http://dbpedia.org/class/yago/WikicatMemberStatesOfTheUnitedNations

**STRAIN** is free text .

**ACCESSION** needs to be found either on ENA or SRA.

**collection_date** needs to be on the YYYY-MM-DD format.

**host** needs to be an acceptable taxon name of common name found at https://www.uniprot.org/taxonomy/ . Common names such as 'pig' 'human' will also work.

**host_disease** needs to be a disease ID found at http://www.disease-ontology.org/ . e.g salmonellosis will be 0060859

**lat** needs to be a float number.

**long** needs to be a float number.

**isol_source** is free text.

All metadata fields are optional except FILE.

**Outputs files**:
The output will be a new modified file with the new modified profile (newProfile.tsv). Alleles id with an * are only present locally.

----------
## 4. Sync local schema with NS schema

This function is to be used to get alleles new alleles from the Nomenclature server. This will allow to
update your local schema.

Basic usage:

`chewie_ns SyncSchema -p my_schema/ --cpu 6`
	
`-p` path to schema folder

`--cpu` number of cpu to use

----------

## 5. Download profiles from a given species for a given schema

This function is to be used to get profiles from the Nomenclature server. Given a species id and a
schema id from that species, get a tsv with all the profiles. Due to the way the nomenclature server is built,
getting a profile is costly to the server, for that reason, downloading all the profiles may take a while. It
is advised to specify a list of genomes you are mostly interested in (-r option) or to keep a file where you update
new profiles (-p option)

Basic usage:

`chewie_ns DownloadProfiles --sp http://chewbbaca.online/app/v1/NS/species/2 --sc 1 --cpu 6`
	
`--sp` species uri

`--sc` schema id

`--cpu` number of cpu to use

`-r` (optional) list of specific genomes to download profiles

`-p` (optional) profile with already downloaded profiles for that schema.

`-t` (optional) private authentication token. If given, will download only the isolates that belong to the user.

The output will be a new file with all the profiles (profiles.tsv)


----------

## 6. Manipulating the profile

This function is to be used to clean a raw output 
file from an allele calling to a PhyloViz readable file. Alleles that are only present locally (*alleles)
will be replaced by a high allele number (e.g. 99999)

Basic usage:

`chewie_ns ExtractCgMLST -i rawDataToClean.tsv -o output_folders`
	
`-i` raw output file from an allele calling

`-o` output folder (created by the script if not existent yet)

`-r` (optional) list of genes to remove, one per line (e.g. the list of gene detected by ParalogPrunning.py)

`-g` (optional) list of genomes to remove, one per line (e.g. list of genomes to be removed selected based on testGenomeQuality results) 

`-p` (optional) minimum percentage of loci presence (e.g 0.95 to get a matrix with the loci that are present in at least 95% of the genomes)

----------

## 7. Send metadata to isolates already on the Nomenclature server

This function is to be used to send metadata to isolates that are already on the Nomenclature server.
Metadata that were already uploaded to the server will not be updated.


Usage:


`chewie_ns SendMetadata -t qweqweasd -m info.tsv --cpu 3`

`-t` private authentication token

`--cpu` Number of cpus to use

<<<<<<< HEAD
`-m` metadata file

Metadata file should be a tsv file with the following headers:
```
FILE	ST	ACCESSION	COUNTRY	STRAIN	collection_date	host	host_disease	lat	long	isol_source
```

Metadata fields with no info should be completed with NA.

**FILE** column must have the **isolate URI**.

**Country** should be one of the present on the DBpedia list:
http://dbpedia.org/class/yago/WikicatMemberStatesOfTheUnitedNations

**STRAIN** is free text .

**ACCESSION** needs to be found either on ENA or SRA.

**collection_date** needs to be on the YYYY-MM-DD format.

**host** needs to be an acceptable taxon name of common name found at https://www.uniprot.org/taxonomy/ . Common names such as 'pig' 'human' will also work.

**host_disease** needs to be a disease ID found at http://www.disease-ontology.org/ . e.g salmonellosis will be 0060859

**lat** needs to be a float number.

**long** needs to be a float number.

**isol_source** is free text.

All metadata fields are optional except FILE.

See the info.tsv file for an example.
=======
----------
## FAQ

### Q: Step 2 is taking hours, will it ever end?  
A: Depending on the variability of the strains used to create the schema and the number 
of CPUs you have selected, the computing time used will vary. The more variable the strains, the more BLAST 
comparisons will be made, meaning more time will be needed for finishing the analysis.

### Q: Step 3 just crashed at 99% after 2 days running, do I need to start over?  
A: chewBBACA should allow you to continue where you stopped, just re-run the same command and you should be prompted to continue the allele call or use the flag --fc.

### Q: I ran all the steps and my cgMLST loci size is smaller than traditional MLST, does this even work?  
A: You probably forgot to eliminate from the analysis genomes responsible for a considerable loss of loci. 
Try to run again step 4, remove some of those genomes and check if the cgMLST loci number rises.

### Q: Can I use a schema from an external source?
A: Yes. Be sure to have a single fasta for each locus and use the "PrepExternalSchema​" function.

### Q: Which species already have a training file?  
A: At the moment:
 - *Acinetobacter baumannii*
 - *Campylobacter jejuni*
 - *Enterococcus faecium*
 - *Escherichia coli*
 - *Haemophilus influenzae*
 - *Legionella pneumophila*
 - *Listeria monocytogenes*
 - *Salmonella enterica enteritidis*
 - *Streptococcus agalactiae*
 - *Staphylococcus aureus*
 - *Staphylococcus haemolyticus*
 - *Yersinia enterocolitica*
 
get them at https://github.com/mickaelsilva/prodigal_training_files
 
### Q: My favorite species has no training file. What can I do?
A: You can propose a new one to be added to the repository or create your own training 
files. To create a training file do:

`prodigal -i myGoldStandardGenome.fna -t myTrainedFile.trn -p single`
 

----------  
  
 
## Citation
Silva M, Machado M, Silva D, Rossi M, Moran-Gilad J, Santos S, Ramirez M, Carriço J. 15/03/2018. M Gen 4(3): [doi:10.1099/mgen.0.000166](doi:10.1099/mgen.0.000166)
>>>>>>> master_cp
