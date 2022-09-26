Installation
============

Install the latest released version using `conda <https://anaconda.org/bioconda/chewbbaca>`_:

::

	conda create -c bioconda -c conda-forge -n chewie chewbbaca=2.8.5 blast=2.9

Install using `pip <https://pypi.org/project/chewBBACA/>`_:

::

	pip3 install chewbbaca


Python dependencies (defined in the `requirements <https://github.com/B-UMMI/chewBBACA/blob/master/CHEWBBACA/requirements.txt>`_ file, should be automatically installed when using conda or pip):

* numpy>=1.14.0
* scipy>=0.13.3
* biopython>=1.70
* plotly>=1.12.9
* SPARQLWrapper>=1.8.0
* requests>=2.2.1
* pandas>=0.22.0

Main dependencies: 

* `BLAST 2.9.0+ <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/>`_
* `Prodigal 2.6.0 <https://github.com/hyattpd/prodigal/releases/>`_
* `MAFFT <https://mafft.cbrc.jp/alignment/software/>`_ (for schema evaluation only)

Installation through conda should take care of all dependencies. If you install through pip you will need to ensure that you have BLAST, Prodigal and MAFFT installed and added to the PATH.
