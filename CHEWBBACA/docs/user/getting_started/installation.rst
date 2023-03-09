Installation
============

Conda
.....

Install the latest released version using `conda <https://anaconda.org/bioconda/chewbbaca>`_:

::

	conda create -c bioconda -c conda-forge -n chewie "chewbbaca=3.1.1"

If you're having issues installing chewBBACA through conda, we recommend that you install
`mamba <https://mamba.readthedocs.io/en/latest/index.html>`_ and run the following command:

::

	mamba create -c bioconda -c conda-forge -n chewie "chewbbaca=3.1.1"

Pip
...

Install using `pip <https://pypi.org/project/chewBBACA/>`_:

::

	pip3 install chewbbaca


Python dependencies
...................

* numpy >=1.23.4
* scipy >=1.9.3
* biopython >=1.78
* plotly >=5.8.0
* SPARQLWrapper >=2.0.0
* requests >=2.27.1
* pandas >=1.5.1

.. note::
	These dependencies are defined in the `requirements <https://github.com/B-UMMI/chewBBACA/blob/master/CHEWBBACA/requirements.txt>`_
	file and should be automatically installed when using conda or pip.

Other dependencies
..................

* `BLAST >=2.9.0 <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/>`_
* `Prodigal >=2.6.3 <https://github.com/hyattpd/prodigal/releases/>`_
* `MAFFT >=7.505 <https://mafft.cbrc.jp/alignment/software/>`_ (for schema evaluation only)

.. important::
	Installation through conda should take care of all dependencies. If you install through
	pip you will need to ensure that you have BLAST, Prodigal and MAFFT installed and added to
	the PATH.
