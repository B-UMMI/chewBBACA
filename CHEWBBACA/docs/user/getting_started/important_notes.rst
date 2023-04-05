Important Notes
===============

- chewBBACA only works with **Python 3** (automatic testing for Python 3.8-3.11
  with GitHub Actions).
- We strongly recommend that users install and use **BLAST 2.9.0+** with chewBBACA, as
  chewBBACA's processes have been extensively tested with that version of BLAST.
- chewBBACA includes Prodigal training files for some species. You can consult the list of
  Prodigal training files that are readily available `here <https://github.com/B-UMMI/chewBBACA/tree/master/CHEWBBACA/prodigal_training_files>`_.
  We strongly recommend using the same Prodigal training file for schema creation and allele calling to ensure consistent results.
- chewBBACA defines an allele as a complete Coding DNA Sequence, with start and stop codons
  according to the `NCBI genetic code table 11 <http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>`_
  (identified using `Prodigal <https://github.com/hyattpd/prodigal/releases/>`_ by default, but with the option to provide FASTA
  files with coding sequences). It will automatically exclude any allele for which the DNA sequence does not contain start or stop
  codons and for which the length is not multiple of three. Alleles that contain ambiguous bases are also excluded.
- Make sure that your FASTA files are UNIX format. If they were created in Linux or MacOS
  systems they should be in the correct format, but if they were created in Windows systems,
  you should do a quick conversion using for example `dos2unix <https://waterlan.home.xs4all.nl/dos2unix.html>`_.
- If you are running chewBBACA in an environment with multiple processes accessing the same schema please use the ``--no-inferred`` option (see :doc:`Allele call </user/modules/AlleleCall>`)
