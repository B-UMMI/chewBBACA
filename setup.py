from setuptools import setup

import CHEWBBACA

VERSION = CHEWBBACA.__version__

setup(
  name = 'chewBBACA',
  packages = ['CHEWBBACA', 'CHEWBBACA.AlleleCall', 'CHEWBBACA.CreateSchema',
              'CHEWBBACA.utils', 'CHEWBBACA.SchemaEvaluator', 'CHEWBBACA.PrepExternalSchema',
              'CHEWBBACA.CHEWBBACA_NS', 'CHEWBBACA.UniprotFinder'],
  version = VERSION,
  description = 'A complete suite for gene-by-gene schema creation and strain identification.',
  author = 'Mickael Silva, Pedro Cerqueira, Rafael Mamede',
  author_email = 'imm-bioinfo@medicina.ulisboa.pt',
  url = 'https://github.com/B-UMMI/chewBBACA', 
  keywords = ['cgMLST', 'bacterial typing', 'nomenclature server'],
  install_requires = ['numpy>=1.23.4', 'scipy>=1.9.3', 'biopython>=1.78',
                      'plotly>=5.8.0', 'SPARQLWrapper>=2.0.0', 'pandas>=1.5.1',
                      'requests>=2.27.1'],
  python_requires = '>=3.7',
  include_package_data = True,
  entry_points={'console_scripts': ["chewBBACA.py = CHEWBBACA.chewBBACA:main",
                                    "chewie = CHEWBBACA.chewBBACA:main"]
                }
)
