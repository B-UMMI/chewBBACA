from setuptools import setup

import CHEWBBACA

VERSION = CHEWBBACA.__version__

setup(
  name = 'chewBBACA',
  packages = ['CHEWBBACA','CHEWBBACA.allelecall','CHEWBBACA.createschema','CHEWBBACA.utils','CHEWBBACA.SchemaEvaluator'],
  version = VERSION,
  description = 'A complete suite for gene-by-gene schema creation and strain identification',
  author = 'Mickael Silva',
  author_email = 'mickaelsilva@medicina.ulisboa.pt',
  url = 'https://github.com/B-UMMI/chewBBACA', 
  keywords = ['cgMLST', 'bacterial typing'],
  install_requires=['numpy>=1.13.1','scipy>=0.13.3','biopython>=1.70','plotly>=1.12.9','SPARQLWrapper>=1.8.0'],
  entry_points={
        "console_scripts": [
            "chewBBACA.py = CHEWBBACA.chewBBACA:main"
        ]
}
)
