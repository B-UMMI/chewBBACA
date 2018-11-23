from setuptools import setup

import CHEWBBACA_NS

VERSION = CHEWBBACA_NS.__version__

setup(
  name = 'chewBBACA_NServer',
  packages = ['CHEWBBACA_NS','CHEWBBACA_NS.allelecall','CHEWBBACA_NS.utils','CHEWBBACA_NS.SchemaEvaluator'],
  version = VERSION,
  description = 'chewBBACA comunicates with the nomenclature server HUUURRR',
  author = 'Mickael Silva',
  author_email = 'mickaelsilva@medicina.ulisboa.pt',
  url = 'https://github.com/B-UMMI/chewBBACA/tree/chewie_NS', 
  keywords = ['cgMLST', 'bacterial typing','nomenclature server'],
  install_requires=['numpy>=1.14.0','scipy>=0.13.3','biopython>=1.70','plotly>=1.12.9','SPARQLWrapper>=1.8.0','requests==2.2.1','pandas>=0.22.0'],
  entry_points={
        "console_scripts": [
            "chewie_ns = CHEWBBACA_NS.chewBBACA:main"
        ]
}
)
