from setuptools import setup

<<<<<<< HEAD
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
=======
import CHEWBBACA

VERSION = CHEWBBACA.__version__

setup(
  name = 'chewBBACA',
  packages = ['CHEWBBACA','CHEWBBACA.allelecall','CHEWBBACA.createschema','CHEWBBACA.utils','CHEWBBACA.SchemaEvaluator','CHEWBBACA.PrepExternalSchema'],
  version = VERSION,
  description = 'A complete suite for gene-by-gene schema creation and strain identification',
  author = 'Mickael Silva',
  author_email = 'mickaelsilva@medicina.ulisboa.pt',
  url = 'https://github.com/B-UMMI/chewBBACA', 
  keywords = ['cgMLST', 'bacterial typing'],
  install_requires=['numpy>=1.14.0','scipy>=0.13.3','biopython>=1.70','plotly>=1.12.9','SPARQLWrapper>=1.8.0','pandas>=0.22.0','requests>=2.2.1'],
  python_requires='>=3.4',
  include_package_data = True,
  entry_points={
        "console_scripts": [
            "chewBBACA.py = CHEWBBACA.chewBBACA:main",
            "chewie = CHEWBBACA.chewBBACA:main"
>>>>>>> master_cp
        ]
}
)
