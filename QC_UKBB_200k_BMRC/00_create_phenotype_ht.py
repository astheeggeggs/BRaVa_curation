# Note, this function is not present in the templates - this is because phenotype QC information may not be present for 
# all datasets. This data is used within the 200k UKB WES pipeline to annotate reported sex.

import hail as hl
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes

parser = argparse.ArgumentParser()
parser.add_argument("--tranche", type=str, default='200k')
args = parser.parse_args()

TRANCHE = args.tranche

hail_init.hail_bmrc_init('logs/hail/hail_export.log', 'GRCh38')

# Inputs
PHENOFILE_AUX = '/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/' + TRANCHE + '/QC_phenotypes.tsv.gz'
PHENOTYPES_TABLE = '/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/' + TRANCHE + '/phenotypes.ht'

pheno_table_aux = hl.import_table(PHENOFILE_AUX, key='ID', types={'ID': hl.tstr}, impute=True, missing=['', 'NA'], force_bgz=True)
pheno_table_aux.write(PHENOTYPES_TABLE, overwrite=True)

ht = hl.read_table(PHENOTYPES_TABLE)
n_samples = ht.count()

print(f'n samples: {n_samples}')
pprint(ht.describe())
