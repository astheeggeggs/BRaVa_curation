import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes

parser = argparse.ArgumentParser()
parser.add_argument("--chr", type=str, default='20')
parser.add_argument("--tranche", type=str, default='200k')
args = parser.parse_args()

TRANCHE = args.tranche
CHR = str(args.chr)

# Inputs:
QC_MT_PREFIX = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/final_mt/10_strict_filtered_chr'
REFERENCE = 'GRCh38'

# Outputs:
QC_VCF_PREFIX = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/final_mt/10_strict_filtered_chr'

hail_init.hail_bmrc_init('logs/hail/hail_export.log', REFERENCE)

for pop in ["AFR", "AMR", "EAS", "EUR", "SAS"]:
	QC_MT = QC_MT_PREFIX + '.' + pop + '.mt'
	mt = hl.read_matrix_table(QC_MT)
	QC_VCF = QC_VCF_PREFIX + '.' + pop + '.vcf.bgz'
	hl.export_vcf(mt, output=QC_VCF)
