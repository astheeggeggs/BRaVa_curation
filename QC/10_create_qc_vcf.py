import hail as hl
import sys

from ukb_utils import hail_init
from ukb_utils import genotypes

# Inputs
QC_MT_PREFIX = sys.argv[1]

# Outputs:
QC_VCF_PREFIX = sys.argv[2]
REFERENCE = 'GRCh38'

hail_init.hail_bmrc_init_local('logs/hail/hail_export.log', REFERENCE)

# Inputs:
QC_MT_PREFIX = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/final_mt/10_strict_filtered_chr'

# Outputs:
QC_VCF_PREFIX = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/final_mt/10_strict_filtered_chr'

for pop in ["AFR", "AMR", "EAS", "EUR", "SAS"]:

	QC_MT = QC_MT_PREFIX + CHR + '.' + pop + '.mt'
	mt = hl.read_matrix_table(QC_MT)
	QC_VCF = QC_VCF_PREFIX + CHR + '.' + pop + '.vcf.bgz'
	hl.export_vcf(mt, output=QC_VCF)
