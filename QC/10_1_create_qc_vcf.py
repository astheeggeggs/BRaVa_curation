import hail as hl
import sys

from ukb_utils import hail_init
from ukb_utils import genotypes

# Inputs
QC_MT_PREFIX = sys.argv[1]
REFERENCE = 'GRCh38'

# Outputs:
QC_VCF_PREFIX = sys.argv[2]

hl.init(default_reference=REFERENCE)

print("Inputs:")
print('QC_MT_PREFIX : prefix for 1000 genomes label specific matrix tables', QC_MT_PREFIX)

print("Outputs:")
print('QC_VCF_PREFIX : prefix for 1000 genomes label specific vcf.bgz files', QC_VCF_PREFIX)

for pop in ["AFR", "AMR", "EAS", "EUR", "SAS"]:
	QC_MT = QC_MT_PREFIX + '.' + pop + '.mt'
	mt = hl.read_matrix_table(QC_MT)
	QC_VCF = QC_VCF_PREFIX + '.' + pop + '.vcf.bgz'
	hl.export_vcf(mt, output=QC_VCF)
