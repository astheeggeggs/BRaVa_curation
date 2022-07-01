import hail as hl
import sys

from ukb_utils import hail_init
from ukb_utils import genotypes

# Inputs
MT = sys.argv[1]
FINAL_SAMPLE_LIST  = sys.argv[2]
FINAL_VARIANT_LIST  = sys.argv[3]
SUPERPOPS = sys.argv[4]

# Outputs:
QC_MT_PREFIX = sys.argv[5]
REFERENCE = 'GRCh38'

hail_init.hail_bmrc_init('logs/hail/hail_export.log', REFERENCE)

# # Inputs:
# MT  = '/well/lindgren/UKBIOBANK/nbaya/wes_' + TRANCHE + '/ukb_wes_qc/data/filtered/ukb_wes_' + TRANCHE + '_filtered_chr' + CHR + '.mt'
# FINAL_SAMPLE_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/09_final_qc.keep.BRaVa.sample_list'
# FINAL_VARIANT_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/variants/08_final_qc.pop.keep.variant_list'
# SUPERPOPS = '/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/superpopulation_labels.tsv'

# # Outputs
# QC_MT_PREFIX = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/final_mt/10_strict_filtered_chr'

ht_final_samples = hl.import_table(FINAL_SAMPLE_LIST, no_header=True, key='f0', delimiter=',', types={'f0': hl.tstr})
ht_final_variants = hl.import_table(FINAL_VARIANT_LIST, types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)
ht_superpops = hl.import_table(SUPERPOPS, impute=True, types={'f0': hl.tstr}).key_by("sample.ID").select("classification_strict")

mt = hl.read_matrix_table(MT)
mt = mt.annotate_cols(pops = ht_superpops[mt.s])
mt = mt.drop('qual', 'info', 'filters')

for pop in ["AFR", "AMR", "EAS", "EUR", "SAS"]:
	mt_tmp = mt.filter_cols((mt.pops.classification_strict == pop))
	# Filter
	mt_tmp = mt_tmp.filter_cols(hl.is_defined(ht_final_samples[mt_tmp.col_key]))
	ht_final_variants_tmp = ht_final_variants.filter(ht_final_variants.pop==pop)
	mt_tmp = mt_tmp.filter_rows(hl.is_defined(ht_final_variants_tmp[mt_tmp.row_key]))
	# Write the result
	QC_MT = QC_MT_PREFIX + '.' + pop + '.mt'
	mt_tmp.write(QC_MT, overwrite=True)
