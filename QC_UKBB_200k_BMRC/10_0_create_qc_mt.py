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

# Inputs
MT  = '/well/lindgren/UKBIOBANK/nbaya/wes_200k/ukb_wes_qc/data/filtered/ukb_wes' + TRANCHE + '_filtered_chr' + str(CHR) + '.mt'
FINAL_SAMPLE_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes' + TRANCHE + '/ukb_wes_qc/data/samples/09_final_qc.keep.BRaVa.sample_list'
FINAL_VARIANT_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes' + TRANCHE + '/ukb_wes_qc/data/variants/08_final_qc.pop.keep.variant_list'
SUPERPOPS = '/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/superpopulation_labels.tsv'
REFERENCE = 'GRCh38'

# Outputs
QC_MT_PREFIX = '/well/lindgren/UKBIOBANK/dpalmer/wes' + TRANCHE + '/ukb_wes_qc/data/final_mt/10_strict_filtered_chr' + str(CHR)

hail_init.hail_bmrc_init('logs/hail/hail_export.log', REFERENCE)

ht_final_samples = hl.import_table(FINAL_SAMPLE_LIST, no_header=True, key='f0', delimiter=',', types={'f0': hl.tstr})
ht_final_variants = hl.import_table(FINAL_VARIANT_LIST, types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)
ht_superpops = hl.import_table(SUPERPOPS, impute=True, types={'f0': hl.tstr}).key_by("sample.ID").select("classification_strict")

# Loop over populations
for pop in ["AFR", "AMR", "EAS", "EUR", "SAS"]:
	mt = hl.read_matrix_table(MT)
	mt = mt.annotate_cols(pops = ht_superpops[mt.s])
	mt = mt.drop('qual', 'info', 'filters')
	mt = mt.filter_cols((mt.pops.classification_strict == pop))
	# Filter
	mt = mt.filter_cols(hl.is_defined(ht_final_samples[mt.col_key]))
	ht_final_variants_tmp = ht_final_variants.filter(ht_final_variants.pop==pop)
	mt = mt.filter_rows(hl.is_defined(ht_final_variants_tmp[mt.row_key]))
	# Write the result
	QC_MT = QC_MT_PREFIX + '.' + pop + '.mt'
	mt.write(QC_MT, overwrite=True)
