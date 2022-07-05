import hail as hl
import sys

from ukb_utils import hail_init
from ukb_utils import genotypes

# Inputs
MT = sys.argv[1]
FINAL_SAMPLE_LIST  = sys.argv[2]
FINAL_VARIANT_LIST  = sys.argv[3]
SUPERPOPS = sys.argv[4]
REFERENCE = 'GRCh38'

# Outputs:
QC_MT_PREFIX = sys.argv[5]

print("Inputs:")
print('MT; input matrix table: ', MT_HARDCALLS)
print('FINAL_SAMPLE_LIST; output .sample_list file from 09_filter_final_sample_qc.r: ', FINAL_SAMPLE_LIST)
print('FINAL_VARIANT_LIST; output .variant_list file from 08_final_variant_qc_filter.r: ', FINAL_VARIANT_LIST)
print('SUPERPOPS; list of 1000G labels to loop over:', SUPERPOPS)

print("Outputs:")
print('QC_MT_PREFIX : prefix for 1000 genomes label specific matrix tables', QC_MT_PREFIX)

hl.init(default_reference=REFERENCE)

# Outputs
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
