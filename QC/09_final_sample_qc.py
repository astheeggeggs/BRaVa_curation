import hail as hl
import sys

from ukb_utils import hail_init
from ukb_utils import genotypes

# # Inputs:
# CHR = '20'
# MT_HARDCALLS = '/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/filtered/ukb_wes_200k_filtered_hardcalls_chr' + CHR + '.mt'
# IMPUTESEX_TABLE = '/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/04_imputesex_'
# SUPERPOPS = "/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/superpopulation_labels.tsv"
# SEXCHECK_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/04_sexcheck.remove.sample_list'
# SAMPLE_LIST_INITIAL_QC = '/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/03_initial_qc.keep.sample_list'
# INITIAL_VARIANT_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/02_prefilter_chr' + CHR + '.keep.variant.ht'
# FINAL_VARIANT_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.pop.keep.variant_list'

# # Outputs:
# SAMPLE_BEFORE_QC_FILE = '/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc_chr' + CHR + '.before'
# SAMPLE_AFTER_QC_FILE = '/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc_chr' + CHR + '.after'

# Inputs
MT_HARDCALLS = sys.argv[1]
IMPUTESEX_TABLE = sys.argv[2]
SUPERPOPS = sys.argv[3]
SEXCHECK_LIST = sys.argv[4]
SAMPLE_LIST_INITIAL_QC = sys.argv[5]
INITIAL_VARIANT_LIST = sys.argv[6]
FINAL_VARIANT_LIST = sys.argv[7]
REFERENCE = 'GRCh38'

# Outputs:
SAMPLE_BEFORE_QC_FILE = sys.argv[8]
SAMPLE_AFTER_QC_FILE = sys.argv[9]

hail_init.hail_bmrc_init('logs/hail/hail_export.log', 'GRCh38')

ht_superpops = hl.import_table(SUPERPOPS, impute=True).key_by("sample.ID").select("classification_strict")
ht_initial_variants = hl.read_table(INITIAL_VARIANT_LIST)
ht_initial_samples = hl.import_table(SAMPLE_LIST_INITIAL_QC, no_header=True, key='f0')
ht_sexcheck_samples = hl.import_table(SEXCHECK_LIST, no_header=True, key='f0')

ht_final_variants = hl.import_table(FINAL_VARIANT_LIST,
	types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)

mt_before = hl.read_matrix_table(MT_HARDCALLS)
mt_before = mt_before.annotate_cols(**ht_superpops[mt_before.s])
mt_before = mt_before.filter_rows(hl.is_defined(ht_initial_variants[mt_before.row_key]))

# Loop over the 1000 genomes labels
for pop in ["AFR", "AMR", "EAS", "EUR", "SAS"]:
	impute_sex_annotations = hl.read_table(IMPUTESEX_TABLE + pop + '.ht')
	mt_before = mt_before.annotate_cols(filter_pop = (mt_before.classification_strict == pop))
	mt_before_pop = mt_before.filter_cols(mt_before.filter_pop == True)
	mt_before_pop = mt_before_pop.filter_cols(hl.is_defined(ht_initial_samples[mt_before_pop.col_key]))
	mt_before_pop = mt_before_pop.filter_cols(~hl.is_defined(ht_sexcheck_samples[mt_before_pop.col_key]))
	mt_before_pop = hl.variant_qc(mt_before_pop, name = 'variant_qc')
	mt_before_pop = mt_before_pop.annotate_rows(
		variant_qc = mt_before_pop.variant_qc.annotate(
			AC=mt_before_pop.variant_qc.AC[1],
			AF = mt_before_pop.variant_qc.AF[1],
			homozygote_count = mt_before_pop.variant_qc.homozygote_count[1]
			)
		)
	mt_before_pop = mt_before_pop.filter_rows(
		(mt_before_pop.variant_qc.AF > 0) & (mt_before_pop.variant_qc.AF < 1)
		)
	mt_before_pop = hl.sample_qc(mt_before_pop, name='sample_qc')
	ht_final_variants_pop = ht_final_variants.annotate(filter_pop = (ht_final_variants.pop == pop))
	ht_final_variants_pop = ht_final_variants_pop.filter(ht_final_variants_pop.filter_pop)
	mt_after_pop = mt_before_pop.filter_rows(hl.is_defined(ht_final_variants[mt_before_pop.row_key]))
	mt_after_pop = hl.sample_qc(mt_after_pop, name='sample_qc')
	mt_before_pop.cols().select("sample_qc").flatten().export(f'{SAMPLE_BEFORE_QC_FILE}.{pop}.samples.tsv.bgz')
	mt_after_pop.cols().select("sample_qc").flatten().export(f'{SAMPLE_AFTER_QC_FILE}.{pop}.samples.tsv.bgz')
