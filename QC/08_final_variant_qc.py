import hail as hl
import sys

from ukb_utils import hail_init
from ukb_utils import genotypes

# Inputs
MT = sys.argv[1]
IMPUTESEX_TABLE = sys.argv[2]
SUPERPOPS = sys.argv[3]
SEXCHECK_LIST = sys.argv[4]
RELATED_SAMPLES = sys.argv[5]
INITIAL_VARIANT_LIST = sys.argv[6]
SAMPLE_LIST_INITIAL_QC = sys.argv[7]
REFERENCE = 'GRCh38'

# Outputs
VARIANT_QC_FILE = sys.argv[8]

hail_init.hail_bmrc_init('logs/hail/hail_export.log', REFERENCE)

# # Inputs:
# TRANCHE = '200k'
# CHR = '20'
# MT = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/filtered/ukb_wes_' + TRANCHE + '_filtered_hardcalls_chr' + CHR + '.mt'
# IMPUTESEX_TABLE = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/04_imputesex_'
# SUPERPOPS = "/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/superpopulation_labels.tsv"
# SEXCHECK_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/04_sexcheck.remove.sample_list'
# RELATED_SAMPLES = '/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/07_king.related.sample_list'
# INITIAL_VARIANT_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/variants/02_prefilter_chr' + CHR + '.keep.variant.ht'
# SAMPLE_LIST_INITIAL_QC = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/03_initial_qc.keep.sample_list'

# # Outputs:
# VARIANT_QC_FILE = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/variants/08_final_qc.variants_chr' + CHR + '_'

print("Inputs:")
print('MT; input matrix table: ', MT)
print('IMPUTESEX_TABLE; input file prefix of .tsvs file to plot imputed sex information output from 06_0_impute_sex_superpopulations.py: ', IMPUTESEX_TABLE)
print('SUPERPOPS; list of 1000G labels to loop over:', SUPERPOPS)
print('SEXCHECK_LIST; list of sex swaps to remove:', SEXCHECK_LIST)
print('RELATED_SAMPLES; input set of related samples output from 07_0_ibd.py, 07_0_pc_relate.py, or 07_0_ukb_relatedness_king.py', RELATED_SAMPLES)
print('INITIAL_VARIANT_LIST; set of variants that pass initial filtering \
	(output from 02_prefilter_variants.py): ', INITIAL_VARIANT_LIST)
print('SAMPLE_LIST_INITIAL_QC; set of initial samples output from 03_01_initial_sample_qc_filter.r: ', SAMPLE_LIST_INITIAL_QC)

print("Outputs:")
print('VARIANT_QC_FILE; output .tsv file variant QC information: ', VARIANT_QC_FILE)

ht_superpops = hl.import_table(SUPERPOPS, impute=True).key_by("sample.ID").select("classification_strict")
ht_related_samples = hl.import_table(RELATED_SAMPLES, no_header=True, key='f0')
ht_initial_variants = hl.read_table(INITIAL_VARIANT_LIST)
ht_sexcheck_samples = hl.import_table(SEXCHECK_LIST, no_header=True, key='f0')
ht_initial_samples = hl.import_table(SAMPLE_LIST_INITIAL_QC, no_header=True, key='f0')

mt = hl.read_matrix_table(MT)
mt = mt.annotate_cols(**ht_superpops[mt.s])
mt = mt.filter_rows(hl.is_defined(ht_initial_variants[mt.row_key]))
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_cols(~hl.is_defined(ht_sexcheck_samples[mt.col_key]))
mt = mt.filter_cols(~hl.is_defined(ht_related_samples[mt.col_key]))

for pop in ["AFR", "AMR", "EAS", "EUR", "SAS"]:
	impute_sex_annotations = hl.read_table(IMPUTESEX_TABLE + pop + '.ht')
	mt = mt.annotate_cols(filter_pop = (mt.classification_strict == pop))
	mt_pop = mt.filter_cols(mt.filter_pop == True)
	print(mt_pop.count())
	mt_pop = mt_pop.annotate_cols(imputesex = impute_sex_annotations[mt_pop.col_key])
	mt_pop = hl.variant_qc(mt_pop, name='variant_qc')
	mt_pop = mt_pop.annotate_rows(
		variant_qc=mt_pop.variant_qc.annotate(AC=mt_pop.variant_qc.AC[1],
		AF=mt_pop.variant_qc.AF[1],
		homozygote_count=mt_pop.variant_qc.homozygote_count[1]))
	mt_pop = mt_pop.annotate_rows(
		variant_qc = mt_pop.variant_qc.annotate(p_value_hwe = hl.case()
		.when(mt_pop.locus.in_autosome(), mt_pop.variant_qc.p_value_hwe)
		.default(hl.agg.filter(mt_pop.imputesex.impute_sex.is_female,
			hl.agg.hardy_weinberg_test(mt_pop.GT).p_value)))
	)
	mt_pop = mt_pop.annotate_rows(
		variant_qc = mt_pop.variant_qc.annotate(het_freq_hwe = hl.case()
		.when(mt_pop.locus.in_autosome(), mt_pop.variant_qc.het_freq_hwe)
		.default(hl.agg.filter(mt_pop.imputesex.impute_sex.is_female,
			hl.agg.hardy_weinberg_test(mt_pop.GT).het_freq_hwe)))
	)
	mt_pop.rows().select('variant_qc').flatten().export(f'{VARIANT_QC_FILE}{pop}.tsv.bgz')

