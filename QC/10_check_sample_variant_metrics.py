import hail as hl
import sys

from ukb_utils import hail_init
from ukb_utils import genotypes

hail_init.hail_bmrc_init('logs/hail/hail_export.log', 'GRCh38')

# TRANCHE="200k"
# CURATED_MT="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/final_mt/10_european.strict_filtered_chr20.mt"
# SAMPLE_QC_FILE="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/10_sample_metrics_for_plotting_chr20.tsv.gz"
# SAMPLE_QC_TARGET_FILE="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/10_sample_target_interval_metrics_for_plotting_chr20.tsv.gz"
# VARIANT_QC_FILE="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/10_variant_metrics_for_plotting_chr20.tsv.gz"
# IMPUTESEX_TABLE="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/04_imputesex.ht"

# Inputs
CURATED_MT = sys.argv[1]

# Outputs
VARIANT_QC_FILE = sys.argv[2]
SAMPLE_QC_FILE = sys.argv[3]
SAMPLE_QC_TARGET_FILE = sys.argv[4]
IMPUTESEX_TABLE = sys.argv[5]

print("Inputs:")
print('CURATED_MT; input matrix table: ', CURATED_MT)

print("Outputs:")
print('VARIANT_QC_FILE; output .tsv file variant QC information: ', VARIANT_QC_FILE)
print('SAMPLE_QC_FILE; output .tsv file variant QC information: ', SAMPLE_QC_FILE)
print('SAMPLE_QC_TARGET_FILE; output .tsv file variant QC information: ', SAMPLE_QC_TARGET_FILE)
print('IMPUTESEX_TABLE; output hail table with imputed sex: ', IMPUTESEX_TABLE)

mt = hl.read_matrix_table(CURATED_MT)
# Filter down to the target intervals
target_path = '/well/lindgren-ukbb/projects/ukbb-11867/nbaya/resources/ref/xgen_plus_spikein.b38.chr_prefix.bed'
target = hl.import_bed(target_path, reference_genome='GRCh38')
mt = mt.annotate_rows(in_target = hl.is_defined(target[mt.locus]))

mt = hl.variant_qc(mt, name='variant_qc')

mt = mt.annotate_rows(
	qc=mt.variant_qc.annotate(AC=mt.variant_qc.AC[1],
	AF=mt.variant_qc.AF[1],
	homozygote_count=mt.variant_qc.homozygote_count[1]))
mt = mt.filter_rows((mt.qc.AF > 0) & (mt.qc.AF < 1))

ht = hl.read_table(IMPUTESEX_TABLE)
mt = mt.annotate_cols(**ht[mt.col_key])

mt = mt.annotate_rows(variant_qc = mt.variant_qc.annotate(p_value_hwe = hl.case()
	.when(mt.locus.in_autosome(), mt.variant_qc.p_value_hwe)
	.default(hl.agg.filter(mt.impute_sex.is_female,
		hl.agg.hardy_weinberg_test(mt.GT).p_value)))
)

mt = hl.sample_qc(mt, name='sample_qc')

# Spit out the resultant files
ht = mt.cols().select('sample_qc').flatten()
ht.export(SAMPLE_QC_FILE)

ht = mt.rows().select('variant_qc', 'in_target').flatten()
ht.export(VARIANT_QC_FILE)

mt = mt.filter_rows(mt.in_target)
mt = hl.sample_qc(mt, name='sample_qc')
ht = mt.cols().select('sample_qc').flatten()
ht.export(SAMPLE_QC_TARGET_FILE)

