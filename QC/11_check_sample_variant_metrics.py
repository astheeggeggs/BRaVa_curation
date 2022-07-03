import hail as hl
import sys

from ukb_utils import hail_init
from ukb_utils import genotypes

hail_init.hail_bmrc_init('logs/hail/hail_export.log', 'GRCh38')

# Inputs
QC_MT_PREFIX = sys.argv[1]

# Outputs
VARIANT_QC_FILE = sys.argv[2]
SAMPLE_QC_FILE = sys.argv[3]
SAMPLE_QC_TARGET_FILE = sys.argv[4]
IMPUTESEX_TABLE = sys.argv[5]

TRANCHE="200k"
CHR=22
QC_MT_PREFIX= "/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/final_mt/10_strict_filtered_chr" + str(CHR)
VARIANT_QC_FILE="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/10_variant_metrics_for_plotting_chr" + str(CHR)
SAMPLE_QC_FILE="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/10_sample_metrics_for_plotting_chr" + str(CHR)
SAMPLE_QC_TARGET_FILE="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/10_sample_target_interval_metrics_for_plotting_chr" + str(CHR)
IMPUTESEX_TABLE="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/04_imputesex_"

print("Inputs:")
print('CURATED_MT; input matrix table: ', CURATED_MT)

print("Outputs:")
print('VARIANT_QC_FILE; output .tsv file variant QC information: ', VARIANT_QC_FILE)
print('SAMPLE_QC_FILE; output .tsv file variant QC information: ', SAMPLE_QC_FILE)
print('SAMPLE_QC_TARGET_FILE; output .tsv file variant QC information: ', SAMPLE_QC_TARGET_FILE)
print('IMPUTESEX_TABLE; output hail table with imputed sex: ', IMPUTESEX_TABLE)

# Grab target intervals
target_path = '/well/lindgren-ukbb/projects/ukbb-11867/nbaya/resources/ref/xgen_plus_spikein.b38.chr_prefix.bed'
target = hl.import_bed(target_path, reference_genome='GRCh38')

# Loop over the 1000 genomes labels
for pop in ["AFR", "AMR", "EAS", "SAS"]:#"EUR", "SAS"]:
	# Read in and annotate
	# QC_MT from the previous script worked.
	QC_MT = QC_MT_PREFIX + '.' + pop + '.mt'
	mt = hl.read_matrix_table(QC_MT)
	mt = mt.annotate_rows(in_target = hl.is_defined(target[mt.locus]))
	mt = hl.variant_qc(mt, name='variant_qc')
	mt = mt.annotate_rows(
		qc=mt.variant_qc.annotate(AC=mt.variant_qc.AC[1],
		AF=mt.variant_qc.AF[1],
		homozygote_count=mt.variant_qc.homozygote_count[1]))
	mt = mt.filter_rows((mt.qc.AF > 0) & (mt.qc.AF < 1))
	# Annotate the sexes
	IMPUTESEX_TABLE_tmp = IMPUTESEX_TABLE + pop + '.ht'
	ht = hl.read_table(IMPUTESEX_TABLE_tmp)
	mt = mt.annotate_cols(**ht[mt.col_key])
	# Annotate the variant information
	mt = mt.annotate_rows(variant_qc = mt.variant_qc.annotate(p_value_hwe = hl.case()
		.when(mt.locus.in_autosome(), mt.variant_qc.p_value_hwe)
		.default(hl.agg.filter(mt.impute_sex.is_female,
			hl.agg.hardy_weinberg_test(mt.GT).p_value)))
	)
	# Annotate the sample information
	mt = hl.sample_qc(mt, name='sample_qc')
	# Spit out the resultant files
	ht = mt.cols().select('sample_qc').flatten()
	ht.export(f'{SAMPLE_QC_FILE}.{pop}.tsv.bgz')
	ht = mt.rows().select('variant_qc', 'in_target').flatten()
	ht.export(f'{VARIANT_QC_FILE}.{pop}.tsv.bgz')
	mt = mt.filter_rows(mt.in_target)
	mt = hl.sample_qc(mt, name='sample_qc')
	ht = mt.cols().select('sample_qc').flatten()
	ht.export(f'{SAMPLE_QC_TARGET_FILE}.{pop}.tsv.bgz')

