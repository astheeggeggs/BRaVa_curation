import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

# 01_load_vcf_filterGT.py

RAW_MT = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/dalio_bipolar_w1_w2/Dalio_W1_W2_GRCh38_exomes.mt'
mt = hl.read_matrix_table(RAW_MT)

# Initial count.
n_raw = mt.count_rows()

print('Initial number of variants:')
pprint(n_raw)

# Number with less that or equal to 6 alleles.
MT_HARDCALLS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT_GRCh38_6_multi.hardcalls.mt'
mt = hl.read_matrix_table(MT_HARDCALLS)
n_6_multi = mt.count_rows()

print('n variants with <= 6 alleles:')
pprint(n_6_multi)

# 02_prefilter_variants.py
INITIAL_VARIANT_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/02_prefilter.keep.variant_list'
INITIAL_VARIANT_QC_FILE  = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/02_prefilter_metrics.tsv'

# Read in the target intervals
TARGET_INTERVALS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/inputs/ice_coding_v1_targets.interval_list'
# Read in the padded target intervals (50bp padding)
PADDED_TARGET_INTERVALS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/inputs/ice_coding_v1_padded_targets.interval_list'

# Low complexity regions in the data.
LCRs = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/inputs/LCR-hs38.bed'

# Import the interval lists for the target intervals.
target_intervals = hl.import_locus_intervals(TARGET_INTERVALS, reference_genome='GRCh38')
# Import the interval lists for the padded target intervals.
padded_target_intervals = hl.import_locus_intervals(PADDED_TARGET_INTERVALS, reference_genome='GRCh38')
# Import the interval lists for the LCRs.
LCR_intervals = hl.import_locus_intervals(LCRs, reference_genome='GRCh38')

# Annotate variants with flag indicating if they are in LCR or failed VQSR.
mt = mt.annotate_rows(fail_VQSR = hl.len(mt.filters) != 0)
mt = mt.annotate_rows(in_LCR = hl.is_defined(LCR_intervals[mt.locus]))
mt = mt.annotate_rows(not_in_target_intervals = ~hl.is_defined(target_intervals[mt.locus]))
mt = mt.annotate_rows(not_in_padded_target_intervals = ~hl.is_defined(padded_target_intervals[mt.locus]))

# Get information about the number of variants that were excluded.
fail_VQSR = mt.filter_rows(mt.fail_VQSR).count_rows()
in_LCR = mt.filter_rows(mt.in_LCR).count_rows()
not_in_target_intervals = mt.filter_rows(mt.not_in_target_intervals).count_rows()
not_in_padded_target_intervals = mt.filter_rows(mt.not_in_padded_target_intervals).count_rows()

print('n variants failing VQSR:')
pprint(fail_VQSR)
print('n variants in low complexity regions:')
pprint(in_LCR)
print('n variants not in padded target intervals:')
pprint(not_in_padded_target_intervals)

# Filter to variants in the autosomes, chrX and Y.
mt_rows = mt.rows()
mt = mt.filter_rows(mt.fail_VQSR | mt.in_LCR | mt.not_in_padded_target_intervals, keep=False)

intervals = [hl.parse_locus_interval(x, reference_genome='GRCh38') for x in ['chr1:START-chr22:END', 'chrX:START-chrX:END', 'chrY:START-chrY:END']]
n_after_filter = mt.count_rows()

mt = hl.filter_intervals(mt, intervals)

n_in_chr = mt.count_rows()
print('n variants in autosomes, chrX, chrY:')
pprint(n_in_chr)

print('n variants removed in contigs outside autosomes, X and Y:')
pprint(n_after_filter - n_in_chr)

# Filter to variants in the autosomes, chrX and Y.
n_initial_variant_list = hl.import_table(INITIAL_VARIANT_LIST).count()
print('Invariant sites after initial variant and genotype filters:')
pprint(n_initial_variant_list)

# Create an initial markdown table
with hl.hadoop_open('gs://dalio_bipolar_w1_w2_hail_02/data/summary_variant_table.tsv', 'w') as f:
	f.write(hl.eval('Initial variants\t' + hl.str(n_raw) + '\n'))
	f.write(hl.eval('Variants with > 6 alleles\t' + hl.str(n_6_multi) + '\n'))
	f.write(hl.eval('Failing VQSR\t' + hl.str(fail_VQSR) + '\n'))
	f.write(hl.eval('In LCRs\t' + hl.str(in_LCR) + '\n'))
	f.write(hl.eval('Outside padded target interval\t' + hl.str(not_in_padded_target_intervals) + '\n'))
	f.write(hl.eval('Variants after initial filter\t' + hl.str(n_after_filter) + '\n'))
	f.write(hl.eval('Variants in contigs outside autosomes, X and Y\t' + hl.str(n_after_filter - n_in_chr) + '\n'))
	f.write(hl.eval('Invariant sites after initial variant and genotype filters\t' + hl.str(n_initial_variant_list) + '\n'))



