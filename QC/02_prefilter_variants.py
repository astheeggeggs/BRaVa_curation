import hail as hl
import sys

# Perform an initial variant filter to remove variants outside the target intervals, 
# within low-complexity regions, and that fail VQSR.

# Inputs
MT_HARDCALLS = sys.argv[1]
TARGET_INTERVALS = sys.argv[2]
PADDED_TARGET_INTERVALS = sys.argv[3]
LCRs = sys.argv[4]
REFERENCE = 'GRCh38'

# Outputs
INITIAL_VARIANT_QC_FILE  = sys.argv[5]
INITIAL_VARIANT_LIST = sys.argv[6]

print("Inputs:")
print('MT_HARDCALLS; input hard calls matrix table: ', MT_HARDCALLS)
print('TARGET_INTERVALS; target intervals file: ', TARGET_INTERVALS)
print('PADDED_TARGET_INTERVALS; padded target intervals file: ', PADDED_TARGET_INTERVALS)
print('LCRs; low-complexity regions intervals file: ', LCRs)

print("Outputs:")
print('INITIAL_VARIANT_QC_FILE; .tsv of variant QC metrics after filtering to variants \
	in the padded target interval, not in a low-complexity region, and that pass VQSR: ',
	INITIAL_VARIANT_QC_FILE)
print('INITIAL_VARIANT_LIST; set of variants that pass this initial filtering: ', INITIAL_VARIANT_LIST)

hl.init(default_reference=REFERENCE)

# Names of .mt files.
# Read in the hard calls matrix table.
mt = hl.read_matrix_table(MT_HARDCALLS)

# Import the interval lists for the target intervals.
target_intervals = hl.import_locus_intervals(TARGET_INTERVALS, reference_genome=REFERENCE)
# Import the interval lists for the padded target intervals.
padded_target_intervals = hl.import_locus_intervals(PADDED_TARGET_INTERVALS, reference_genome=REFERENCE)
# Import the interval lists for the LCRs.
LCR_intervals = hl.import_locus_intervals(LCRs, reference_genome=REFERENCE)

# Annotate variants with flag indicating if they are in LCR or failed VQSR.
mt = mt.annotate_rows(fail_VQSR = hl.len(mt.filters) != 0) # DEV: Note that if GATK was not used, VQSR information will not be available.
mt = mt.annotate_rows(in_LCR = hl.is_defined(LCR_intervals[mt.locus]))
mt = mt.annotate_rows(not_in_target_intervals = ~hl.is_defined(target_intervals[mt.locus]))
mt = mt.annotate_rows(not_in_padded_target_intervals = ~hl.is_defined(padded_target_intervals[mt.locus]))

# Get information about the number of variants that were excluded.
fail_VQSR = mt.filter_rows(mt.fail_VQSR).count_rows()
in_LCR = mt.filter_rows(mt.in_LCR).count_rows()
not_in_target_intervals = mt.filter_rows(mt.not_in_target_intervals).count_rows()
not_in_padded_target_intervals = mt.filter_rows(mt.not_in_padded_target_intervals).count_rows()

print('n variants failing VQSR:')
print(fail_VQSR)
print('n variants in low complexity regions:')
print(in_LCR)
print('n variants not in target intervals:')
print(not_in_target_intervals)
print('n variants not in padded target intervals:')
print(not_in_padded_target_intervals)

# Export variant annotations, and variant keytable.
mt_rows = mt.rows()
mt_rows.select(mt_rows.fail_VQSR, mt_rows.in_LCR, mt_rows.not_in_padded_target_intervals).export(INITIAL_VARIANT_QC_FILE)
mt = mt.filter_rows(mt.fail_VQSR | mt.in_LCR | mt.not_in_padded_target_intervals, keep=False)

intervals = [hl.parse_locus_interval(x, reference_genome=REFERENCE) for x in ['chr1:START-chr22:END', 'chrX:START-chrX:END', 'chrY:START-chrY:END']]
mt = hl.filter_intervals(mt, intervals)

# Filter out the invariant rows.
mt = hl.variant_qc(mt, name='qc')
mt = mt.filter_rows((mt.qc.AF[0] > 0.0) & (mt.qc.AF[0] < 1.0))

mt_rows_filter = mt.rows().select().export(INITIAL_VARIANT_LIST)

n_variants = hl.import_table(INITIAL_VARIANT_LIST).count()

print('n variants after initial filter:')
print(n_variants)
