import hail as hl
import sys

# Create sample QC metrics restricted and not restricted (target plus padding) 
# to the target intervals.

# Inputs
MT = sys.argv[1]
INITIAL_VARIANT_LIST = sys.argv[2]
TARGET_INTERVALS = sys.argv[3]
REFERENCE = 'GRCh38'

# Outputs
INITIAL_SAMPLE_QC_FILE  = sys.argv[4]

print("Inputs:")
print('MT_HARDCALLS; input matrix table: ', MT_HARDCALLS)
print('INITIAL_VARIANT_LIST; input initial list of variants: ', INITIAL_VARIANT_LIST)
print('TARGET_INTERVALS; target intervals file: ', TARGET_INTERVALS)

print("Outputs:")
print('INITIAL_SAMPLE_QC_FILE; output variant QC file: ', INITIAL_SAMPLE_QC_FILE)

hl.init(default_reference=REFERENCE)

variants_to_filter = hl.import_table(INITIAL_VARIANT_LIST,
	types={'locus':hl.tlocus(reference_genome=REFERENCE), 'alleles':hl.tarray(hl.tstr)})
variants_to_filter = variants_to_filter.key_by(locus=variants_to_filter.locus, alleles=variants_to_filter.alleles)

mt = hl.read_matrix_table(MT)
mt = mt.filter_rows(hl.is_defined(variants_to_filter[mt.row_key]))
mt = mt.annotate_cols(gq = hl.agg.stats(mt.GQ), dp = hl.agg.stats(mt.DP))
mt = hl.sample_qc(mt, name='qc_padded_target')

# Import the target interval lists.
target_intervals = hl.import_locus_intervals(TARGET_INTERVALS, reference_genome=REFERENCE)
mt = mt.annotate_rows(not_in_target_intervals = ~hl.is_defined(target_intervals[mt.locus]))
mt = mt.filter_rows(mt.not_in_target_intervals, keep=False)
mt = hl.sample_qc(mt, name='qc_target')

mt.cols().select('qc_padded_target', 'qc_target', 'gq', 'dp').flatten().export(output=INITIAL_SAMPLE_QC_FILE)
