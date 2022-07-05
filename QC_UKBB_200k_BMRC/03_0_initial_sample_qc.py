import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes

# Create sample QC metrics restricted and not restricted (target plus padding) 
# to the target intervals.

parser = argparse.ArgumentParser()
parser.add_argument("--chr", type=str, default='20')
parser.add_argument("--tranche", type=str, default='200k')
args = parser.parse_args()

TRANCHE = args.tranche
CHR = str(args.chr)

# Inputs
MT = '/well/lindgren/UKBIOBANK/nbaya/wes_' + TRANCHE + '/ukb_wes_qc/data/filtered/ukb_wes_' + TRANCHE + '_filtered_chr' + CHR + '.mt'
INITIAL_VARIANT_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/variants/02_prefilter_chr' + CHR +'.keep.variant.ht'
TARGET_INTERVALS = '/well/lindgren-ukbb/projects/ukbb-11867/nbaya/resources/ref/xgen_plus_spikein.b38.chr_prefix.bed'

# Outputs
INITIAL_SAMPLE_QC_FILE = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/03_chr' + CHR + '_initial_sample_qc.tsv.bgz'

print("Inputs:")
print('MT_HARDCALLS; input matrix table: ', MT_HARDCALLS)
print('INITIAL_VARIANT_LIST; input initial list of variants: ', INITIAL_VARIANT_LIST)
print('TARGET_INTERVALS; target intervals file: ', TARGET_INTERVALS)

print("Outputs:")
print('INITIAL_SAMPLE_QC_FILE; output variant QC file: ', INITIAL_SAMPLE_QC_FILE)

hail_init.hail_bmrc_init('logs/hail/hail_export.log', REFERENCE)

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
