import hail as hl
import sys

# Inputs
MT_HARDCALLS = sys.argv[1]
INITIAL_SAMPLES = sys.argv[2]
INITIAL_VARIANT_LIST = sys.argv[3]
HIGH_LD_INTERVALS = sys.argv[4]
REFERENCE = 'GRCh38'

# Outputs
PLINK_FILES = sys.argv[5]

print("Inputs:")
print('MT_HARDCALLS; input hard calls matrix table: ', MT_HARDCALLS)
print('INITIAL_SAMPLES; target intervals file: ', TARGET_INTERVALS)
print('INITIAL_VARIANT_LIST; padded target intervals file: ', PADDED_TARGET_INTERVALS)
print('HIGH_LD_INTERVALS; set of high-LD regions for removal: ', LCRs)

print("Outputs:")
print('PLINK_FILES; prefix for high-quality common variant plink files: ', PLINK_FILES)

hl.init(default_reference=REFERENCE)

ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
ht_initial_variants = hl.import_table(INITIAL_VARIANT_LIST, types={'locus':hl.tlocus(reference_genome=REFERENCE), 'alleles':hl.tarray(hl.tstr)})
ht_initial_variants = ht_initial_variants.key_by(ht_initial_variants.locus, ht_initial_variants.alleles)

high_LD_intervals = hl.import_locus_intervals(HIGH_LD_INTERVALS, reference_genome=REFERENCE)

mt = hl.read_matrix_table(MT_HARDCALLS)
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.annotate_rows(in_high_LD = hl.is_defined(high_LD_intervals[mt.locus]))

mt = mt.filter_rows(hl.is_defined(ht_initial_variants[mt.row_key]) & (~mt.in_high_LD))
mt = mt.filter_rows(mt.locus.in_x_nonpar() | mt.locus.in_autosome_or_par())
mt = hl.variant_qc(mt, name='qc')
mt = mt.filter_rows(
	(mt.qc.AF[0] > 0.01) & 
	(mt.qc.AF [0]< 0.99) & 
	((mt.qc.call_rate > 0.98) | mt.locus.in_x_nonpar() | mt.locus.in_x_par())
	).persist()

mt.count()

for x in range(1,23):

	mt_chr = hl.filter_intervals(
		mt, [hl.parse_locus_interval(hl.eval('chr' + hl.str(x)), reference_genome=REFERENCE)])
	n_chr = mt_chr.count_rows()

	print('\nn variants in chr')
	print(x)
	print(n_chr)

	hl.export_plink(mt_chr, PLINK_FILES + '.chr' + str(x))

mt_chr = hl.filter_intervals(mt, [hl.parse_locus_interval('chrX', reference_genome=REFERENCE)])
n_chr = mt_chr.count_rows()

print('\nn variants in chrX')
print(n_chr)

hl.export_plink(mt_chr, PLINK_FILES + '.chr' + 'X')
