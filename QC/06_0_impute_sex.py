import hail as hl
import sys

# Inputs
MT_HARDCALLS = sys.argv[1]
SAMPLE_LIST_INITIAL_QC= sys.argv[2]
PRUNED_CHRX_VARIANTS = sys.argv[3]
PHENOTYPES_TABLE = sys.argv[4]
REFERENCE = 'GRCh38'

# Outputs
IMPUTESEX_TABLE = sys.argv[5]
IMPUTESEX_FILE = sys.argv[6]
Y_NCALLED = sys.argv[7]

print("Inputs:")
print('MT_HARDCALLS; input hard calls matrix table: ', MT_HARDCALLS)
print('SAMPLE_LIST_INITIAL_QC; set of initial samples output from 03_01_initial_sample_qc_filter.r: ', SAMPLE_LIST_INITIAL_QC)
print('PRUNED_CHRX_VARIANTS; set of LD-pruned high-quality variants on the X: ', PRUNED_CHRX_VARIANTS)
print('PHENOTYPES_TABLE; a collection of sample annotations, which includes self-assigned sex or gender: ', PHENOTYPES_TABLE)

print("Outputs:")
print('IMPUTESEX_TABLE; output .tsv file to plot imputed sex information : ', IMPUTESEX_TABLE)
print('IMPUTESEX_FILE; output .tsv file to use to count mismatches at this step for a summary table: ', IMPUTESEX_FILE)
print('Y_NCALLED; output .tsv file of calls on the Y: ', Y_NCALLED)

hl.init(default_reference=REFERENCE)

ht_initial_samples = hl.import_table(SAMPLE_LIST_INITIAL_QC, no_header=True, key='f0')
ht_pruned_chrx_variants = hl.import_table(PRUNED_CHRX_VARIANTS, no_header=True)
sample_annotations = hl.read_table(PHENOTYPES_TABLE)

ht_pruned_chrx_variants = ht_pruned_chrx_variants.annotate(**hl.parse_variant(ht_pruned_chrx_variants.f0, reference_genome=REFERENCE))
ht_pruned_chrx_variants = ht_pruned_chrx_variants.key_by(ht_pruned_chrx_variants.locus, ht_pruned_chrx_variants.alleles)

mt = hl.read_matrix_table(MT_HARDCALLS)
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_rows(hl.is_defined(ht_pruned_chrx_variants[mt.row_key]))

n = mt.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

imputed_sex = hl.impute_sex(mt.GT, female_threshold=0.6, male_threshold=0.6)
mt = mt.annotate_cols(phenotype = sample_annotations[mt.s])
mt = mt.annotate_cols(impute_sex = imputed_sex[mt.s])

mt.cols().select('impute_sex', 'phenotype').flatten().export(IMPUTESEX_FILE)
mt.cols().write(IMPUTESEX_TABLE, overwrite=True)

# Determine non-missing allele count on the y.
mt = hl.read_matrix_table(MT_HARDCALLS)
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_rows(mt.locus.in_y_nonpar() | mt.locus.in_y_par())
mt = hl.sample_qc(mt, name='qc')

mt_cols = mt.cols()
mt_cols.select(n_called=mt_cols.qc.n_called).export(Y_NCALLED)
