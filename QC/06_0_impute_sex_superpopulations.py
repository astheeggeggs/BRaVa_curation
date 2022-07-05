import hail as hl
import sys

# Inputs
# Plink files - either from sequencing or genotype data
BIM_X = sys.argv[1]
BED_X = sys.argv[2]
FAM_X = sys.argv[3]
BIM_Y = sys.argv[4]
BED_Y = sys.argv[5]
FAM_Y = sys.argv[6]
INITIAL_SAMPLES = sys.argv[7]
PRUNED_CHRX_VARIANTS = sys.argv[8]
SUPERPOPS = sys.argv[9]

# Outputs
IMPUTESEX_TABLE = sys.argv[10]
IMPUTESEX_FILE = sys.argv[11]
Y_NCALLED = sys.argv[12]

REFERENCE = 'GRCh37' # Note that the reference genome could be different for genotyping data vs sequencing data due to age of technology - watch out!

print("Inputs:")
print("BIM_X; X chromrome plink bim file: ", BIM_X)
print("BED_X; X chromrome plink bim file: ", BED_X)
print("FAM_X; X chromrome plink bim file: ", FAM_X)
print("BIM_Y; Y chromrome plink bim file: ", BIM_Y)
print("BED_Y; Y chromrome plink bim file: ", BED_Y)
print("FAM_Y; Y chromrome plink bim file: ", FAM_Y)
print("INITIAL_SAMPLES; initial samples following QC: ", INITIAL_SAMPLES)
print("PRUNED_CHRX_VARIANTS; pruned chromosome X variants for sex imputation: ", PRUNED_CHRX_VARIANTS)
print("SUPERPOPS; 1000G labels file: ", SUPERPOPS)

print("Outputs:")
print("IMPUTESEX_TABLE; imputed sex hail table output prefix: ", IMPUTESEX_TABLE)
print("IMPUTESEX_FILE; imputed sex tsv.bgz file output prefix: ", IMPUTESEX_FILE)
print("Y_NCALLED; calls on y tsv.bgz file output prefix: ", Y_NCALLED)

hail_init(default_reference=REFERENCE)

ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
ht_pruned_chrx_variants = hl.import_table(PRUNED_CHRX_VARIANTS, no_header=True)
ht_pruned_chrx_variants = ht_pruned_chrx_variants.transmute(rsid=ht_pruned_chrx_variants.f0).key_by('rsid')

# Read in the plink file
mt = hl.import_plink(bed=BED_X, bim=BIM_X, fam=FAM_X, reference_genome=REFERENCE)
mt = mt.key_rows_by(mt.rsid)
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_rows(hl.is_defined(ht_pruned_chrx_variants[mt.row_key]))
mt = mt.key_rows_by(mt.locus, mt.alleles)
n = mt.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

ht_superpops = hl.import_table(SUPERPOPS, impute=True).key_by("sample.ID").select("classification_strict")
ht_pruned_chrx_variants = hl.import_table(PRUNED_CHRX_VARIANTS, no_header=True)

mt = mt.annotate_cols(pops = ht_superpops[mt.s])

for pop in ["AFR", "AMR", "EAS", "EUR", "SAS"]:
	mt_tmp = mt.filter_cols(mt.pops.classification_strict == pop)
	n = mt_tmp.count()
	print('n samples:')
	print(n[1])
	print('n variants:')
	print(n[0])
	imputed_sex = hl.impute_sex(mt_tmp.GT, female_threshold=0.2, male_threshold=0.8)
	mt_tmp = mt_tmp.annotate_cols(impute_sex = imputed_sex[mt_tmp.s])
	IMPUTESEX_FILE_tmp = IMPUTESEX_FILE + pop + '.tsv.bgz'
	IMPUTESEX_TABLE_tmp = IMPUTESEX_TABLE + pop + '.ht'
	mt_tmp.cols().select('impute_sex').flatten().export(IMPUTESEX_FILE_tmp)
	mt_tmp.cols().write(IMPUTESEX_TABLE_tmp, overwrite=True)

# Now, look on the Y chromosome, and determine non-missing allele count on the y.
mt = hl.import_plink(bed=BED_Y, bim=BIM_Y, fam=FAM_Y, reference_genome=REFERENCE)
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_rows(mt.locus.in_y_nonpar() | mt.locus.in_y_par())
mt = hl.sample_qc(mt, name='qc')

mt_cols = mt.cols()
mt_cols.select(n_called=mt_cols.qc.n_called).export(Y_NCALLED)
