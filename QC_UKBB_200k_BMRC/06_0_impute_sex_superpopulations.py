import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes

parser = argparse.ArgumentParser()
parser.add_argument("--tranche", type=str, default='200k')
args = parser.parse_args()

TRANCHE = args.tranche

hail_init.hail_bmrc_init_local('logs/hail/hail_export.log', 'GRCh38')

# Genotypes

# Inputs
# UKB genotype calls
ukb_bed_X="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_cal_chrX_v2.bed"
ukb_bim_X="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_snp_chrX_v2.bim"
ukb_bed_Y="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_cal_chrY_v2.bed"
ukb_bim_Y="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_snp_chrY_v2.bim"
ukb_fam="/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"
INITIAL_SAMPLES = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/03_initial_qc.keep.sample_list'
PHENOTYPES_TABLE = '/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/' + TRANCHE + '/phenotypes.ht'
PRUNED_CHRX_VARIANTS = '/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_snp_chrX_pruned.prune.in'
SUPERPOPS = "/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/superpopulation_labels.tsv"

# Outputs
IMPUTESEX_TABLE = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/04_imputesex_'
IMPUTESEX_FILE = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/04_imputesex_'
Y_NCALLED = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/04_ycalled.tsv.bgz'

ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
ht_pruned_chrx_variants = hl.import_table(PRUNED_CHRX_VARIANTS, no_header=True)
ht_pruned_chrx_variants = ht_pruned_chrx_variants.transmute(rsid=ht_pruned_chrx_variants.f0).key_by('rsid')
sample_annotations = hl.read_table(PHENOTYPES_TABLE)

# Read in the plink file
mt = hl.import_plink(bed=ukb_bed_X, bim=ukb_bim_X, fam=ukb_fam, reference_genome='GRCh37')
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
	mt_tmp = mt_tmp.annotate_cols(phenotype = sample_annotations[mt_tmp.s])
	mt_tmp = mt_tmp.annotate_cols(impute_sex = imputed_sex[mt_tmp.s])
	IMPUTESEX_FILE_tmp = IMPUTESEX_FILE + pop + '.tsv.bgz'
	IMPUTESEX_TABLE_tmp = IMPUTESEX_TABLE + pop + '.ht'
	mt_tmp.cols().select('impute_sex', 'phenotype').flatten().export(IMPUTESEX_FILE_tmp)
	mt_tmp.cols().write(IMPUTESEX_TABLE_tmp, overwrite=True)

# Now, look on the Y chromosome, and determine non-missing allele count on the y.
mt = hl.import_plink(bed=ukb_bed_Y, bim=ukb_bim_Y, fam=ukb_fam, reference_genome='GRCh37')
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_rows(mt.locus.in_y_nonpar() | mt.locus.in_y_par())
mt = hl.sample_qc(mt, name='qc')

mt_cols = mt.cols()
mt_cols.select(n_called=mt_cols.qc.n_called).export(Y_NCALLED)
