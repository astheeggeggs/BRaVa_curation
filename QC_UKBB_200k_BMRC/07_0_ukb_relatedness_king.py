import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes

parser = argparse.ArgumentParser()
parser.add_argument("--tranche", type=str, default='200k')
args = parser.parse_args()

TRANCHE = args.tranche

# Inputs
INITIAL_SAMPLES='/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/03_initial_qc.keep.sample_list'
KING_RELATEDS_FILE="/well/lindgren-ukbb/projects/ukbb-11867/DATA/QC/ukb1186_rel_s488366.dat"
SAMPLE_LIST_RELATED = '/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/07_king.related.sample_list'

hail_init.hail_bmrc_init('logs/hail/hail_export.log', 'GRCh38')

KING_KINSHIP_THRESHOLD = 0.08838835

# Outputs
SAMPLE_LIST_RELATED = sys.argv[3]

print("Inputs:")
print('INITIAL_SAMPLES; initial samples following QC: ', INITIAL_SAMPLES)
print('KING_RELATEDS_FILE; UK Biobank provided file of up to 3rd degree relateds: ', KING_RELATEDS_FILE)

print("Outputs:")
print('SAMPLE_LIST_RELATED; output file of related samples for exclusion: ', SAMPLE_LIST_RELATED)

ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')

# Read in the KING file and merge with the initial set of samples
ht_related_samples_king = hl.import_table(
	"/well/lindgren-ukbb/projects/ukbb-11867/DATA/QC/ukb1186_rel_s488366.dat",
	delimiter=' ', impute=True, types = {"ID1":hl.tstr, "ID2":hl.tstr}, key="ID1"
	)

ht_related_samples_king = ht_related_samples_king.annotate(
	ID1_exome = hl.is_defined(ht_initial_samples[ht_related_samples_king.key]))
ht_related_samples_king = ht_related_samples_king.key_by('ID2')
ht_related_samples_king = ht_related_samples_king.annotate(
	ID2_exome = hl.is_defined(ht_initial_samples[ht_related_samples_king.key]))

ht_related_samples_king = ht_related_samples_king.filter((ht_related_samples_king.ID1_exome) & (ht_related_samples_king.ID2_exome))
pairs=ht_related_samples_king.filter(ht_related_samples_king['Kinship'] > KING_KINSHIP_THRESHOLD)
related_samples_to_remove = hl.maximal_independent_set(pairs.ID1, pairs.ID2, False)
related_samples_to_remove.export(SAMPLE_LIST_RELATED)
