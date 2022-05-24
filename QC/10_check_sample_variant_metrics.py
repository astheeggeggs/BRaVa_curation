import hail as hl
from ukb_utils import hail_init
from ukb_utils import genotypes

hail_init.hail_bmrc_init('logs/hail/hail_export.log', 'GRCh38')

# Inputs
CURATED_MT = sys.argv[1]

# Outputs
VARIANT_QC_FILE = sys.argv[2]
SAMPLE_QC_FILE = sys.argv[3]

print("Inputs:")
print('MT; input matrix table: ', MT)

print("Outputs:")
print('VARIANT_QC_FILE; output .tsv file variant QC information: ', VARIANT_QC_FILE)
print('SAMPLE_QC_FILE; output .tsv file variant QC information: ', SAMPLE_QC_FILE)

mt = hl.read_matrix_table(CURATED_MT)

mt = hl.variant_qc(mt, name='variant_qc')
mt = hl.sample_qc(mt, name='sample_qc')

# Spit out the resultant files

ht = mt.cols().select('sample_qc').flatten()
ht.export(SAMPLE_QC_FILE)
ht = mt.rows().select('variant_qc').flatten()
ht.export(VARIANT_QC_FILE)
