import hail as hl
import sys

# Inputs
INPUT_VCF = sys.argv[1]
REFERENCE = 'GRCh38'

# Outputs
OUTPUT_MATRIX_TABLE = sys.argv[2]

print("Inputs:")
print('INPUT_VCF; vcf input: ', INPUT_VCF)

print("Outputs:")
print('RAW_MT; output matrix table: ', RAW_MT)

hl.init(default_reference=REFERENCE)

mt = hl.import_vcf(INPUT_VCF, reference_genome=REFERENCE, force_bgz=True, find_replace=('nul', '.'))
mt.write(output=RAW_MT, overwrite=True)
