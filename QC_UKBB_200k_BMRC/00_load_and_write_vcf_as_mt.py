# Note that in our filtering steps, we did run this exact script to create a raw matrix table to build off of, 
# but something very similar as part of Nik's pipeline. We include this simple script for guidance.

import hail as hl
import sys

parser.add_argument("--chr", type=str, default='20')
parser.add_argument("--tranche", type=str, default='200k')
parser.add_argument("--raw_mt", type=str)
parser.add_argument("--raw_vcf", type=str)
args = parser.parse_args()

# Inputs
TRANCHE = args.tranche
CHR = str(args.chr)
RAW_VCF = args.raw_vcf + '.' + CHR + '.vcf.bgz'

# Outputs
RAW_MT = args.raw_mt + '.' + CHR + '.mt'

print("Inputs:")
print('INPUT_VCF; vcf input: ', INPUT_VCF)

print("Outputs:")
print('RAW_MT; output matrix table: ', RAW_MT)

hl.init(default_reference=REFERENCE)

mt = hl.import_vcf(RAW_VCF, reference_genome=REFERENCE, force_bgz=True, find_replace=('nul', '.'))
mt.write(output=RAW_MT, overwrite=True)
