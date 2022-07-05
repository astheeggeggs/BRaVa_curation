# Note that in our filtering steps, we did not perform GT filtering in hail at this stage. This was performed 
# earlier in Nik's pipeline, but we include the steps here for clarity - genotype filtering should be applied!

import hail as hl
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--chr", type=str, default='20')
parser.add_argument("--tranche", type=str, default='200k')
parser.add_argument("--raw_mt", type=str)
args = parser.parse_args()

# Inputs
TRANCHE = args.tranche
CHR = str(args.chr)
RAW_MT = args.raw_mt
RAW_MT = RAW_MT + '.' + CHR + '.mt'
N_PARTITIONS = 512
REFERENCE = 'GRCh38'

# Outputs
MT  = '/well/lindgren/UKBIOBANK/nbaya/wes_' + TRANCHE + '/ukb_wes_qc/data/filtered/ukb_wes_' + TRANCHE + '_filtered_chr' + CHR + '.mt'
MT_HARDCALLS = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/filtered/ukb_wes_' + TRANCHE + '_filtered_hardcalls_chr' + CHR + '.mt'

print("Inputs:")
print('RAW_MT; matrix table input: ', RAW_MT)

print("Outputs:")
print('MT; output matrix table: ', MT)
print('MT_HARDCALLS; output hard calls matrix table: ', MT_HARDCALLS)

hail_init.hail_bmrc_init('logs/hail/hail_export.log', REFERENCE)

mt = hl.read_matrix_table(RAW_MT)

# Count before splitting multi-allelics.
n = mt.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

# Remove variants with a large number of alleles.
mt = mt.filter_rows(mt.alleles.length() <= 6)

n = mt.count_rows()

print('n variants not more than 6 alleles:')
print(n)

mt = hl.split_multi_hts(mt)

mt = mt.filter_entries(
    hl.is_defined(mt.GT) &
    (
        (mt.GT.is_hom_ref() & 
            (
                (mt.GQ < 20) |
                (mt.DP < 10)
        	)
        ) |
        (mt.GT.is_het() & 
        	( 
                (((mt.AD[0] + mt.AD[1]) / mt.DP) < 0.8) | 
                ((mt.AD[1] / mt.DP) < 0.2) | 
                (mt.PL[0] < 20) |
                (mt.DP < 10)
        	)
        ) |
        (mt.GT.is_hom_var() & 
        	(
                ((mt.AD[1] / mt.DP) < 0.8) |
                (mt.PL[0] < 20) |
                (mt.DP < 10)
        	)
        )
    ),
    keep = False
)

mt.write(MT, overwrite=True)
mt = mt.checkpoint(MT, overwrite=True)
mt = hl.read_matrix_table(MT)
mt.select_entries(mt.GT).repartition(N_PARTITIONS).write(MT_HARDCALLS, overwrite=True)

mt = hl.read_matrix_table(MT_HARDCALLS)
n = mt.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])
