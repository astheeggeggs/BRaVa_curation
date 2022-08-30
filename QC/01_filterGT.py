import hail as hl
import sys

# Inputs
RAW_MT = sys.argv[1]
N_PARTITIONS = 512
REFERENCE = 'GRCh38'

# Outputs
MT = sys.argv[2]
MT_HARDCALLS = sys.argv[3]

print("Inputs:")
print('RAW_MT; matrix table input: ', RAW_MT)

print("Outputs:")
print('MT; output matrix table: ', MT)
print('MT_HARDCALLS; output hard calls matrix table: ', MT_HARDCALLS)

hl.init(default_reference=REFERENCE)

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

mt = mt.checkpoint(MT, overwrite=True)
mt = hl.read_matrix_table(MT)
mt.select_entries(mt.GT).repartition(N_PARTITIONS).write(MT_HARDCALLS, overwrite=True)

mt = hl.read_matrix_table(MT_HARDCALLS)
n = mt.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])
