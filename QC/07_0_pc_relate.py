import hail as hl
import sys

# Inputs
MT_HARDCALLS = sys.argv[1]
INITIAL_SAMPLES = sys.argv[2]
PRUNED_VARIANTS = sys.argv[3]
PHENOTYPES_TABLE = sys.argv[4]
REFERENCE = 'GRCh38'
PC_RELATE_KINSHIP_THRESHOLD = 0.08838835

# Outputs
PC_RELATE_OUTPUT = sys.argv[5]
SAMPLE_LIST_RELATED = sys.argv[6]

print("Inputs:")
print('MT_HARDCALLS; input hard calls matrix table: ', MT_HARDCALLS)
print('INITIAL_SAMPLES; set of initial samples output from 03_01_initial_sample_qc_filter.r: ', INITIAL_SAMPLES)
print('PRUNED_VARIANTS; set of LD-pruned high-quality variants in the autosomes: ', PRUNED_VARIANTS)

print("Outputs:")
print('PC_RELATE_OUTPUT; output .tsv file with pc-relate information for plotting: ', PC_RELATE_OUTPUT)

hl.init(default_reference=REFERENCE)

ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
ht_pruned_variants = hl.import_table(PRUNED_VARIANTS, no_header=True)

ht_pruned_variants = ht_pruned_variants.annotate(**hl.parse_variant(ht_pruned_variants.f0, reference_genome=REFERENCE))
ht_pruned_variants = ht_pruned_variants.key_by(ht_pruned_variants.locus, ht_pruned_variants.alleles)

mt = hl.read_matrix_table(MT_HARDCALLS)
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_rows(hl.is_defined(ht_pruned_variants[mt.row_key]))
mt = mt.repartition(128).persist()

n = mt.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

# Not tested, if this is slow or doesn't work well - consider passing the projected PCs from step 05_estimate_superpopulation.r
pc_rel = hl.pc_relate(dataset.GT, 0.01, k=4, statistics='kin')
pairs = pc_rel.filter(pc_rel['kin'] > PC_RELATE_KINSHIP_THRESHOLD)

related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, False)
result = mt.cols().select()
result.filter(hl.is_defined(related_samples_to_remove[mt.col_key]), keep=True).export(SAMPLE_LIST_RELATED)

pc_relate_table = pc_relate_table.flatten()
pc_relate_table.export(PC_RELATE_OUTPUT)
