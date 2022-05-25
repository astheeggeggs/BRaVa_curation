#!/usr/bin/env bash
#$ -N hail_shell
#$ -cwd
#$ -o /well/lindgren/dpalmer/logs/hail.log
#$ -e /well/lindgren/dpalmer/logs/hail.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 30
#$ -q short.qf@@short.hgf
#$ -t 1-23

set -o errexit
set -o nounset

module purge
source /well/lindgren/dpalmer/ukb_utils/bash/qsub_utils.sh
source /well/lindgren/dpalmer/ukb_utils/bash/hail_utils.sh

module use -a /apps/eb/testing/${MODULE_CPU_TYPE}/modules/all
module load Hail/0.2.93-foss-2021b

_mem=$( get_hail_memory )
new_spark_dir='/well/lindgren/UKBIOBANK/dpalmer/tmp'
export PYSPARK_SUBMIT_ARGS="--conf spark.local.dir=${new_spark_dir} --conf spark.executor.heartbeatInterval=1000000 --conf spark.network.timeout=1000000  --driver-memory ${_mem}g --executor-memory ${_mem}g pyspark-shell"
export PYTHONPATH="${PYTHONPATH-}:/well/lindgren/dpalmer/ukb_utils/python:/well/lindgren/dpalmer:/well/lindgren/dpalmer/ukb_common/src"

chr=$(get_chr ${SGE_TASK_ID})

export HAIL_TMP_DIR="/well/lindgren/UKBIOBANK/dpalmer"

TRANCHE='200k'
MT="/well/lindgren/UKBIOBANK/dpalmer/wes_${TRANCHE}/ukb_wes_qc/data/filtered/ukb_wes_${TRANCHE}_filtered_hardcalls_chr${chr}.mt"
IMPUTESEX_TABLE="/well/lindgren/UKBIOBANK/dpalmer/wes_${TRANCHE}/ukb_wes_qc/data/samples/04_imputesex_"
SUPERPOPS="/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/superpopulation_labels.tsv"
SEXCHECK_LIST="/well/lindgren/UKBIOBANK/dpalmer/wes_${TRANCHE}/ukb_wes_qc/data/samples/04_sexcheck.remove.sample_list"
RELATED_SAMPLES="/well/lindgren/UKBIOBANK/dpalmer/wes_${TRANCHE}/ukb_wes_qc/data/samples/06_king.related.sample_list"
INITIAL_VARIANT_LIST="/well/lindgren/UKBIOBANK/dpalmer/wes_${TRANCHE}/ukb_wes_qc/data/variants/02_prefilter_chr${chr}.keep.variant.ht"
SAMPLE_LIST_INITIAL_QC="/well/lindgren/UKBIOBANK/dpalmer/wes_${TRANCHE}/ukb_wes_qc/data/samples/03_initial_qc.keep.sample_list"
VARIANT_QC_FILE="/well/lindgren/UKBIOBANK/dpalmer/wes_${TRANCHE}/ukb_wes_qc/data/variants/08_final_qc.variants_chr${chr}_"

python 08_final_variant_qc.py ${MT} ${IMPUTESEX_TABLE} ${SUPERPOPS} ${SEXCHECK_LIST} ${RELATED_SAMPLES} ${INITIAL_VARIANT_LIST} ${SAMPLE_LIST_INITIAL_QC} ${VARIANT_QC_FILE}
print_update "Finished running Hail for chr${chr}" "${SECONDS}"
