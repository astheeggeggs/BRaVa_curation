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

TRANCHE="200k"

# Inputs:
QC_MT_PREFIX="/well/lindgren/UKBIOBANK/dpalmer/wes_${TRANCHE}/ukb_wes_qc/data/final_mt/10_strict_filtered_chr"

# Outputs
QC_VCF_PREFIX="/well/lindgren/UKBIOBANK/dpalmer/wes_${TRANCHE}/ukb_wes_qc/data/final_mt/10_strict_filtered_chr"

python 10_create_qc_mt.py ${QC_MT_PREFIX} ${QC_VCF_PREFIX}
print_update "Finished running Hail for chr${chr}" "${SECONDS}"
