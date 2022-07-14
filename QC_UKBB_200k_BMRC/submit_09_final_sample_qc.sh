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

module load Anaconda3/2020.07
module load java/1.8.0_latest
source activate hail-new
_mem=$( get_hail_memory )
new_spark_dir=/well/lindgren/dpalmer/tmp/spark_test/
export PYSPARK_SUBMIT_ARGS="--conf spark.local.dir=${new_spark_dir} --conf spark.executor.heartbeatInterval=1000000 --conf spark.network.timeout=1000000  --driver-memory ${_mem}g --executor-memory ${_mem}g pyspark-shell"
export PYTHONPATH="${PYTHONPATH-}:/well/lindgren/dpalmer/ukb_utils/python:/well/lindgren/dpalmer:/well/lindgren/dpalmer/ukb_common/src"

chr=$(get_chr ${SGE_TASK_ID})

export HAIL_TMP_DIR="/well/lindgren/UKBIOBANK/dpalmer/logs"
export HAIL_SPARK_DIR="/well/lindgren/UKBIOBANK/dpalmer/logs/spark"

TRANCHE="200k"

python 09_final_sample_qc.py --chr ${chr} --tranche ${TRANCHE}
print_update "Finished running Hail for chr${chr}" "${SECONDS}"

