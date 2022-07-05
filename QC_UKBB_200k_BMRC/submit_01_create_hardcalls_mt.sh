#!/usr/bin/env bash
#$ -N hail_shell
#$ -cwd
#$ -o /well/lindgren/dpalmer/logs/hail.log
#$ -e /well/lindgren/dpalmer/logs/hail.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 10
#$ -q short.qe@@short.hge
#$ -t 1-24

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

python 01_create_hardcalls_mt.py --chr ${chr}
print_update "Finished running Hail for chr${chr}" "${SECONDS}"
