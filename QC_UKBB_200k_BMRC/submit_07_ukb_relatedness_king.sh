#!/usr/bin/env bash
#$ -N hail_shell
#$ -cwd
#$ -o /well/lindgren/dpalmer/logs/hail.log
#$ -e /well/lindgren/dpalmer/logs/hail.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 8
#$ -q short.qe@@short.hge

set -o errexit
set -o nounset

cd /well/lindgren/dpalmer/SAIGE_gene_munging/QC_scripts
module purge
source /well/lindgren/dpalmer/ukb_utils/bash/qsub_utils.sh
source /well/lindgren/dpalmer/ukb_utils/bash/hail_utils.sh

module load Anaconda3/2020.07
module load java/1.8.0_latest
source activate hail-new
_mem=$( get_hail_memory )
new_spark_dir=/well/lindgren/dpalmer/tmp/
export PYSPARK_SUBMIT_ARGS="--conf spark.local.dir=${new_spark_dir} --conf spark.executor.heartbeatInterval=1000000 --conf spark.network.timeout=1000000  --driver-memory ${_mem}g --executor-memory ${_mem}g pyspark-shell"
# export PYSPARK_SUBMIT_ARGS="${PYSPARK_SUBMIT_ARGS} --conf spark.executor.extraClassPath='/usr/lib64/atlas'" # Added to avoid error: "undefined symbol: cblas_dgemv"
export PYTHONPATH="${PYTHONPATH-}:/well/lindgren/dpalmer/ukb_utils/python:/well/lindgren/dpalmer:/well/lindgren/dpalmer/ukb_common/src"

export HAIL_TMP_DIR="/well/lindgren/dpalmer"
export HAIL_SPARK_DIR="/well/lindgren/dpalmer/spark"

module load OpenBLAS/0.3.1-GCC-7.3.0-2.30
export LD_PRELOAD=/apps/eb/skylake/software/OpenBLAS/0.3.1-GCC-7.3.0-2.30/lib/libopenblas.so

python 06_0_pc_relate.py
print_update "Finished running Hail "${SECONDS}"


spark.executor.extraLibraryPath





module purge
source /well/lindgren/dpalmer/ukb_utils/bash/qsub_utils.sh
source /well/lindgren/dpalmer/ukb_utils/bash/hail_utils.sh

cd /well/lindgren/dpalmer/SAIGE_gene_munging/QC_scripts/
module use -a /apps/eb/testing/${MODULE_CPU_TYPE}/modules/all
module load Hail/0.2.93-foss-2021b

_mem=$( get_hail_memory )
new_spark_dir=/well/lindgren/dpalmer/tmp/
export PYSPARK_SUBMIT_ARGS="--conf spark.local.dir=${new_spark_dir} --conf spark.executor.heartbeatInterval=1000000 --conf spark.network.timeout=1000000  --driver-memory ${_mem}g --executor-memory ${_mem}g pyspark-shell"
export PYTHONPATH="${PYTHONPATH-}:/well/lindgren/dpalmer/ukb_utils/python:/well/lindgren/dpalmer:/well/lindgren/dpalmer/ukb_common/src"
export PYSPARK_SUBMIT_ARGS="${PYSPARK_SUBMIT_ARGS} --conf spark.executor.extraClassPath='/apps/eb/2020b/skylake/software/OpenBLAS/0.3.18-GCC-11.2.0/lib64/libopenblas.so'"
export LD_PRELOAD="${EBROOTOPENBLAS}/lib/libopenblas.so"

export HAIL_TMP_DIR="/well/lindgren/dpalmer"
export HAIL_SPARK_DIR="/well/lindgren/dpalmer/spark"

python pc_relate_minimal.py
# 06_0_pc_relate.py


hl.init(spark_conf={"spark.executor.extraClassPath": "/foo/bar/baz.so"})
hl.init(spark_conf={"SparkContext.setLogLevel": "DEBUG"})
hl.init(spark_conf={"spark.executor.extraJavaOptions": "-Dlog4j.debug=true"})