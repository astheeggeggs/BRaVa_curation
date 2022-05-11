#!/bin/bash

PLINK_FILES=$1

for i in `seq 1 22`;
do
	plink --bfile ${PLINK_FILES}.chr${i} \
	--indep 50 5 2 --out ${PLINK_FILES}.pruned.chr${i}
done 

plink --bfile ${PLINK_FILES}.chrX --indep 50 5 2 --out ${PLINK_FILES}.pruned.chrX

cp ${PLINK_FILES}.chr1.prune.in ${PLINK_FILES}.keep.variant_list

for chr in `seq 2 22`;
do
	cat ${PLINK_FILES}.chr${chr}.prune.in >> ${PLINK_FILES}.keep.variant_list
done
