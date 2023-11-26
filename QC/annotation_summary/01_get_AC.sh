# author: Duncan Palmer

in_dir="/mnt/project/Barney/wes/sample_filtered"
out_dir="/Duncan/counts"

for anc in AFR AMR EAS EUR SAS; do
   for CHR in {{1..22},X}; do
       in_file="${in_dir}/ukb_wes_450k.qced.chr${CHR}"
       dx run app-swiss-army-knife -icmd="
         awk '{print 0, \$1}' /mnt/project/brava/inputs/ancestry_sample_ids/qced_${anc}_sample_IDs.txt > tmp.txt; \
         plink --keep tmp.txt \
         --freqx gz --bfile ${in_file}; \
         mv plink.frqx.gz plink.frqx.chr${CHR}.${anc}.gz; \
         rm tmp.txt
         " \
     --instance-type mem2_ssd1_v2_x4 \
     --folder=".${out_dir}" \
     --name AC_chr${CHR} \
     --priority normal -y
   done
done
