#!/bin/bash

# This is https://github.com/BRaVa-genetics/BRaVa_curation/blob/main/QC/annotation_summary/02_merged_annots_and_AC.r
rscript_remote="Duncan/counts/scripts/02_merged_annots_and_AC.r"

in_dir="/mnt/project/Duncan/long_annotations"
out_dir="/Duncan/annotation_counts"
dx mkdir -p ${out_dir}

for anc in AFR AMR EAS EUR SAS; do
   for CHR in {{1..22},X}; do
      out_prefix="ukb_wes_450k.${anc}.chr${CHR}"
      AC_path=/mnt/project/Duncan/counts/plink.frqx.chr${CHR}.${anc}.gz
      vep_processed_long_path="${in_dir}/ukb_wes_450k.july.qced.brava_common_rare.v7.chr${CHR}.saige.txt.long.csv.gz"
      dx run app-swiss-army-knife \
        -iimage_file="/docker/rsuite.tar.gz"\
        -icmd="
           Rscript /mnt/project/${rscript_remote} \
            --AC_path ${AC_path} \
            --vep_spliceAI_processed ${vep_processed_long_path} \
            --out ${out_prefix} \
            --spliceAI_bins
          "\
        --instance-type mem1_ssd1_v2_x8 \
        --folder=".${out_dir}" \
        --priority normal \
        --name annotations_counts -y
  done
done
