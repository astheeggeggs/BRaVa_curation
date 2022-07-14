library(bigsnpr)
library(data.table)
library(dplyr)
library(ggplot2)
library(randomForest)
library(ggrastr)

# Note, there is no step 4 in our WES 200k pipeline because we do not need to create a high quality plink file
# of variants - we already have the genotype data! These will allow us to perform PCA more accurately than 
# with the exome data alone. We also have the KING relatedness coefficients from UK Biobank, so don't need 
# thinned autosomal variants for IBD estimation either.

# Plotting functions:
source('utils/pretty_plotting.r')
# Thresholds and plotting file locations defined in r_options.r
source("utils/r_options.r")

# Function to sanity check the random matrix theory projections
project_onto_ref_PCs <- function(bed.ref, obj.bed, n_PCs = 20, n_PC_plot = 8,
  outdir = "/well/lindgren/dpalmer/ukb_get_EUR/", filename = "output", strand_flip=FALSE) {

  test <- bed_projectPCA(bed.ref, obj.bed, k = n_PCs, strand_flip = FALSE)#, ncores = nb_cores())

  plot(test$obj.svd.ref)
  ggsave(paste0(outdir, filename, "_scree.pdf"), width = 13, height = 7)
  plot(test$obj.svd, type = "loadings", loadings = 1:n_PCs, coeff = 0.4)
  ggsave(paste0(outdir, filename, "_loadings.pdf"), width = 13, height = 7)

  PC.ref <- predict(test$obj.svd.ref)
  proj1 <- test$simple_proj
  proj2 <- test$OADP_proj
  save("test", "PC.ref", "proj1", "proj2", file=paste0(outdir, filename, "_data.Rdata"))

  # shrinkage coefficients
  shrinkage <- unname(sapply(1:n_PCs, function(k) {
    MASS::rlm(proj2[, k] ~ proj1[, k] + 0, maxit=100)$coef
  }))
  print(round(shrinkage, 2))

  ind <- seq(1, min(20e3, nrow(proj1)))

  plot_grid(plotlist = lapply(1:floor(n_PC_plot/2), function(k) {
    k1 <- 2 * k - 1
    k2 <- 2 * k
    plot_grid(
      qplot(PC.ref[, k1], PC.ref[, k2], size = I(2)) +
        geom_point(aes(proj1[ind, k1], proj1[ind, k2]), color = "red") +
        theme_bigstatsr(0.5) +
        labs(x = paste0("PC", k1), y = paste0("PC", k2)),
      qplot(PC.ref[, k1], PC.ref[, k2], size = I(2)) +
        geom_point(aes(proj2[ind, k1], proj2[ind, k2]), color = "blue") +
        theme_bigstatsr(0.5) +
        labs(x = paste0("PC", k1), y = paste0("PC", k2)),
      scale = 0.95
    )
  }), nrow = floor(sqrt(n_PCs)))
  ggsave(file=paste0(outdir, filename, "_shrinkage.pdf"), width=15, height=15)
  return(round(shrinkage, 2))

}

# Generate the fam file for all people, so that that it can be read by plink.
ukb_fam_lindgren <- "/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"
ukb_fam <- "/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb11867_cal_chr1_v2_s488363_for_plink.fam"
fwrite(fread(ukb_fam) %>% mutate(V6=-9), file=ukb_fam, sep=' ', col.names=FALSE)

# First, generate cleaned up genotype data
write(sapply(1:22, function(chr) {
  c(paste0("/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_cal_chr", chr, "_v2.bed"),
    paste0("/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_snp_chr", chr, "_v2.bim"),
    ukb_fam)
}), '/well/lindgren/dpalmer/ukb_get_EUR/plink_ukb_filepaths.txt', ncolumns = 3)

# If this dies, it's likely due to memory. Requires at least 10 slots on qe.
snp_plinkQC(
  plink.path = '/well/lindgren/dpalmer/plink',
  prefix.in = '/well/lindgren/dpalmer/ukb_get_EUR/plink_ukb_filepaths.txt',
  file.type = "--merge-list",
  prefix.out = "/well/lindgren/UKBIOBNK/dpalmer/ukb_get_EUR/data/ukb_autosomes_combined",
  geno = 0.01,
  autosome.only = TRUE
)

obj.bed <- bed('/well/lindgren/UKBIOBNK/dpalmer/ukb_get_EUR/data/ukb_autosomes_combined.bed')
bed.ref <- bed(download_1000G("/well/lindgren/dpalmer/ukb_get_EUR/data"))

outdir <- '/well/lindgren/dpalmer/ukb_get_EUR/'
filename <- 'ukb_projected_to_1kg_PCs'
n_PCs <- 10
save_figures <- TRUE
perform_plotting <- TRUE
creating_new_EUR_def <- TRUE

project_onto_ref_PCs(bed.ref, obj.bed, filename=filename, outdir=outdir, n_PCs=n_PCs)

# Grab 1000G information containing populations
# download.file("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped",
#               destfile = "/well/lindgren/UKBIOBANK/dpalmer/1kg_for_EUR_assign/20130606_g1k.ped")
dt_pop <- fread("/well/lindgren/UKBIOBANK/dpalmer/1kg_for_EUR_assign/20130606_g1k.ped") %>% 
  mutate(sample.ID = `Individual ID`, population=Population) %>% 
  select(sample.ID, population) %>%
  # Determine the superpopulations
  mutate(super.population = case_when(
    population %in% c("CHB", "JPT", "CHS", "CDX", "KHV") ~ "EAS",
    population %in% c("CEU", "TSI", "FIN", "GBR", "IBS") ~ "EUR",
    population %in% c("YRI", "LWK", "MAG", "MSL", "ESN", "ASW", "ACB", "GWD") ~ "AFR",
    population %in% c("AMR", "PUR", "CLM", "PEL", "MXL") ~ "AMR",
    population %in% c("GIH", "PJL","BEB", "STU", "ITU") ~ "SAS",
    TRUE ~ "other"
  )
)

load(paste0(outdir, filename, "_data.Rdata"))
PC.ref <- data.table(PC.ref)
names(PC.ref) <- paste0("PC", seq(1,n_PCs))

dt <- cbind(data.table(bed.ref$fam), PC.ref)
setkey(dt, "sample.ID")
setkey(dt_pop, "sample.ID")
dt_1kg <- merge(dt, dt_pop)

proj2 <- data.table(proj2)
names(proj2) <- paste0("PC", seq(1,n_PCs))
dt_ukb <- cbind(data.table(obj.bed$fam), proj2)

# Now, need to grab the population labels to train the random forest.
dt_train = dt_1kg %>%
  select(c(super.population, population, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))

dt_predict = dt_ukb %>%
  select(c(sample.ID, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))

PCs_to_use <- paste0('PC', seq(1,4))

# Determine a classifier.
set.seed(160487)
T_RF <- 0.99
rf <- randomForest(x=dt_train %>% select(PCs_to_use), y=as.factor(as.character(dt_train$super.population)), ntree=10000)
rf_probs <- predict(rf, dt_predict %>% select(PCs_to_use), type='prob')

rf_probs_plot <- melt(data.table(rf_probs), measure.vars=c("AFR", "AMR", "EAS", "EUR", "SAS"),
  variable.name="population", value.name="prob")

p <- ggplot(rf_probs_plot, aes(x=prob, color=factor(population))) + geom_density(stat="density") + 
  theme_classic() + labs(color='1000G label', x='Random Forest Assignment Probability', y="Density") +
  geom_vline(xintercept=T_RF, linetype='dashed')
ggsave(paste0(PLOTS, '05_random_forest_assignment_density', '.pdf'), p, width=160, height=90, units='mm')
ggsave(paste0(PLOTS, '05_random_forest_assignment_density', '.jpg'), p, width=160, height=90, units='mm', dpi=500)
p <- ggplot(rf_probs_plot, aes(x=prob, color=factor(population))) + geom_density(stat="density") + 
  theme_classic() + labs(color='1000G label', x='Random Forest Assignment Probability', y="Density") +
  xlim((T_RF - 0.05), NA) + geom_vline(xintercept=T_RF, linetype='dashed')
ggsave(paste0(PLOTS, '05_random_forest_assignment_density_zoom', '.pdf'), p, width=160, height=90, units='mm')
ggsave(paste0(PLOTS, '05_random_forest_assignment_density_zoom', '.jpg'), p, width=160, height=90, units='mm', dpi=500)
p <- ggplot(rf_probs_plot, aes(x=prob, color=factor(population))) + geom_density(stat="ecdf") + 
  theme_classic() + labs(color='1000G label', x='Cumulative Proportion', y="Cumulative Density") +
  geom_vline(xintercept=T_RF, linetype='dashed')
ggsave(paste0(PLOTS, '05_random_forest_assignment_cum_density', '.pdf'), p, width=160, height=90, units='mm')
ggsave(paste0(PLOTS, '05_random_forest_assignment_cum_density', '.jpg'), p, width=160, height=90, units='mm', dpi=500)
p <- ggplot(rf_probs_plot, aes(x=prob, color=factor(population))) + geom_density(stat="ecdf") + 
  theme_classic() + labs(color='1000G label', x='Cumulative Proportion', y="Cumulative Density") +
  xlim((T_RF - 0.49), NA)  + geom_vline(xintercept=T_RF, linetype='dashed')
  # coord_cartesian(xlim = c((T_RF - 0.05), NA)) + geom_vline(xintercept=T_RF, linetype='dashed')
ggsave(paste0(PLOTS, '05_random_forest_assignment_cum_density_zoom', '.pdf'), p, width=160, height=90, units='mm')
ggsave(paste0(PLOTS, '05_random_forest_assignment_cum_density_zoom', '.jpg'), p, width=160, height=90, units='mm', dpi=500)

check_thres <- function(row, threshold) {
  return(!any(row > threshold))
}

unsure <- apply(rf_probs, 1, check_thres, T_RF)
classification <- as.character(predict(rf, dt_predict %>% select(PCs_to_use)))
dt_predict$classification_loose <- as.factor(classification)
classification[unsure] <- 'unsure'
dt_predict$classification_strict <- as.factor(classification)

dt_predict <- dt_predict %>% mutate(sample.ID = as.character(sample.ID)) %>% select(sample.ID, classification_strict, classification_loose, starts_with("PC"))
dt_predict <- data.table(dt_predict)
setkeyv(dt_predict, c("sample.ID", paste0("PC", seq(1,10))))
setkeyv(dt_1kg, c("sample.ID", paste0("PC", seq(1,10))))
dt_classify <- merge(dt_predict, dt_1kg, all=TRUE) %>% select(sample.ID, classification_strict, classification_loose, starts_with("PC"), super.population, population)

# Write the result to file.
fwrite(dt_classify, file = "/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/superpopulation_labels.tsv", sep='\t')

PCs <- c(1,3,5)

if (perform_plotting)
{
    dt_plot <- dt_classify %>% filter(!is.na(classification_loose))
    levels(dt_plot$classification_loose) <- c(levels(dt_plot$classification_loose) ,"1000 genomes")
    dt_plot$classification_loose[is.na(dt_plot$classification_loose)] <- "1000 genomes"

    for (i in PCs)
    {
        aes <- aes_string(x=paste0('PC', i), y=paste0('PC', i+1), color='classification_loose')
        p <- create_pretty_scatter(dt_plot, aes,
          save_figure=save_figures, file=paste0(PLOTS,'05_PC',i,'_PC',i+1,'_classify_loose'), n_x_ticks=5,
          x_label=paste0('Principal Component ',i),
          y_label=paste0('Principal Component ', i+1)
          )
        
        p <- p + geom_point(data=dt_1kg, mapping=aes_string(x=paste0('PC',i), y=paste0('PC',i+1), shape='super.population'),
          inherit.aes=FALSE, show.legend=FALSE)
        # ggsave(file=paste0(PLOTS,'05_PC',i,'_PC',i+1,'_classify_loose_1kg_labelled.pdf'), width=160, height=90, units='mm')
        ggsave(file=paste0(PLOTS,'05_PC',i,'_PC',i+1,'_classify_loose_1kg_labelled.jpg'), width=160, height=90, units='mm', dpi=500)

        aes <- aes_string(x=paste0('PC',i), y=paste0('PC',i+1), color='classification_strict')
        p <- create_pretty_scatter(dt_plot, aes,
          save_figure=FALSE, file=paste0(PLOTS,'00_PC',i,'_PC',i+1,'_classify_strict'), n_x_ticks=5,
          x_label=paste0('Principal Component ',i),
          y_label=paste0('Principal Component ', i+1)
          )
        p <- p + geom_point(data=dt_1kg, mapping=aes_string(x=paste0('PC',i), y=paste0('PC',i+1), shape='super.population'),
          inherit.aes=FALSE, show.legend=FALSE)
        # ggsave(file=paste0(PLOTS,'05_PC',i,'_PC',i+1,'_classify_strict_1kg_labelled.pdf'), width=160, height=90, units='mm')
        ggsave(file=paste0(PLOTS,'05_PC',i,'_PC',i+1,'_classify_strict_1kg_labelled.jpg'), width=160, height=90, units='mm', dpi=500)

        p <- create_pretty_scatter(dt_plot %>% 
          filter(classification_strict != "unsure") %>% 
          mutate(classification_strict = as.factor(as.character(classification_strict))),
          aes, save_figure=FALSE, n_x_ticks=5,
          x_label=paste0('Principal Component ',i),
          y_label=paste0('Principal Component ', i+1)
          )
        # ggsave(file=paste0(PLOTS,'05_PC',i,'_PC',i+1,'_classify_strict_1kg_labelled.pdf'), width=160, height=90, units='mm')
        ggsave(file=paste0(PLOTS,'05_PC',i,'_PC',i+1,'_classify_strict_1kg_labelled.jpg'), width=160, height=90, units='mm', dpi=500)

    }
}
