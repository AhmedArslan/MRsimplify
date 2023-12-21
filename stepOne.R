library(TwoSampleMR)
library(ieugwasr)
set_plink()

#Step-1: READ Exposure data:
args = commandArgs(trailingOnly=TRUE)
exposures <- read.table(args[1], header = TRUE, sep = "\t")
SE <- TwoSampleMR:::get_se(exposures$BETA, exposures$P)
exposures$se <- SE
colnames(exposures) <- c("SNP", "chr", "position", "effect_allele", "other_allele", "beta", "Phenotype", "pval", "se", "eaf", "samplesize")
exposures_frmt <- format_data(exposures, type = "exposure")
exposures_clump_local <- ld_clump(dplyr::tibble(rsid=exposures$SNP, pval=exposures$pval, id=exposures$Phenotype), plink_bin = genetics.binaRies::get_plink_binary(), bfile = "/path/to/LD/EUR")
exposures_clump_local_preH$pval <- NULL
exposures_clump_local_preH$id <- NULL
rm(exposures_clump_local)
#keep count of exposures snps
exposures_snps_count <- exposures_clump_local_preH %>% count(exposure); exposures_snps_count <- arrange(exposures_snps_count, desc(n))
