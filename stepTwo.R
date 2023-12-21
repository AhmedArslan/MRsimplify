
library(TwoSampleMR)
library(stringr)
library(tidyverse)
library(LDlinkR)
library(ggplot2)
library(ieugwasr)
library(dplyr)
library(gwasvcf)
set_plink()

start_time <- Sys.time()

#step-2: READ outcome data:
  print("reading formatted exposure dataset... ")
  args = commandArgs(trailingOnly=TRUE)
  #remove the pattern
  file_name <- str_remove(string = basename(args[1]), pattern = '\\.txt')
  #read outcome
  print("reading outcome dataset... ")
  outCome_args <- read_outcome_data(filename=args[1],snp_col = "SNP", snps = NULL, sep="\t", beta_col = "beta", se_col = "se", effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "pval", chr_col = "chr", pos_col = "position", eaf_col = "eaf", samplesize="samplesize")
  #extract exposure SNPs present in outcome
  outcome_clump_semi <- semi_join(outCome_args, exposures_clump_local_preH, by = "SNP")
  #Exposure SNPs not present in outomce
  exposures_snps_not <- anti_join(exposures_clump_local_preH, outCome_args,  by = "SNP")
  print("clumping outcome SNPs... ")
  #identify local proxy snps
  LDproxy <- get_ld_proxies(rsid = unique(exposures_snps_not$SNP), bfile =  "/path/to/LD/EUR", tag_kb=10000, tag_r2=0.8, tag_nsnp=50000, threads = 4)
  #SNP proxy file cleaning before running next steps:
  LDproxy_reshuffle <- LDproxy[,c(3,7,1,2,8,4,10)]
  LDproxy_reshuffle$distance <- LDproxy$BP_A - LDproxy$BP_B
  alleles <- str_c(LDproxy$A1, '=' ,LDproxy$B1)
  alleles_B <- str_c(LDproxy$A2, '=' ,LDproxy$B2)
  alleles_C <- str_c(alleles, ',', alleles_B)
  LDproxy_reshuffle$correlated_alleles <- alleles_C
  LDproxy_reshuffle$dprime <- "1"
  colnames(LDproxy_reshuffle) <- c("snp", "proxy_snp", "chr", "pos", "alleles", "maf", "rsq", "distance", "correlated_alleles", "dprime")
  proxy_outcome <- left_join(
      LDproxy_reshuffle, outCome_args, by = c("proxy_snp" = "SNP")

    ) %>%
      separate(correlated_alleles, c("target_a1.outcome", "proxy_a1.outcome", "target_a2.outcome", "proxy_a2.outcome"), sep = ",|=") %>%
      filter(!is.na(chr.outcome)) %>%
      arrange(snp, -rsq, abs(distance)) %>%
      group_by(snp) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(
         proxy.outcome = TRUE,
         target_snp.outcome = snp,
         proxy_snp.outcome = proxy_snp, 
      ) %>% 
      mutate(
           new_effect_allele.outcome = case_when(
            proxy_a1.outcome == effect_allele.outcome & proxy_a2.outcome == other_allele.outcome ~ target_a1.outcome,
            proxy_a2.outcome == effect_allele.outcome & proxy_a1.outcome == other_allele.outcome ~ target_a2.outcome,
            TRUE ~ NA_character_
         ), 
          new_other_allele.outcome = case_when(
            proxy_a1.outcome == effect_allele.outcome & proxy_a2.outcome == other_allele.outcome ~ target_a2.outcome,
            proxy_a2.outcome == effect_allele.outcome & proxy_a1.outcome == other_allele.outcome ~ target_a1.outcome,
            TRUE ~ NA_character_
         ), 
         effect_allele.outcome = new_effect_allele.outcome, 
         other_allele.outcome = new_other_allele.outcome
      ) %>%
      dplyr::select(-proxy_snp, -chr, -pos, -alleles, -maf, -distance, -rsq, -dprime,  
             -new_effect_allele.outcome, -new_other_allele.outcome) %>%
      relocate(target_a1.outcome, proxy_a1.outcome, target_a2.outcome, proxy_a2.outcome, .after = proxy_snp.outcome) %>%
      rename(SNP = snp) %>%
      relocate(SNP, .after = samplesize.outcome)
    # Merge outcome and proxy outcomes
    outcome_dat_munge <- bind_rows(
      outcome_clump_semi, proxy_outcome
    ) %>% 
      arrange(chr.outcome, pos.outcome)

  rm(alleles)
  rm(alleles_C)
  rm(alleles_B)
  #step-3: Harmanise data
  print("harmonising outcome and exposure SNPs... ")
  harmonise_munge <- harmonise_data(exposure_dat = exposures_clump_local_preH, outcome_dat = outcome_dat_munge)
  print(dim(harmonise_munge))
  #step-4: MR
  print("performaing two sample MR... ") 
  mr_results <- mr(harmonise_munge)
  #step-5: sensitivity analyses:
  print("performaing sensitivity tests... ")
  mr_heterogeneity_munge <- mr_heterogeneity(harmonise_munge)
  mr_pleiotropy_munge <- mr_pleiotropy_test(harmonise_munge)
  res_single_munge <- mr_singlesnp(harmonise_munge, all_method = "mr_two_sample_ml")
  res_loo_munge <- mr_leaveoneout(harmonise_munge) 
  resPRESSO_munge <-run_mr_presso(harmonise_munge)
  #Combine all results
  print(paste("writing TwoSampleMR results for ", file_name, sep= ''))
  all_res_munge <-  combine_all_mrresults(mr_results, mr_heterogeneity_munge, mr_pleiotropy_munge, res_single_munge, ao_slc=FALSE)
  write.table(all_res_munge, file = paste('all_res_', file_name, '.txt', sep='') ,sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

  #ploting on local machine: (end dir have to change)
  print("creating scatter plots... ") 
  p1 <- mr_scatter_plot(mr_results, harmonise_munge)
  for (i in p1) { ggsave(i, file = paste(unique(c(i$data$exposure)),".p1.pdf", sep= ""), width = 9, height = 7)}
  print("creating funnel plots... ") 
  p4 <- mr_funnel_plot(res_single_munge)
  for (i in p4) { ggsave(i, file = paste(unique(c(i$data$exposure)),".p4.pdf", sep= ""), width = 9, height = 7)}

end_time <- Sys.time(); r <- end_time - start_time; r1 <- str_remove(string = r, pattern = 'Time difference of')
print(paste('TwoSampleMR for ', file_name, ' is finished in ', r1, sep= ''))
