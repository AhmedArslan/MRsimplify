# MRsimplify
TwosampleMR and MultivariableMR perform all steps with simple command(s) without prior knowledge or having to go through lengthy boring protocols

## Step-1:
  ** Download full gwas summary stat: **
  From either GWASCatalog or individual publications with necessary information: SNP, CHR, POS, A1 (effect_allele), A2 (other_allele), BETA, SE, Phenotype, Pval, EAF (effect_allele Freq), Samplesize. 
  
  ** Install required R library: ** 
  TwoSampleMR, stringr, tidyverse, LDlinkR, ggplot2, ieugwasr, dplyr, gwasvcf

## Step-2: 
  
  Filter Exposure data of above columns with pval (<5e-08) whereas Outcome data should be full length summary stats files without pval threshold

## Step-3: read exposure data:
  
  Rscript --vanilla stepOne.r < exposure file >

## Step-4: read outcome data and perform MR:
 
  Rscript --vanilla stepTwo.r < outcome file >


## Explanation:

Step-3:
  (i) read exposure data, (ii) perform SNP clumping and (iii) prepare data for step-4

Step-4:
  (i) read outcome data, (ii) get proxy SNP(s); (iii) hormonise SNPs, (iv) perform twosampleMR, (v) sensitivity tests (heterogeneity, pleiotropy, singlesnp, leaveoneout, presso), (vi) write output files containing all the results. (vii) Generate visulations in terms og scatter plots and funnel plots. 

  contact: <ahmed.arslan@ulb.be> or leave comment in issues page. 
