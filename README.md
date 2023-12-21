# MRsimplify

TwosampleMR and MultivariableMR perform all steps with simple command(s) without prior knowledge or having to go through lengthy boring protocols
### Internal steps in MRsimplify
 exposure data: read; format data; SNP clumping
  outcome data read; proxy SNP 
   Harmonize 
    MR
     Sensitivity tests (heterogeneity, pleiotropy, singlesnp, leaveoneout, presso)
      Visualization (scatter plots and funnel plots)

## Step-1: installtion..   
  **Install required R library:** 
   [TwoSampleMR](https://github.com/mrcieu/TwoSampleMR), [stringr](https://stringr.tidyverse.org), [tidyverse](https://www.tidyverse.org/packages/), [LDlinkR](https://cran.r-project.org/web/packages/LDlinkR/vignettes/LDlinkR.html), [ggplot2](https://ggplot2.tidyverse.org), [ieugwasr](https://mrcieu.github.io/ieugwasr/index.html), [dplyr](https://dplyr.tidyverse.org), [gwasvcf](https://github.com/MRCIEU/gwasvcf).
  
  **Download and install:** 
   _R codes_ (stepOne.r and  stepTwo.r) and before running add the path to local LD reference panel on line-12 and line-28, respectively. 
   _The LD reference panel_ contains information on 5 super-populations (EUR = European; EAS = East Asian; AMR = Admixed American; SAS = South Asian; AFR = African). The LD reference panel can be downloaded from: [mrcieu](http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz) (correctly supporting only GRCh37/hg19 reference built).

 **Download full gwas summary stat:**
   From either [GWASCatalog](https://www.ebi.ac.uk/gwas/) or individual publications with necessary information: SNP, CHR, POS, A1 (effect_allele), A2 (other_allele), BETA, SE, Phenotype, Pval, EAF (effect_allele Freq), samplesize. (NOTE: all data must have GRCh37 corrdinates for smooth processing and reliable results)

## Step-2: formate data..
  
  Filter exposure data with mentioned above columns by pval (<5e-08) whereas outcome data should be full length summary stats files _without_ pval threshold.

## Step-3: read exposure data..
  
  Rscript --vanilla stepOne.r _exposure file_

## Step-4: read outcome data and perform twosampleMR..
 
  Rscript --vanilla stepTwo.r _outcome file_


### Explanation:

**Step-3:**
  (i) read exposure data, (ii) perform SNP clumping and (iii) prepare data for step-4

**Step-4:**
  (i) read outcome data, (ii) get proxy SNP(s); (iii) hormonise SNPs, (iv) perform twosampleMR, (v) sensitivity tests (heterogeneity, pleiotropy, singlesnp, leaveoneout, presso), (vi) write output files containing all the results. (vii) Generate visulations in terms of scatter plots and funnel plots. 

**contact:** <ahmed.arslan@ulb.be> or leave comments in issues page. 
