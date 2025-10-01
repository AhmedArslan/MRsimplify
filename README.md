# MRsimplify

TwosampleMR and MultivariableMR, perform all steps with simple command(s) without prior knowledge or having to go through lengthy boring protocols
### Internal steps in MRsimplify
 
 A:  _exposure data_: (i) read exposure data, (ii) perform SNP clumping and (iii) store data.
 
 B:  _outcome data_: (i) read outcome data, (ii) get proxy SNP(s)
 
 C:  _harmonise_ 
 
 D:  _MR_
 
 E:  _sensitivity tests_ (heterogeneity, pleiotropy, singlesnp, leaveoneout, MR-PRESSO)
 
 F:  _visualization_ (scatter plot, forest plot, leaveoneout and funnel plot)
 
 G:  compile all _results_ into a file.

## Step-1: installation..   
  **Install required R library:** 
   [TwoSampleMR](https://github.com/mrcieu/TwoSampleMR), [stringr](https://stringr.tidyverse.org), [tidyverse](https://www.tidyverse.org/packages/), [LDlinkR](https://cran.r-project.org/web/packages/LDlinkR/vignettes/LDlinkR.html), [ggplot2](https://ggplot2.tidyverse.org), [ieugwasr](https://mrcieu.github.io/ieugwasr/index.html), [dplyr](https://dplyr.tidyverse.org), [gwasvcf](https://github.com/MRCIEU/gwasvcf).
  
  **Download and install:** 
   * _R codes_ (MRsimplify.r) (before running add the path to (i) plink executable (line 9), (ii) local LD reference panel on line-20 and line-38). 
   * The _LD reference panel_ can be downloaded from [here](http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz) (currently supporting GRCh37/hg19 genome built).
   * The _LD reference panel_ contains information of 5 super-populations (EUR = European; EAS = East Asian; AMR = Admixed American; SAS = South Asian; AFR = African).

 **Download full gwas summary stat:**
   * From either [GWASCatalog](https://www.ebi.ac.uk/gwas/) or individual publications with necessary information: SNP, CHR, POS, A1 (effect_allele), A2 (other_allele), BETA, SE, Phenotype, Pval, EAF (effect_allele Freq), samplesize. (NOTE: all data must have GRCh37 coordinates for smooth processing and reliable results.)
   * In case genomic coordinates change required, [MungeSumstats](https://al-murphy.github.io/MungeSumstats/articles/MungeSumstats.html) can be used.

## Step-2: formate data..
  
  * Filter exposure data with above mentioned columns by pval (recommended: p<5e-08) whereas outcome data should be full length summary stats files _without_ pval threshold.
  * Note: To save time, (it is recommended to) include data of different exposure(s) into one file, however in all TwosampleMR subsequent steps each exposure-outcome MR is computed separately.

 ## Step-3: perform MR..
  
   Rscript --vanilla MRsimplify.r _exposure_ _outcome_

 ## Step-4: MR analysis results..

   a folder will be geneated with _outcome_ name containing all the results including _sensitivity tests_ plus all the _visualizations_
 
_________________________________________________________________________


## Caution: things to consider to perform a successful MR analysis...

1) exposure and outcome variants have same (i) _genomic positions_ and,(ii) _A1_ and _A2_ alleles

2) (i) LD reference penal must be upto date (1k_v3) and (ii) using same population as of used for the generation of exposure and outcome data

_________________________________________________________________________

 **additional readings:**
   https://mrcieu.github.io/TwoSampleMR/index.html
   
 **citation:** If you find repo useful please cite Exposure of Early Growth Traits Genetics and Childhood Disorders is Causally Associated with the Gallbladder Outcomes: A Mendelian Randomization study (2025) [https://www.researchsquare.com/article/rs-6234473/v1]

 **contact:** <ahmed.arslan@ulb.be> or leave comments in issues page. 
