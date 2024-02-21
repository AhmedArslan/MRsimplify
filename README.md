# MRsimplify

TwosampleMR and MultivariableMR, perform all steps with simple command(s) without prior knowledge or having to go through lengthy boring protocols
### Internal steps in MRsimplify
 
 A:  _exposure data_: (i) read exposure data, (ii) perform SNP clumping and (iii) store data.
 
 B:  _outcome data_: (i) read outcome data, (ii) get proxy SNP(s)
 
 C:  _harmonise_ 
 
 D:  _MR_
 
 E:  _sensitivity tests_ (heterogeneity, pleiotropy, singlesnp, leaveoneout, MR-PRESSO)
 
 F:  _visualization_ (scatter plots and funnel plots)
 
 G:  compile all _results_ into a file.

## Step-1: installation..   
  **Install required R library:** 
   [TwoSampleMR](https://github.com/mrcieu/TwoSampleMR), [stringr](https://stringr.tidyverse.org), [tidyverse](https://www.tidyverse.org/packages/), [LDlinkR](https://cran.r-project.org/web/packages/LDlinkR/vignettes/LDlinkR.html), [ggplot2](https://ggplot2.tidyverse.org), [ieugwasr](https://mrcieu.github.io/ieugwasr/index.html), [dplyr](https://dplyr.tidyverse.org), [gwasvcf](https://github.com/MRCIEU/gwasvcf).
  
  **Download and install:** 
   * _R codes_ (stepOne.r and  stepTwo.r) and before running add the path to local LD reference panel on line-12 and line-28, respectively. 
   * The _LD reference panel_ can be downloaded from [here](http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz) (currently supporting GRCh37/hg19 genome built).
   * The _LD reference panel_ contains information of 5 super-populations (EUR = European; EAS = East Asian; AMR = Admixed American; SAS = South Asian; AFR = African).

 **Download full gwas summary stat:**
   * From either [GWASCatalog](https://www.ebi.ac.uk/gwas/) or individual publications with necessary information: SNP, CHR, POS, A1 (effect_allele), A2 (other_allele), BETA, SE, Phenotype, Pval, EAF (effect_allele Freq), samplesize. (NOTE: all data must have GRCh37 coordinates for smooth processing and reliable results.)
   * In case genomic coordinates change required, [MungeSumstats](https://neurogenomics.github.io/MungeSumstats/articles/MungeSumstats.html) can be used.

## Step-2: formate data..
  
  * Filter exposure data with above mentioned columns by pval (<5e-08) whereas outcome data should be full length summary stats files _without_ pval threshold.
  * Note: To save time, (it is recommended to) include data of different exposure(s) into one file, however in all TwosampleMR subsequent steps each exposure-outcome MR is computed separately.

## Step-3: read exposure data.. 

(_internal step: A_)
  
  Rscript --vanilla stepOne.r _exposure file_

## Step-4: read outcome data and perform TwosampleMR.. 

(_internal steps: B - G_)
 
  Rscript --vanilla stepTwo.r _outcome file_

_________________________________________________________________________

 **additional readings:**
   https://mrcieu.github.io/TwoSampleMR/index.html
   
 **citation:** If you find repo useful please cite the link while manuscript is in preparation. 

 **contact:** <ahmed.arslan@ulb.be> or leave comments in issues page. 
