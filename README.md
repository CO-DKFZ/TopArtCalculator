# TopArtCalculator

The TOP-ART score is a computational biomarker to assess the level of 
homologous-recombination deficiency (HRD) in tumor cells (cit. TOP-ART manuscript).
The current implementation as an R-package can be used to calculate the TOP-ART score
from bulk NGS whole-genome sequencing (WGS) or whole-exome sequencing (WES) data. 

## Installation

```r
library(devtools) # required to install R-package from github
install_github("https://github.com/CO-DKFZ/TopArtCalculator") 
```

## Prerequisites 

Before TOP-ART calculator can be used, the following upstream steps need to be performed:

* somatic copynumber calling (somatic copynumber abberations, sCNAs)
  - Annotation of **total copynumber of segments**
  - Info whether segment has **LOH**
* somatic SNV + Indel calling (somatic small variants)
  - Annotation of **allele-frequency**, **alternative base**, **reference base**
  - Annotation of **gene name** the variant is located in (in line with gene model used for all annotations of input data)
* germline SNV + Indel calling (germline small variants)
  - Annotation of **allele-frequency**, **alternative base**, **reference base**
  - Annotation of **gene name** the variant is located in (in line with gene model used for all annotations of input data)
  - Information about **Pathogenicity** in ACMG classification (1-5; benign, likely benign, uncertain significance, likely pathogenic, pathogenic) (https://github.com/NagaComBio/CharGer)
* calling of mutational signature 3 ([Alexandrov et al.](https://www.sciencedirect.com/science/article/pii/S0959437X13001639?via%3Dihub))
  - Info about confidence interval (Recommendation: [YAPSA](https://bioconductor.org/packages/YAPSA.html) package ([Hübschmann et al.](https://onlinelibrary.wiley.com/doi/10.1002/gcc.22918))
  
## Example

The tool mainly relies on the [ZygosityPredictor](https://bioconductor.org/packages/ZygosityPredictor.html) R-package ([Rheinnecker et al.](https://academic.oup.com/bioinformaticsadvances/article/4/1/vbae017/7601458)).
Some inputs are therefore inherited and require the same format as if used
for ZygosityPredictor itself. The following example uses the minimum required input to run.
If applicable, additional input files can be provided (RNA-seq, haploblocks, SNPs, see [ZygosityPredictor vignette](https://bioconductor.org/packages/release/bioc/vignettes/ZygosityPredictor/inst/doc/Usage.html)) to
increase chances of successful haplotype phasing of variants.


```r
   
 ## we use the example dataset from the Bioconductor package ZygosiytPredictor

 library(ZygosityPredictor)
 
 bamfile <- system.file("extdata", "ZP_example.bam", 
                        package = "ZygosityPredictor")
                        
 ## load ZygosityPredictor variant datasets
 data("GR_GERM_SMALL_VARS")
 data("GR_SCNA")
 data("GR_SOM_SMALL_VARS")
 data("GR_GENE_MODEL")
 
 ## create mutational signatures table as it is produced by bioconductor package
 ## YAPSA
 yapsaTable <- data.frame(
   sig=c("AC1","AC3","AC8"),
   exposure=c(4566, 2123, 598),
   lower=c(3899, 1902, 430),
   upper=c(5698, 2455, 789)
 )
 
 library(TopArtCalc)
 ## calculate TOP-ART score
 topart_run <- calculate_TOP_ART_score(
   seqMethod="WGS", ## can be WGS or WES
   HRD=10, ## count from sCNA calling
   LST=4, ## count from sCNA calling
   purity=0.9, ## 1 == 100 % tumor cells in bulk sample
   ploidy=2, ## overall ploidy of the tumor cells
   sex="female",
   ## files
   bamDna=bamfile, ## path to sequence alignment in .bam.gz format
   ## granges
   somCna=GR_SCNA,   # meta cols required: tcn, cna_type
   somSmallVars=GR_SOM_SMALL_VARS,   # meta cols required: gene, af, ref, alt
   germSmallVars=GR_GERM_SMALL_VARS,  # meta cols required: gene, af, ref, alt, ACMG_class
   geneModel=GR_GENE_MODEL,     # meta cols required: gene
   yapsaTable=yapsaTable
   )
   
   ## the output contains an overview tibble and an output message:
   cat(topart_run$output_message)
   
```




