RHOGE
================
Nick Mancuso
2018-02-27

RHOGE is an R package that computes estimates the genome-wide genetic correlation between two complex traits (diseases) as a function of predicted gene expression effect on trait (\\rho\_{ge}). Given output from two transcriptome-wide association studies, RHOGE estimates the mediating effect of predicted gene expression and estimates the correlation of effect sizes. This approach is extended to a bi-directional regression that provides putative causal directions between traits with non-zero \\rho\_{ge}.

This approach is described in:

> [Integrating Gene Expression with Summary Association Statistics to Identify Genes Associated with 30 Complex Traits](https://doi.org/10.1016/j.ajhg.2017.01.031)
> Nicholas Mancuso, Huwenbo Shi, Pagé Goddard, Gleb Kichaev, Alexander Gusev, Bogdan Pasaniuc.
> American Journal of Human Genetics. 2017.

Installation
============

Bioconductor + devtools is the most straightforward way to install RHOGE. To do this open an R terminal and enter

``` r
source("http://bioconductor.org/biocLite.R")
biocLite("devtools")    # only if devtools not yet installed
biocLite("bogdanlab/RHOGE")
```

Example
=======

The following example computes \\rho\_{ge} between BMI and triglycerides, as well as putative causal directions.

``` r
library(RHOGE)

# example BMI TWAS results from FUSION
data(bmi)
head(bmi)
```

    ## # A tibble: 6 x 22
    ##   FILE   ID      CHR     P0     P1    HSQ BEST.GWAS.ID BEST.GWAS.Z EQTL.ID
    ##   <chr>  <chr> <int>  <dbl>  <dbl>  <dbl> <chr>              <dbl> <chr>  
    ## 1 /u/ho… LOC6…     1 7.63e5 7.95e5 0.0979 rs17160824          1.96 rs3094…
    ## 2 /u/ho… AGRN      1 9.56e5 9.91e5 0.0655 rs17160824          1.96 rs9442…
    ## 3 /u/ho… C1or…     1 1.02e6 1.05e6 0.0917 rs17160824          1.96 rs3766…
    ## 4 /u/ho… SCNN…     1 1.22e6 1.23e6 0.170  rs9660180           4.89 rs1126…
    ## 5 /u/ho… MXRA8     1 1.29e6 1.30e6 0.0503 rs9660180           4.89 rs1739…
    ## 6 /u/ho… LOC1…     1 1.33e6 1.34e6 0.0648 rs9660180           4.89 rs2649…
    ## # ... with 13 more variables: EQTL.R2 <dbl>, EQTL.Z <dbl>,
    ## #   EQTL.GWAS.Z <dbl>, NSNP <int>, NWGT <int>, MODEL <chr>,
    ## #   MODELCV.R2 <dbl>, MODELCV.PV <dbl>, TWAS.Z <dbl>, TWAS.P <dbl>,
    ## #   PERM.PV <dbl>, PERM.N <int>, PERM.ANL_PV <dbl>

``` r
# example triglyceride TWAS results from FUSION
data(tg)
head(tg)
```

    ## # A tibble: 6 x 22
    ##   FILE   ID      CHR     P0     P1    HSQ BEST.GWAS.ID BEST.GWAS.Z EQTL.ID
    ##   <chr>  <chr> <int>  <dbl>  <dbl>  <dbl> <chr>              <dbl> <chr>  
    ## 1 /u/ho… LOC6…     1 7.63e5 7.95e5 0.0979 rs4314833           2.32 rs3094…
    ## 2 /u/ho… AGRN      1 9.56e5 9.91e5 0.0655 rs4314833           2.31 rs9442…
    ## 3 /u/ho… C1or…     1 1.02e6 1.05e6 0.0917 rs11552172         -2.36 rs3766…
    ## 4 /u/ho… SCNN…     1 1.22e6 1.23e6 0.170  rs6604981           2.44 rs1126…
    ## 5 /u/ho… MXRA8     1 1.29e6 1.30e6 0.0503 rs13303010          2.60 rs1739…
    ## 6 /u/ho… LOC1…     1 1.33e6 1.34e6 0.0648 rs13303010          2.89 rs2649…
    ## # ... with 13 more variables: EQTL.R2 <dbl>, EQTL.Z <dbl>,
    ## #   EQTL.GWAS.Z <dbl>, NSNP <int>, NWGT <int>, MODEL <chr>,
    ## #   MODELCV.R2 <dbl>, MODELCV.PV <dbl>, TWAS.Z <dbl>, TWAS.P <dbl>,
    ## #   PERM.PV <dbl>, PERM.N <int>, PERM.ANL_PV <dbl>

``` r
# Estimate rho_ge genome-wide for BMI and Triglyerides and approximate sample sizes
ge_cor_res <- rhoge.gw(bmi, tg, 14000, 91000)
head(ge_cor_res)
```

    ##       RHOGE         SE    TSTAT  DF            P
    ## 1 0.2168836 0.05869152 3.695313 461 0.0002459941

``` r
bidir_res <- rhoge.bd(bmi, tg, 14000, 91000, p1 = 0.05 / nrow(bmi), p2 = 0.05 / nrow(tg))
head(bidir_res)
```

    ##     ESTIMATE         SE      TSTAT       DF            P             TEST
    ## 1  0.5981401 0.09416694  6.3519115 35.00000 2.662674e-07 Trait1 -> Trait2
    ## 2 -0.1068047 0.20526312 -0.5203307 35.00000 6.061092e-01 Trait2 -> Trait1
    ## 3         NA         NA  3.1215369 50.51057 2.973945e-03             DIFF

Currently, only [FUSION](https://github.com/gusevlab/fusion_twas) style output is supported.
