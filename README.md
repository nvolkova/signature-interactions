# signature-interactions
Codes and data for the manuscript "Mutational signatures are jointly shaped by DNA damage and repair"
(to be submitted in June 2019)

Human analysis part requires a matrix with TCGA variants, and a metadata matrix with sample names, cancer types, and FPKM expression values for DNA repair genes.

Session information for the analyses:

```
R version 3.5.1 (2018-07-02)

Matrix products: default

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base   

other attached packages:
 [1] greta_0.2.3                 reshape2_1.4.3              ggplot2_3.1.1               VariantAnnotation_1.30.0   
 [5] Rsamtools_2.0.0             Biostrings_2.52.0           XVector_0.24.0              SummarizedExperiment_1.14.0
 [9] DelayedArray_0.10.0         BiocParallel_1.18.0         matrixStats_0.54.0          Biobase_2.44.0             
[13] GenomicRanges_1.36.0        GenomeInfoDb_1.20.0         IRanges_2.18.0              S4Vectors_0.22.0           
[17] BiocGenerics_0.30.0         dndscv_0.0.1.0              bayesplot_1.6.0             RColorBrewer_1.1-2

loaded via a namespace (and not attached):
[1] httr_1.4.0               bit64_0.9-7              jsonlite_1.6             assertthat_0.2.1         blob_1.1.1              
[6] BSgenome_1.52.0          GenomeInfoDbData_1.2.1   progress_1.2.1           pillar_1.4.0             RSQLite_2.1.1           
[11] lattice_0.20-38          glue_1.3.1               reticulate_1.12          digest_0.6.18            colorspace_1.4-1        
[16] Matrix_1.2-17            plyr_1.8.4               XML_3.98-1.19            pkgconfig_2.0.2          biomaRt_2.40.0          
[21] zlibbioc_1.30.0          purrr_0.3.2              scales_1.0.0             whisker_0.3-2            tibble_2.1.1            
[26] withr_2.1.2              GenomicFeatures_1.36.0   lazyeval_0.2.2           magrittr_1.5             crayon_1.3.4            
[31] memoise_1.1.0            MASS_7.3-51.4            tools_3.6.1              prettyunits_1.0.2        hms_0.4.2               
[36] stringr_1.4.0            munsell_0.5.0            AnnotationDbi_1.46.0     ade4_1.7-13              compiler_3.6.1          
[41] rlang_0.3.4              RCurl_1.95-4.12          ggridges_0.5.1           rstudioapi_0.10          bitops_1.0-6            
[46] base64enc_0.1-3          labeling_0.3             gtable_0.3.0             DBI_1.0.0                R6_2.4.0                
[51] GenomicAlignments_1.20.0 tfruns_1.4               dplyr_0.8.0.1            tensorflow_1.13.1        rtracklayer_1.40.2      
[56] seqinr_3.4-5             bit_1.1-14               stringi_1.4.3            Rcpp_1.0.1               tidyselect_0.2.5        
[61] coda_0.19-2 
```
