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
[17] BiocGenerics_0.30.0        

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.1               lattice_0.20-38          prettyunits_1.0.2        assertthat_0.2.1         digest_0.6.18           
 [6] R6_2.4.0                 plyr_1.8.4               coda_0.19-2              RSQLite_2.1.1            httr_1.4.0              
[11] pillar_1.4.0             tfruns_1.4               zlibbioc_1.30.0          rlang_0.3.4              GenomicFeatures_1.36.0  
[16] progress_1.2.1           lazyeval_0.2.2           rstudioapi_0.10          whisker_0.3-2            blob_1.1.1              
[21] Matrix_1.2-17            reticulate_1.12          stringr_1.4.0            RCurl_1.95-4.12          bit_1.1-14              
[26] biomaRt_2.40.0           munsell_0.5.0            compiler_3.6.1           rtracklayer_1.40.2       base64enc_0.1-3         
[31] pkgconfig_2.0.2          tensorflow_1.13.1        tidyselect_0.2.5         tibble_2.1.1             GenomeInfoDbData_1.2.1  
[36] XML_3.98-1.19            crayon_1.3.4             dplyr_0.8.0.1            withr_2.1.2              GenomicAlignments_1.20.0
[41] bitops_1.0-6             grid_3.6.1               jsonlite_1.6             gtable_0.3.0             DBI_1.0.0               
[46] magrittr_1.5             scales_1.0.0             stringi_1.4.3            tools_3.6.1              bit64_0.9-7             
[51] BSgenome_1.52.0          glue_1.3.1               purrr_0.3.2              hms_0.4.2                AnnotationDbi_1.46.0    
[56] colorspace_1.4-1         memoise_1.1.0  
```
