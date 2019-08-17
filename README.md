# signature-interactions
Codes and data for the manuscript "Mutational signatures are jointly shaped by DNA damage and repair"
(to be submitted in June 2019)

Human analysis part requires a matrix with TCGA variants, and a metadata matrix with sample names, cancer types, and FPKM expression values for DNA repair genes.

Session information for the analyses:

R version 3.5.1 (2018-07-02)

Matrix products: default

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base   

other attached packages:
 [1] reshape2_1.4.3              ggplot2_3.1.1               VariantAnnotation_1.30.0    Rsamtools_2.0.0            
 [5] Biostrings_2.52.0           XVector_0.24.0              SummarizedExperiment_1.14.0 DelayedArray_0.10.0        
 [9] BiocParallel_1.18.0         matrixStats_0.54.0          Biobase_2.44.0              GenomicRanges_1.36.0       
[13] GenomeInfoDb_1.20.0         IRanges_2.18.0              S4Vectors_0.22.0            BiocGenerics_0.30.0 
[17] greta_0.3.0.9002

loaded via a namespace (and not attached):
[1] progress_1.2.1           tidyselect_0.2.5         purrr_0.3.2              lattice_0.20-38          colorspace_1.4-1        
[6] rtracklayer_1.40.2       GenomicFeatures_1.36.0   blob_1.1.1               XML_3.98-1.19            rlang_0.3.4             
[11] pillar_1.4.0             withr_2.1.2              glue_1.3.1               DBI_1.0.0                bit64_0.9-7             
[16] GenomeInfoDbData_1.2.1   plyr_1.8.4               stringr_1.4.0            zlibbioc_1.30.0          munsell_0.5.0           
[21] gtable_0.3.0             memoise_1.1.0            biomaRt_2.40.0           AnnotationDbi_1.46.0     Rcpp_1.0.1              
[26] scales_1.0.0             BSgenome_1.52.0          bit_1.1-14               hms_0.4.2                digest_0.6.18           
[31] stringi_1.4.3            dplyr_0.8.0.1            grid_3.6.1               tools_3.6.1              bitops_1.0-6            
[36] magrittr_1.5             RCurl_1.95-4.12          lazyeval_0.2.2           RSQLite_2.1.1            tibble_2.1.1            
[41] crayon_1.3.4             pkgconfig_2.0.2          Matrix_1.2-17            prettyunits_1.0.2        assertthat_0.2.1        
[46] httr_1.4.0               rstudioapi_0.10          R6_2.4.0                 GenomicAlignments_1.20.0 compiler_3.6.1          
[51] codetools_0.2-15         listenv_0.7.0            future_1.12.0            jsonlite_1.6             coda_0.19-2         
[56] tfruns_1.4               whisker_0.3-2            reticulate_1.11.1        parallel_3.5.1           base64enc_0.1-3
[61] tensorflow_1.10          globals_0.12.4

