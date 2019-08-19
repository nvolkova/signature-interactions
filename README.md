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
 [1] deconstructSigs_1.8.0                   org.Hs.eg.db_3.8.2                      xlsx_0.6.1                             
 [4] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 GenomicFeatures_1.36.0                  AnnotationDbi_1.46.0                   
 [7] BSgenome.Hsapiens.UCSC.hg19_1.4.0       BSgenome_1.52.0                         rtracklayer_1.40.2                     
[10] dndscv_0.0.1.0                          bayesplot_1.6.0                         RColorBrewer_1.1-2                     
[13] greta_0.2.3                             reshape2_1.4.3                          ggplot2_3.1.1                          
[16] VariantAnnotation_1.30.0                Rsamtools_2.0.0                         Biostrings_2.52.0                      
[19] XVector_0.24.0                          SummarizedExperiment_1.14.0             DelayedArray_0.10.0                    
[22] BiocParallel_1.18.0                     matrixStats_0.54.0                      Biobase_2.44.0                         
[25] GenomicRanges_1.36.0                    GenomeInfoDb_1.20.0                     IRanges_2.18.0                         
[28] S4Vectors_0.22.0                        BiocGenerics_0.30.0                   

loaded via a namespace (and not attached):
 [1] httr_1.4.0               bit64_0.9-7              jsonlite_1.6             assertthat_0.2.1         xlsxjars_0.6.1          
 [6] blob_1.1.1               GenomeInfoDbData_1.2.1   progress_1.2.1           pillar_1.4.0             RSQLite_2.1.1           
[11] lattice_0.20-38          glue_1.3.1               reticulate_1.12          digest_0.6.18            colorspace_1.4-1        
[16] Matrix_1.2-17            plyr_1.8.4               XML_3.98-1.19            pkgconfig_2.0.2          biomaRt_2.40.0          
[21] zlibbioc_1.30.0          purrr_0.3.2              scales_1.0.0             whisker_0.3-2            tibble_2.1.1            
[26] withr_2.1.2              lazyeval_0.2.2           magrittr_1.5             crayon_1.3.4             memoise_1.1.0           
[31] MASS_7.3-51.4            tools_3.6.1              prettyunits_1.0.2        hms_0.4.2                stringr_1.4.0           
[36] munsell_0.5.0            ade4_1.7-13              compiler_3.6.1           rlang_0.3.4              RCurl_1.95-4.12         
[41] ggridges_0.5.1           rstudioapi_0.10          bitops_1.0-6             base64enc_0.1-3          labeling_0.3            
[46] gtable_0.3.0             DBI_1.0.0                R6_2.4.0                 GenomicAlignments_1.20.0 tfruns_1.4              
[51] dplyr_0.8.0.1            tensorflow_1.13.1        seqinr_3.4-5             bit_1.1-14               rJava_0.9-11            
[56] stringi_1.4.3            Rcpp_1.0.1               tidyselect_0.2.5         coda_0.19-2
```
