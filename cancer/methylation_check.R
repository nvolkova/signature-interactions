# methylation analysis

library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
annotation.table = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

genes <- openxlsx::read.xlsx('Supplementary Table 4. TCGA samples with mutations in DNA pathways.xlsx', sheet = 2)

cpgs <- list()
for (j in as.character(genes$Gene[genes$CORE])) {
  cpgs[[j]] <- rownames(annotation.table)[grep(j,annotation.table$UCSC_RefGene_Name)]
}

methylated <- list()
for (n in tmp$Gene) {
  methylated[[n]] <- data.frame()
}

# This needs methylation tables from TCGA which already have only necessary genes included
# they are called 'CANCERTYPE_methylation_DNA_repair_core.RData'
PATHTOMETH='/path/to/methylation/tables/per/cancer/'
f <- list.files(PATHTOMETH)
cancernames <- sapply(f, function(x) unlist(strsplit(x,split = '[_]'))[1])

no.tests <- NULL
for (can in cancernames) {

  load(paste0(PATHTOMETH,can,'_methylation_DNA_repair_core.RData'))

  for (n in intersect(names(methyl), colnames(metadata))) {
  
    expr <- metadata[,n][match(substr(colnames(methyl[[n]]),1,12),substr(metadata$tumour,1,12))]
    if (sum(!is.na(expr)) == 0) next
    methyl[[n]] <- methyl[[n]][sapply(rownames(methyl[[n]]), function(x) sum(!is.na(methyl[[n]][x,])) > 0),]
    relevant_cpgs <- rownames(methyl[[n]])[sapply(rownames(methyl[[n]]),
                                                  function(x) cor(expr,
                                                                  as.numeric(methyl[[n]][x,]),
                                                                  use = 'complete.obs')) < -0.5]

    if (length(relevant_cpgs) > 0) {
      betas <- apply(methyl[[n]][relevant_cpgs,,drop = F],2,median, na.rm = T)
      unmethyl <- which(betas < 0.1)
      unm_mean <- mean(expr[unmethyl], na.rm = T)
      unm_sd <- sd(expr[unmethyl], na.rm = T)
      zcores <- (expr - unm_mean) / unm_sd
      meth <- which(betas > 0.2)
      met_mean <- mean(zcores[meth], na.rm = T)
      if (is.na(met_mean)) next
      while (length(meth) > 3 & (met_mean > qnorm(0.05) | 
             mean(expr[meth], na.rm = T) / mean(expr[unmethyl], na.rm = T) > 0.5)) {
        meth <- meth[-which.min(betas[meth])]
        met_mean <- mean(zcores[meth], na.rm = T)
      }
      if (length(meth)>2 & sum(!is.na(zcores[meth]))>2 & sum(!is.na(zcores[unmethyl]))>2) {
        no.tests <- c(no.tests,t.test(zcores[meth],zcores[unmethyl], alternative = 'less')$p.value)
        if (met_mean < qnorm(0.05) & t.test(zcores[meth],zcores[unmethyl], alternative = 'less')$p.value < 0.05 & 
            max(betas[meth]) > 0.5 & mean(expr[meth], na.rm = T) / mean(expr[unmethyl], na.rm = T) < 0.5)
        methylated[[n]] <- rbind(methylated[[n]],
                               cbind(names(betas[meth]), betas[meth]))
      }
    }
  
  }
  
}
save(methylated, file = 'methylation.list.RData')
