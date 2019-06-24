read_ce_vcf <- function(file) {
  out <- tryCatch(
    {
      readVcf(file, genome="WBcel235") 
    },
    error=function(cond) {
      message(paste("\nBad file in: ", file))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(paste("\nFile caused a warning: ", file))
      message("Here's the original warning message:")
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    }
  )    
  return(out)
}

genoToCoverage <- function(vcf){
  rowSums(sapply(grep("F|R.Z",names(geno(vcf)), value=TRUE), function(field){
    geno(vcf)[[field]]
  }, simplify="array"), dims=2) # number of all reads covering a particular position (substitution)
}

altCounts <- function(vcf){
  altAllele <- as.character(unlist(alt(vcf)))
  alleleCounts <- sapply(grep("F|R.Z",names(geno(vcf)), value=TRUE), function(field){
    geno(vcf)[[field]]
  }, simplify="array")
  t(sapply(seq_along(altAllele), function(i){
    alleleCounts[i,"TUMOUR",paste0(c("F","R"),altAllele[i], "Z")]
  }))
}


read_ce_table <- function(file, ...) {
  out <- tryCatch(
    {
      read.delim(file, ...) 
    },
    error=function(cond) {
      message(paste("\nBad file in: ", file))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(paste("\nFile caused a warning: ", file))
      message("Here's the original warning message:")
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    }
  )    
  return(out)
}