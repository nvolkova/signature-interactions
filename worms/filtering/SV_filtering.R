####################################################
##### Filtering of new BRASS vcf's    ##############
##### Hinxton, EMBL-EBI, Sept 2017    ##############
##### N. Volkova, nvolkova@ebi.ac.uk  ##############
####################################################

library(tidyr)
library(rtracklayer)
library(VariantAnnotation)
source('../../useful_functions.R')


# 1. Upload the variants and perform QC

data <- openxlsx::read.xlsx(xlsxFile = "Supplementary Table 1. Sample description for C.elegans experiments.xlsx", sheet = 2, cols = 1:8)
data$Sample <- as.character(data$Sample)
data$Genotype <- as.character(data$Genotype)
CD2Mutant <- sapply(1:nrow(data), function(i) {
  if (data$Type[i] == 'mut.acc.') return(paste0(data$Genotype[i],':',data$Generation[i]))
  return(paste0(data$Genotype[i],':',data$Mutagen[i],':',data$Drug.concentration[i]))
})
names(CD2Mutant) <- data$Sample
CD2Mutant <- CD2Mutant[sort(names(CD2Mutant))]

# Upload the VCFs
library(VariantAnnotation)
VCFPATH='/path/to/SV/VCFs'
# # FILTERING
raw.delly.vcf <- sapply(names(CD2Mutant), function(x) read_ce_vcf(paste0(VCFPATH,x,'.vcf')))
raw.delly.vcf <- raw.delly.vcf[!is.na(raw.delly.vcf)]
delly.vcf <- sapply(raw.delly.vcf, function(vcf) 
  vcf[geno(vcf)[['DV']][,1]>=10 & # variant support in test
        geno(vcf)[['DV']][,2]<1 &  # variant support in control
        geno(vcf)[['DR']][,1]+geno(vcf)[['DV']][,1] < 150 & # coverage
        geno(vcf)[['DR']][,2]+geno(vcf)[['DV']][,2] < 150 & # coverage
        granges(vcf)$FILTER=='PASS'])  # quality filter
barplot(sapply(delly.vcf,length))

# Remove MtDNA variants
delly.vcf <- lapply(delly.vcf, function(vcf) {
  if (length(vcf)==0) return(vcf)
  return(vcf[seqnames(vcf) %in% c("I","II","III","IV","V","X")])
})
# Make tables
delly.tables <- lapply(delly.vcf, function(vcf) {
  tmp <- data.frame(CHR1 = seqnames(granges(vcf)),
                    POS1 = start(granges(vcf)),
                    CHR2 = info(vcf)$CHR2,
                    POS2 = info(vcf)$END,
                    READS = info(vcf)$PE,
                    TYPE = info(vcf)$SVTYPE)
  rownames(tmp) <- names(vcf)
  return(tmp)
})


# 2. Filter out telomeric stuff etc

# Filter out telomeric stuff
# upload genome and get chromosome sizes
url <- 'ftp://ftp.ensembl.org/pub/release-96/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa.gz'
download.file(url = url, destfile = 'Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz', method = "auto", quiet=FALSE)
WBcel235 <- readDNAStringSet("Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz") # get worm reference genome
chr_sizes <- width(WBcel235)
names(chr_sizes) <- c("I","II","III","IV","MtDNA","V","X")
genome_size = sum(as.numeric(chr_sizes))

# telomeres
telomere <- matrix(1,nrow=6,ncol=4,dimnames=list(c("I","II","III","IV","V","X"),c("l.start","l.end","r.start","r.end")))
telomere[,2] <- c(436,400,200,47500,700,5000)
telomere[,3] <- c(15060000,15277300,13781900,
                  17492000,20923800,17718700)
telomere[,4] <- chr_sizes[-5]

telomere.delly <- list()

# get rid of variants intersecting with telomeres
for (x in names(delly.tables)) {
  bad <- NULL
  if (nrow(delly.tables[[x]])==0 || is.na(delly.tables[[x]])) next
  for (j in 1:nrow(delly.tables[[x]])) {
    if (delly.tables[[x]][j,2] < telomere[as.character(delly.tables[[x]][j,1]),2] || 
        delly.tables[[x]][j,4] > telomere[as.character(delly.tables[[x]][j,3]),3])
      bad <- c(bad,j)
  }
  if (length(bad)>0) {
    telomere.delly[[x]] <- delly.tables[[x]][bad,]
    delly.tables[[x]] <- delly.tables[[x]][-bad,]
    print(c("telomeres",length(bad),x))
  }
}
dimensions_filt <- sapply(delly.tables,nrow)


# 3. Deduplicate breakpoints

# Special treatment for MA samples: they may be related across generations so shall not be compared against each other within same genotype

mut.acc <- intersect(names(CD2Mutant)[is.na(data$Mutagen[match(names(CD2Mutant),data$Sample)])],names(delly.tables))
all_genotypes <- data$Genotype.new[match(mut.acc, data$Sample)]
all_generations <- data$Generation[match(mut.acc, data$Sample)]

allSV <- data.frame()
for (x in names(delly.tables)) {
  if(nrow(delly.tables[[x]]) == 0) next
  allSV <- rbind(allSV,
                 data.frame(delly.tables[[x]],
                            name = x,
                            stringsAsFactors = F))
}

allSV[,1] <- as.character(allSV[,1]) 
allSV[,3] <- as.character(allSV[,3]) 
allSV[,2] <- as.numeric(allSV[,2]) 
allSV[,4] <- as.numeric(allSV[,4]) # 41643

allSV.mut <- allSV[!(allSV$name %in% mut.acc),] # 34039

delly.tables.dedup <- list()
for (worm in names(delly.tables)) {
  
  if (worm %in% mut.acc) {
    worm_genotype <- data$Genotype.new[match(worm, data$Sample)]
    if (all_generations[match(worm,mut.acc)] > 1) 
      worms_to_compare <- setdiff(mut.acc[all_genotypes != worm_genotype | all_generations < 2],worm)
    else worms_to_compare <- setdiff(mut.acc,worm)
    allSV.compare <- allSV[allSV$name %in% worms_to_compare,]
    reference <- rbind(allSV.mut, allSV.compare)
    delly.tables.dedup[[worm]] <- filter.breakpoints(SV = delly.tables[[worm]], reference = reference)
  } else {
    delly.tables.dedup[[worm]] <- filter.breakpoints(SV = delly.tables[[worm]], reference = allSV[allSV$name != worm,])
  }
  print(worm)
}
dimensions_dedup <- sapply(delly.tables.dedup,nrow)
plot(dimensions_dedup,cex=0.2,main="numbers of SVs after deduplication of breakpoints")

# count 'em all
delly.filtered.bp.counts <- list()
for (i in seq_along(delly.tables.dedup)) {
  rearr.counts <- vector("numeric",5)
  names(rearr.counts) <- c("BND","DEL","INV","DUP","ALL")
  if (nrow(delly.tables.dedup[[i]]) < 1) {
    delly.filtered.bp.counts[[i]] <- rearr.counts
    next
  }
  if (is.na(delly.tables.dedup[[i]])) {
    delly.filtered.bp.counts[[i]] <- NA
    next
  }
  rearr.counts[1] <- sum(delly.tables.dedup[[i]]$TYPE == 'BND')
  rearr.counts[2] <- sum(delly.tables.dedup[[i]]$TYPE == 'DEL')
  rearr.counts[3] <- sum(delly.tables.dedup[[i]]$TYPE == 'INV')
  rearr.counts[4] <- sum(delly.tables.dedup[[i]]$TYPE == 'DUP')
  rearr.counts[5] <- sum(rearr.counts[1:4])
  delly.filtered.bp.counts[[i]] <- rearr.counts
}
names(delly.filtered.bp.counts) <- names(delly.tables.dedup)

# visualize class sizes
barplot(colSums(do.call("rbind",delly.filtered.bp.counts)))
# compare to non-bp deduplication
barplot(colSums(do.call("rbind",delly.filtcounts)))

# size filter 
barplot(sapply(delly.tables.dedup, nrow))
for (i in 1:length(delly.tables.dedup)) {
  ind <- which(as.character(delly.tables.dedup[[i]][,1])==as.character(delly.tables.dedup[[i]][,3]) & 
                 (as.numeric(delly.tables.dedup[[i]][,4])-as.numeric(delly.tables.dedup[[i]][,2])<400))
  if (length(ind)>0)
    delly.tables.dedup[[i]] <- delly.tables.dedup[[i]][-ind,,drop=F] 
}
barplot(sapply(delly.tables.dedup, nrow))


# 4. Now make sure that DELs and TDs are actually DELs and TDs (get them a p-value)

# Will need the paths to BigWig and bamstats files for each of the samples which have any deletions or tandem duplications.
# bamstats files were acuired as 
#  samtools idxstats sample.bam > sample.stat.dat

PATHTOBW <- 'path/to/bw/files'
PATHTOBAMSTATS <- 'path/to/bam/stats/files'

ref.reads1 <- 44773320 # overall number of mapped reads in reference genome (CD0001b)
ref.reads2 <- 31427593 # overall number of mapped reads in reference genome (CD0850b)

# Check TD and DEL, deduplicate
set1 <- names(CD2Mutant)[1:(match('CD0683b',names(CD2Mutant))-1)]
for (worm in names(delly.tables.dedup)[794:length(delly.tables.dedup)]) {
  if (worm %in% set1) {
    ref.reads <- ref.reads1
    ref.worm <- 'CD0001b'
  }
  else {
    ref.reads <- ref.reads2
    ref.worm <- 'CD0850b'
  }
   # upload stats files to get the number of mapped reads in the sample of interest
  if (nrow(delly.tables.dedup[[worm]]) == 0) next
   file <- paste(PATHTOBAMSTATS,worm,".stat.dat",sep="")
   alt.reads <- sum(read.table(file)$V3[1:6])
   r.ratio <- ref.reads / alt.reads
   # get the deletions
   dels <- delly.tables.dedup[[worm]][delly.tables.dedup[[worm]]$TYPE=="DEL",]
   if (nrow(dels)!=0) {
     goodness <- sapply(1:nrow(dels),function(j) {
       big.ranges <- c(as.numeric(as.character(dels[j,2])),as.numeric(as.character(dels[j,4])))
       chr <- as.character(dels[j,1])
       if (big.ranges[2]-big.ranges[1]<100000) {
         unit <- 10**round(log10(big.ranges[2]-big.ranges[1])-1)
         if (big.ranges[2]+unit > telomere[chr,'r.start']) big.ranges[2] <- telomere[chr,'r.start'] - unit
         if (big.ranges[1]-unit < telomere[chr,'l.end']) big.ranges[1] <- telomere[chr,'l.end'] + unit
         big.ranges <- GRanges(seqnames = chr, ranges = IRanges(start=big.ranges[1]-unit, end=big.ranges[2]+unit))
         x <- seq(from=start(big.ranges), to=end(big.ranges), by = unit/10)
         region <- rtracklayer::import(paste0(PATHTOBW,'/',worm,'.merged.bw'),which=big.ranges)
         region.ref <- rtracklayer::import(paste0(PATHTOBW,'/',ref.worm,".merged.bw"),which=big.ranges)
         bins <- GRanges(seqnames=dels[j,1],
                         ranges=IRanges(start=x[-length(x)]+1,end=x[-1]),
                         seqinfo=seqinfo(region))
         numvar <- mcolAsRleList(x=region,varname="score")
         numvar.ref <- mcolAsRleList(x=region.ref,varname="score")
         numvar.ref <- numvar.ref[names(numvar)]
         points <- as.numeric(binnedAverage(bins,numvar,varname="score",na.rm=T)$score) + 0.1
         points[is.na(points)] <- 1
         points.ref <- as.numeric(binnedAverage(bins,numvar.ref,varname="score")$score) + 0.1
         points.ref[is.na(points.ref)] <- 1
         
         x <- x[-length(x)]
         logfold <- log(points / points.ref * r.ratio)
         normal = c(1:10,(length(logfold)-9):length(logfold))
         odd = c(10:(length(logfold)-10))
         if (length(which(is.nan(points)))>0) {
           print(c("NaN points!!!",file))
           points.ref <- points.ref[-which(points=="NaN")]
           normal <- setdiff(normal,which(points=="NaN"))
           odd <- setdiff(odd,which(points=="NaN"))
           x <- x[-which(points=="NaN")]
           points <- points[-which(points=="NaN")]
         }
         if (length(which(is.na(points.ref)))>0) {
           print(c("NA points in ref!!!",file))
           normal <- setdiff(normal,which(is.na(points.ref)))
           odd <- setdiff(odd,which(is.na(points.ref)))
           x <- x[-which(is.na(points.ref))]
           points <- points[-which(is.na(points.ref))]
           points.ref <- points.ref[-which(is.na(points.ref))]
         }
         logfold <- log(points / points.ref * r.ratio)
         if (length(logfold)==0) return(NA)
         return(wilcox.test(logfold[normal], logfold[odd], alternative="greater")$p.value)
       }
       else {
         unit <- 10**round(log10(big.ranges[2]-big.ranges[1]) - 1)
         if (big.ranges[2]+unit > telomere[chr,'r.start']) big.ranges[2] <- telomere[chr,'r.start'] - unit
         if (big.ranges[1]-unit < telomere[chr,'l.end']) big.ranges[1] <- telomere[chr,'l.end'] + unit
         big.ranges <- GRanges(seqnames = chr, ranges = IRanges(start=big.ranges[1]-unit, end=big.ranges[2]+unit))
         x <- seq(from=start(big.ranges), to=end(big.ranges), by = unit/10)
         region <- rtracklayer::import(paste0(PATHTOBW,'/',worm,'.merged.bw'),which=big.ranges)
         region.ref <- rtracklayer::import(paste0(PATHTOBW,'/',ref.worm,".merged.bw"),which=big.ranges)
         bins <- GRanges(seqnames=dels[j,1],
                         ranges=IRanges(start=x[-length(x)]+1,end=x[-1]),
                         seqinfo=seqinfo(region))
         numvar <- mcolAsRleList(x=region,varname="score")
         numvar.ref <- mcolAsRleList(x=region.ref,varname="score")
         numvar.ref <- numvar.ref[names(numvar)]
         points <- as.numeric(binnedAverage(bins,numvar,varname="score",na.rm=T)$score) + 0.1
         points[is.na(points)] <- 1
         points.ref <- as.numeric(binnedAverage(bins,numvar.ref,varname="score")$score) + 0.1
         points.ref[is.na(points.ref)] <- 1
         if (length(which(is.nan(points)))>0) {
           print(c("NaN points!!!",file))
           points.ref <- points.ref[-which(points=="NaN")]
           points <- points[-which(points=="NaN")]
         }
         if (length(which(is.nan(points.ref)))>0) {
           print(c("NaN points in ref!!!",file))
           points <- points[-which(points.ref=="NaN")]
           points.ref <- points.ref[-which(points.ref=="NaN")]
         }
         logfold <- log(points / points.ref * r.ratio)
         return(wilcox.test(logfold,mu=0,alternative = 'less')$p.value) 	
       }
     })
   }
   else goodness <- NULL
   del.td.pvalue <- rep(NA,nrow(delly.tables.dedup[[worm]]))
   del.td.pvalue[delly.tables.dedup[[worm]][,'TYPE']=="DEL"] <- as.numeric(goodness)
   delly.tables.dedup[[worm]]$del.td.pvalue <- del.td.pvalue
   
   # get the deletions
   tds <- delly.tables.dedup[[worm]][delly.tables.dedup[[worm]]$TYPE=="DUP",]
   if (nrow(tds)!=0) {
     goodness <- sapply(1:nrow(tds),function(j) {
       big.ranges <- c(as.numeric(as.character(tds[j,2])),as.numeric(as.character(tds[j,4])))
       chr <- as.character(tds[j,1])
       if (big.ranges[2]-big.ranges[1]<100000) {
         unit <- 10**round(log10(big.ranges[2]-big.ranges[1])-1)
         if (big.ranges[2]+unit > telomere[chr,'r.start']) big.ranges[2] <- telomere[chr,'r.start'] - unit
         if (big.ranges[1]-unit < telomere[chr,'l.end']) big.ranges[1] <- telomere[chr,'l.end'] + unit
         big.ranges <- GRanges(seqnames = chr, ranges = IRanges(start=big.ranges[1]-unit, end=big.ranges[2]+unit))
         x <- seq(from=start(big.ranges), to=end(big.ranges), by = unit/10)
         region <- rtracklayer::import(paste0(PATHTOBW,'/',worm,'.merged.bw'),which=big.ranges)
         region.ref <- rtracklayer::import(paste0(PATHTOBW,'/',ref.worm,".merged.bw"),which=big.ranges)
         bins <- GRanges(seqnames=tds[j,1],
                         ranges=IRanges(start=x[-length(x)]+1,end=x[-1]),
                         seqinfo=seqinfo(region))
         numvar <- mcolAsRleList(x=region,varname="score")
         numvar.ref <- mcolAsRleList(x=region.ref,varname="score")
         numvar.ref <- numvar.ref[names(numvar)]
         points <- as.numeric(binnedAverage(bins,numvar,varname="score",na.rm=T)$score) + 0.1
         points[is.na(points)] <- 1
         points.ref <- as.numeric(binnedAverage(bins,numvar.ref,varname="score")$score) + 0.1
         points.ref[is.na(points.ref)] <- 1
         x <- x[-length(x)]
         logfold <- log(points / points.ref * r.ratio)
         normal = c(1:10,(length(logfold)-9):length(logfold))
         odd = c(10:(length(logfold)-10))
         if (length(which(is.nan(points)))>0) {
           print(c("NaN points!!!",file))
           points.ref <- points.ref[-which(points=="NaN")]
           normal <- setdiff(normal,which(points=="NaN"))
           odd <- setdiff(odd,which(points=="NaN"))
           x <- x[-which(points=="NaN")]
           points <- points[-which(points=="NaN")]
         }
         if (length(which(is.na(points.ref)))>0) {
           print(c("NA points in ref!!!",file))
           normal <- setdiff(normal,which(is.na(points.ref)))
           odd <- setdiff(odd,which(is.na(points.ref)))
           x <- x[-which(is.na(points.ref))]
           points <- points[-which(is.na(points.ref))]
           points.ref <- points.ref[-which(is.na(points.ref))]
         }
         logfold <- log(points / points.ref * r.ratio)
         if (length(logfold)==0) return(NA)
         return(wilcox.test(logfold[normal], logfold[odd], alternative="less")$p.value)
       }
       else {
         unit <- 10**round(log10(big.ranges[2]-big.ranges[1]) - 1)
         if (big.ranges[2]+unit > telomere[chr,'r.start']) big.ranges[2] <- telomere[chr,'r.start'] - unit
         if (big.ranges[1]-unit < telomere[chr,'l.end']) big.ranges[1] <- telomere[chr,'l.end'] + unit
         big.ranges <- GRanges(seqnames = chr, ranges = IRanges(start=big.ranges[1]-unit, end=big.ranges[2]+unit))
         x <- seq(from=start(big.ranges), to=end(big.ranges), by = unit/10)
         region <- rtracklayer::import(paste0(PATHTOBW,'/',worm,'.merged.bw'),which=big.ranges)
         region.ref <- rtracklayer::import(paste0(PATHTOBW,'/',ref.worm,".merged.bw"),which=big.ranges)
         bins <- GRanges(seqnames=tds[j,1],
                         ranges=IRanges(start=x[-length(x)]+1,end=x[-1]),
                         seqinfo=seqinfo(region))
         numvar <- mcolAsRleList(x=region,varname="score")
         numvar.ref <- mcolAsRleList(x=region.ref,varname="score")
         numvar.ref <- numvar.ref[names(numvar)]
         points <- as.numeric(binnedAverage(bins,numvar,varname="score",na.rm=T)$score)+0.1
         points[is.na(points)] <- 1
         points.ref <- as.numeric(binnedAverage(bins,numvar.ref,varname="score")$score)+0.1
         points.ref[is.na(points.ref)] <- 1
         if (length(which(is.nan(points)))>0) {
           print(c("NaN points!!!",file))
           points.ref <- points.ref[-which(points=="NaN")]
           points <- points[-which(points=="NaN")]
         }
         if (length(which(is.nan(points.ref)))>0) {
           print(c("NaN points in ref!!!",file))
           points <- points[-which(points.ref=="NaN")]
           points.ref <- points.ref[-which(points.ref=="NaN")]
         }
         logfold <- log(points / points.ref * r.ratio)
         return(wilcox.test(logfold,mu=0,alternative='greater')$p.value) 	
       }
     })
   }
   else goodness <- NULL
   if (!('del.td.pvalue' %in% colnames(delly.tables.dedup[[worm]]))) {
     del.td.pvalue <- rep(NA,nrow(delly.tables.dedup[[worm]]))
     del.td.pvalue[delly.tables.dedup[[worm]][,'TYPE']=="DUP"] <- as.numeric(goodness)
     delly.tables.dedup[[worm]]$del.td.pvalue <- del.td.pvalue
   }
   else {
     delly.tables.dedup[[worm]]$del.td.pvalue[delly.tables.dedup[[worm]][,'TYPE']=="DUP"] <- as.numeric(goodness)
   }
   print(worm)
}


# p-values - check NAs
for (worm in names(delly.tables.dedup)) {
  dels <- which(delly.tables.dedup[[worm]]$TYPE=="DUP")
  tds <- which(delly.tables.dedup[[worm]]$TYPE=="DEL")
  if (length(dels)>0) {
    bad.del <- which(is.na(delly.tables.dedup[[worm]][dels,7]))
    if (length(bad.del)>0) delly.tables.dedup[[worm]] <- delly.tables.dedup[[worm]][-dels[bad.del],]
  }
  if (length(tds)>0) {
    bad.td <- which(is.na(delly.tables.dedup[[worm]][tds,7]))
    if (length(bad.td)>0) delly.tables.dedup[[worm]] <- delly.tables.dedup[[worm]][-tds[bad.td],]
  }
}
delly.tables <- delly.tables.dedup[sapply(delly.tables.dedup,nrow)>0]
barplot(sapply(delly.tables.dedup,nrow)) # no difference

# do multiple testing correction - Benjamini-Hochberg
pval.del <- NULL
pval.td <- NULL
for (x in names(delly.tables.dedup)) {
  pval.del <- c(pval.del,delly.tables.dedup[[x]][delly.tables.dedup[[x]]$TYPE=="DEL",7])
  pval.td <- c(pval.td,delly.tables.dedup[[x]][delly.tables.dedup[[x]]$TYPE=="DUP",7])
}
pval.del.adj <- p.adjust(pval.del,method="BH")
pval.td.adj <- p.adjust(pval.td,method="BH")
for (x in names(delly.tables.dedup)) {
  td.length <- length(which(delly.tables.dedup[[x]]$TYPE=="DUP"))
  del.length <- length(which(delly.tables.dedup[[x]]$TYPE=="DEL"))
  if (del.length>0) {
    delly.tables.dedup[[x]]$del.td.pvalue[delly.tables.dedup[[x]]$TYPE=="DEL"] <- pval.del.adj[1:del.length]
    pval.del.adj <- pval.del.adj[-c(1:del.length)]
  }
  if (td.length>0) {
    delly.tables.dedup[[x]]$del.td.pvalue[delly.tables.dedup[[x]]$TYPE=="DUP"] <- pval.td.adj[1:td.length]
    pval.td.adj <- pval.td.adj[-c(1:td.length)]
  }
}

PATHTOCLUST='/path/where/to/write/clustered/tables'
# CLUSTER
source("classification_script.R")
for (worm in names(delly.tables.dedup)) {
  d <- delly.tables.dedup[[worm]]
  if (nrow(d) == 0) next
  output_file = paste(PATHTOCLUST,'/',worm,".clust_mat",sep="")
  d[,1] = as.character(d[,1]) # chomosomes
  d[,3] = as.character(d[,3]) # chomosomes
  d[,2] = as.numeric(as.character(d[,2]))
  d[,4] = as.numeric(as.character(d[,4]))
  d[,6] = as.character(d[,6])
  if (nrow(d) == 1) {
    ct = 1
    # Get the footprint info
    pos = c(d[,2], d[,4])
    chrs = c(d[,1], d[,3])
    res = get_footprints(pos, chrs)
    footprint_idx = sprintf("%s.chr%s.%s", c(ct, ct), chrs, res$footprint_idx)
    footprint_bounds = res$footprint_bounds
    
    write.table(
      data.frame(
        d,
        clust = ct,
        clust_size = sapply(ct, function(x) sum(x == ct)),
        fp_l = footprint_idx[1:nrow(d)],
        fp_h = footprint_idx[-(1:nrow(d))]
      ),
      output_file,
      quote = F,
      sep = "\t"
    )
  } else {
    ct = clust_rgs_new(d)$cutree
    
    # Get the footprint info
    pos = c(d[,2], d[,4])
    chrs = c(d[,1], d[,3])
    res = get_footprints(pos, chrs)
    footprint_idx = sprintf("%s.chr%s.%s", c(ct, ct), chrs, res$footprint_idx)
    footprint_bounds = res$footprint_bounds
    
    write.table(
      data.frame(
        d,
        clust = ct,
        clust_size = sapply(ct, function(x) sum(x == ct)),
        fp_l = footprint_idx[1:nrow(d)],
        fp_h = footprint_idx[-(1:nrow(d))]
      ),
      output_file,
      quote = F,
      sep = "\t"
    )
  }
}

delly.tables.dedup <- delly.tables.dedup[sapply(delly.tables.dedup,nrow)>0]

# reading reclustered SVs
SVclust <- list()
for (worm in names(delly.tables.dedup))
{
  if(nrow(delly.tables.dedup[[worm]]) > 0) {
  file = paste(PATHTOCLUST,'/',worm,".clust_mat",sep="")
  SVclust[[worm]] <- read.table(file=file, header=T, sep = "\t")
  }
}

# add sample name to all SV tables
for (worm in names(SVclust)){
  SVclust[[worm]] <- cbind(SVclust[[worm]][,1:10],Sample=as.character(rep(worm,nrow(SVclust[[worm]]))),
                           del.td.pvalue = delly.tables.dedup[[worm]]$del.td.pvalue)
  
}

# Assessing types of clusters
rearr.count.final.dedup <- list()
for (i in seq_along(SVclust)) {
  SVclust[[i]] <- cbind(SVclust[[i]],clust.type=as.character(rep("some",nrow(SVclust[[i]]))))
  SVclust[[i]]$clust.type <- as.character(SVclust[[i]]$clust.type)
  rearr.counts <- vector("numeric",8)
  names(rearr.counts) <- c("TD","DEL","INV","COMPLEX","TRSL","INTCHR","FOLDBACK","MOVE") # TRSL = copypaste, MOVE = TRSL with deletion
  for (j in unique(SVclust[[i]]$clust)) {
    which(SVclust[[i]]$clust==j) -> clust_ind
    bp.types <- as.character(SVclust[[i]]$TYPE[clust_ind])
    
    # 2 INTERSECTING pairs of breakpoints DEL, DEL, TD - translocation
    # any >2 pairs - complex
    if (length(clust_ind)>2) {
      if (length(clust_ind)==3 & 
          (length(which(bp.types=="DEL"))==2) & 
          ("DUP" %in% bp.types)) {
        rearr.counts["MOVE"] = rearr.counts["MOVE"] + 1
        SVclust[[i]]$clust.type[clust_ind] <- rep("MOVE",length(clust_ind))        
      } else {
        rearr.counts["COMPLEX"] = rearr.counts["COMPLEX"] + 1
        SVclust[[i]]$clust.type[clust_ind] <- rep("COMPLEX",length(clust_ind))
      }
    }
    
    # 2 pairs of breakpoints: inversions, interchromosomal
    if (length(clust_ind)==2) {
      if (length(setdiff(c("INV","INV"),bp.types))==0) {
        rearr.counts["INV"] = rearr.counts["INV"] + 1
        SVclust[[i]]$clust.type[clust_ind] <- c("INV","INV")
      }
      else if (length(setdiff(bp.types, c("BND","BND")))==0) {
        SVclust[[i]]$clust.type[clust_ind] <- c("INTCHR","INTCHR")
        rearr.counts["INTCHR"] = rearr.counts["INTCHR"] + 1
      }
      else if (length(setdiff(bp.types,c("DUP","DEL")))==0) {
        dist1 <- SVclust[[i]]$POS1[clust_ind][2]-SVclust[[i]]$POS1[clust_ind][1]
        dist2 <- SVclust[[i]]$POS2[clust_ind][2]-SVclust[[i]]$POS2[clust_ind][1]
        dist <- SVclust[[i]]$POS2[clust_ind][2]-SVclust[[i]]$POS1[clust_ind][1]
        if ((abs(dist1)<(0.1*dist) & bp.types[which.max(SVclust[[i]]$POS2[clust_ind])]=="DUP") || 
            (abs(dist2)<(0.2*dist) & bp.types[1]=="DUP") ||
            (abs(dist1)<(0.1*dist) & abs(dist2)<(0.1*dist))) {
          SVclust[[i]]$clust.type[clust_ind] <- rep("TRSL",2)
          rearr.counts["TRSL"] = rearr.counts["TRSL"] + 1
        } else {
          SVclust[[i]]$clust.type[clust_ind] <- c("COMPLEX","COMPLEX")
          rearr.counts["COMPLEX"] = rearr.counts["COMPLEX"] + 1
        }
      }
      else if (length(which(bp.types=="DEL"))==2) {
        if (length(which(SVclust[[i]][clust_ind,12]>0.05))==0) {
          SVclust[[i]]$clust.type[clust_ind] <- c("COMPLEX","COMPLEX")
          rearr.counts["COMPLEX"] = rearr.counts["COMPLEX"] + 1
        }
        else if (length(which(SVclust[[i]][clust_ind,12]>0.05))==1) {
          SVclust[[i]]$clust.type[clust_ind] <- c("DEL","DEL")
          rearr.counts["DEL"] = rearr.counts["DEL"] + 1
        }
      }
      else if (length(which(bp.types=="DUP"))==2) {
        if (length(which(SVclust[[i]][clust_ind,12]>0.05))==0) {
          SVclust[[i]]$clust.type[clust_ind] <- c("COMPLEX","COMPLEX")
          rearr.counts["COMPLEX"] = rearr.counts["COMPLEX"] + 1
        }
        else if (length(which(SVclust[[i]][clust_ind,12]>0.05))==1) {
          SVclust[[i]]$clust.type[clust_ind] <- c("TD","TD")
          rearr.counts["TD"] = rearr.counts["TD"] + 1
        }
      }
      else {
        rearr.counts["COMPLEX"] = rearr.counts["COMPLEX"] + 1
        SVclust[[i]]$clust.type[clust_ind] <- c("COMPLEX","COMPLEX")
      }
    }
    
    if (length(clust_ind)==1) {
      if (bp.types=="DUP"  && length(which(SVclust[[i]]$del.td.pvalue[clust_ind]>0.05))==0) {
        SVclust[[i]]$clust.type[clust_ind] <- "TD"
        rearr.counts["TD"] = rearr.counts["TD"] + 1
      }
      if (bp.types=="DEL" && length(which(SVclust[[i]]$del.td.pvalue[clust_ind]>0.05))==0)
      {
        rearr.counts["DEL"] = rearr.counts["DEL"] + 1
        SVclust[[i]]$clust.type[clust_ind] <- "DEL"
      }
      if (bp.types=="INV")
      {
        rearr.counts["FOLDBACK"] = rearr.counts["FOLDBACK"] + 1
        SVclust[[i]]$clust.type[clust_ind] <- "FOLDBACK"
      }
      if (bp.types=="BND") {
        rearr.counts["INTCHR"] = rearr.counts["INTCHR"] + 1
        SVclust[[i]]$clust.type[clust_ind] <- "INTCHR"
      }
    }
    
  }
  print(i)
  rearr.count.final.dedup[[i]] <- rearr.counts
}
names(rearr.count.final.dedup) <- names(SVclust)

# visualize class sizes
colSums(do.call("rbind",delly.filtcounts))
#   TD      DEL      INV  COMPLEX     TRSL   INTCHR FOLDBACK  MOVE    ALL 
#  150     167       99       128       32      125       90    1     792
colSums(do.call("rbind",rearr.count.final.dedup))
#  TD      DEL      INV  COMPLEX     TRSL   INTCHR FOLDBACK     MOVE    ALL
# 767      604      257      320      100      337      265       32   2682

delly.filtcounts <- rearr.count.final.dedup
SVclust -> delly.SVclust

delly.tables.dedup <- delly.tables.dedup[sapply(delly.tables.dedup,nrow)>0]
for (worm in names(delly.tables.dedup)) {
  if ('del.td.pvalue' %in% colnames(delly.tables.dedup[[worm]])) next
  delly.tables.dedup[[worm]]$del.td.pvalue <- NA
}

SVclust <- SVclust[sapply(SVclust,nrow)>0]
for (worm in names(SVclust)) {
  if ('del.td.pvalue' %in% colnames(SVclust[[worm]])) next
  SVclust[[worm]]$del.td.pvalue <- NA
}

save(delly.SVclust, delly.tables.dedup, delly.filtcounts, file='filtered_SV.RData')

filt.delly.vcf <- lapply(names(delly.vcf), function(z) delly.vcf[[z]][rownames(delly.tables.dedup[[z]])])
names(filt.delly.vcf) <- names(delly.vcf)

save(delly.vcf, file = 'Deduplicated_and_filtered_DELLY_vcfs.Rds')

