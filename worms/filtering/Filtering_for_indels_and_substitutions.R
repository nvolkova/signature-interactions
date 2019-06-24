#Doublechecking and speeding up filtering of worm data for mutational signatures project. Potentially also proper spectra acquisition and clustering of them.

## 0. Preparations

#Libraries
library(VariantAnnotation)
library(BSgenome)
library(deepSNV)
library(Biostrings)
library(dplyr)
worm_ref_genome <- "BSgenome.Celegans.UCSC.ce11"
library(worm_ref_genome, character.only = TRUE)

#Paths
DATAPATHSUBS = '/path/to/caveman/vcfs/'
DATAPATHINDEL = '/path/to/pindel/vcfs/'
PROCDATAPATH = "/path/to/filtered/data"
source("read_ce_vcf.R")
source("../../useful_functions.R")

# 1. UPLOADING DATA

## 1.0 Prepare the names and descriptions of samples

data <- openxlsx::read.xlsx(xlsxFile = "Supplementary Table 1. Sample description for C.elegans experiments.xlsx", sheet = 2, cols = 1:8)
data$Sample <- as.character(data$Sample)
data$Genotype <- as.character(data$Genotype)
CD2Mutant <- sapply(1:nrow(data), function(i) {
  if (data$Type[i] == 'mut.acc.') return(paste0(data$Genotype[i],':',data$Generation[i]))
  return(paste0(data$Genotype[i],':',data$Mutagen[i],':',data$Drug.concentration[i]))
})
names(CD2Mutant) <- data$Sample
CD2Mutant <- CD2Mutant[sort(names(CD2Mutant))]
# mutation accumulation samples
mut.acc <- names(CD2Mutant)[is.na(data$Mutagen)]

#Upload the genome (just in case)
url <- 'ftp://ftp.ensembl.org/pub/release-96/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa.gz'
download.file(url = url, destfile = 'Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz', method = "auto", quiet=FALSE)
WBcel235 <- readDNAStringSet("Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz") # get worm reference genome
chr_lens <- width(WBcel235)
names(chr_lens) <- c("I","II","III","IV","MtDNA","V","X")

## 1.1 Upload substitutions

files <- paste0(DATAPATHSUBS,"/",names(CD2Mutant),".caveman_c.sub.vcf.gz",sep="")
vcfs <- sapply(files, read_ce_vcf)
names(vcfs) <- names(CD2Mutant)

# Handle `NA` and check the dimensions
if (length(which(sapply(vcfs,length)==0))>0) {
  print('This samples have zero length of subs vcf:')
  print(names(which(sapply(vcfs,length)==0)))
  vcfs[which(sapply(vcfs,length)==0)] <- NA
}
barplot(sapply(vcfs,length),main="Number of substitutions before filtering")

## 1.2 Upload indels

files <- paste0(DATAPATHINDEL,"/",names(CD2Mutant),".pindel.annot.vcf.gz",sep="")
pindels <- sapply(files, read_ce_vcf)
names(pindels) <- names(CD2Mutant)

# Make empty vcfs `NA` and check the dimensions
if (length(which(sapply(pindels,length)==0))>0) {
  print('This samples have zero length of indels vcf:')
  print(names(which(sapply(pindels,length)==0)))
  pindels[which(sapply(pindels,length)==0)] <- NA
}
barplot(sapply(pindels,length),main="Number of indels before filtering")

## 1.3 Structural variants

# Calling performed with DELLY-SV
# go to SV_filtering.R


# 2. FILTERING

## 2.1 Filtering indels

# Normal panel + very slight quality filter
normal_panel = names(CD2Mutant)[grep("N2:1",CD2Mutant[1:ind1])][1:6]
normal_erroneous <- pindels[normal_panel]
normal_erroneous <- sapply(normal_erroneous, function(vcf)
  vcf[geno(vcf)$PR[,"NORMAL"] + geno(vcf)$NR[,"NORMAL"] < 150 &     # maximal coverage on control
        geno(vcf)$PR[,"NORMAL"] + geno(vcf)$NR[,"NORMAL"] >= 10 &                # minimal coverage in control
        geno(vcf)$PR[,"TUMOUR"] + geno(vcf)$NR[,"TUMOUR"] >= 10 &                # minimal coverage in test
        geno(vcf)$PP[,"TUMOUR"] + geno(vcf)$NP[,"TUMOUR"] > 1])

# All the others put through quality filter
filtered_indels = lapply(pindels[-match(normal_panel,names(pindels))], function(vcf) {
  
  if (length(vcf)==0)
    return(vcf)
  if (is.na(vcf))
    return(NA)
  else return(vcf[geno(vcf)$PR[,"NORMAL"] + geno(vcf)$NR[,"NORMAL"] < 150 &     # maximal coverage on control 150
                    geno(vcf)$PR[,"NORMAL"] + geno(vcf)$NR[,"NORMAL"] >= 10 &                # minimal coverage in control
                    geno(vcf)$PR[,"TUMOUR"] + geno(vcf)$NR[,"TUMOUR"] >= 10 &                # minimal coverage in test
                    geno(vcf)$PP[,"NORMAL"] + geno(vcf)$NP[,"NORMAL"] <= 1 &                 # <=1 reads reporting variant in control
                    (geno(vcf)$PU[,"TUMOUR"] + geno(vcf)$NU[,"TUMOUR"]) / (geno(vcf)$PR[,"TUMOUR"] + geno(vcf)$NR[,"TUMOUR"]) >= 0.2 & # at least 20% of reads in test worm reporting the variant
                    geno(vcf)$PP[,"TUMOUR"] + geno(vcf)$NP[,"TUMOUR"] >= 5 &                  # at least 5 reads reporting the variant in test sample
                    geno(vcf)$PP[,"TUMOUR"] >=1 & geno(vcf)$NP[,"TUMOUR"] >= 1]) #&             # at least 1 read reporting variant in both directions
})
filtered_indels[normal_panel] <- normal_erroneous[normal_panel]
filtered_indels <- filtered_indels[sort(names(filtered_indels))]

# Repeat filter 
filtered_indels_rep_filt = lapply(filtered_indels[-match(normal_panel,names(filtered_indels))], function(vcf) {
  
  if (length(vcf)==0)
    return(vcf)
  if (is.na(vcf))
    return(vcf)
  return(vcf[info(vcf)$REP <= 18 ])                                                    # less than 18 repeat units at site (after that detection is bad)
})
barplot(sapply(filtered_indels[-match(normal_panel,names(filtered_indels))], length), main = 'Number of indels after QC')
barplot(sapply(filtered_indels_rep_filt,length), main = 'Number of indels after repeat filter')
filtered_indels[names(filtered_indels_rep_filt)] <- filtered_indels_rep_filt

# Duplicates removal: first - for mutagen exposed samples as they should be unrelated

# there will be special tratment for these
mut.acc <- mut.acc[mut.acc %in% names(filtered_indels)]
# fct to produce comprehensive indel labels
rename_indels <- function(vcf) {
  if (length(vcf)!=0 && !is.na(vcf))
    paste0(seqnames(granges(vcf)),":",start(granges(vcf)),"_",as.character(granges(vcf)$REF),"/",as.character(unlist(granges(vcf)$ALT)))
}
# merge stuff together and duplicate normal panel
full_MUT_indellist <- do.call("c",sapply(filtered_indels,rename_indels))
full_duplicates <- unique(full_MUT_indellist[duplicated(full_MUT_indellist)])

# fct that matches processed indel names
deduplicate <- function(vcf, ref) {
  if (length(vcf)==0) return(vcf)
  if (is.na(vcf)) return(vcf)
  variants_in_vcf <- rename_indels(vcf)
  bad <- variants_in_vcf %in% ref
  return(vcf[!bad])
}
# deduplicate
MA_indels_dedup <- lapply(filtered_indels[!(names(filtered_indels) %in% mut.acc)], function = deduplicate, ref = full_duplicates)

# Mutation accumulation samples
mut.acc <- setdiff(mut.acc, normal_panel)

all_genotypes <- data$Genotype.new[match(mut.acc, data$Sample)]
all_generations <- data$Generation[match(mut.acc, data$Sample)]

MA_indels_dedup <- list()
for (worm in mut.acc) {
  
  worm_genotype <- data$Genotype.new[match(worm, data$Sample)]
  if (all_generations[match(worm,mut.acc)] > 1) {
    worms_to_compare <- mut.acc[all_genotypes != worm_genotype | all_generations < 2]
    if (!(worm %in% worms_to_compare)) worms_to_compare <- c(worm,worms_to_compare)
    indellist <- c(do.call("c",sapply(filtered_indels[worms_to_compare],rename_indels)),MUT_indellist)
    duplicates <- unique(indellist[duplicated(indellist)])
    MA_indels_dedup[[worm]] <- deduplicate(vcf = filtered_indels[[worm]], ref = duplicates)
  }
  else {
    MA_indels_dedup[[worm]] <- deduplicate(vcf = filtered_indels[[worm]], ref = full_duplicates)
  }
  print(c(worm, CD2Mutant[worm], length(MA_indels_dedup[[worm]])))
  
}

# Merge them back together
# Size filter - smaller than 400
indels_dedup <- c(MUT_indels_dedup, MA_indels_dedup)
indels_dedup <- sapply(indels_dedup, function(vcf) {
  if (length(vcf)==0 || is.na(vcf)) return(vcf)
  return(vcf[info(vcf)$LEN<400])
})
# Visualize
barplot(sapply(indels_dedup[sort(names(indels_dedup))],length),main="Numbers of indels after deduplication")

# Cleanup
rm(MUT_indellist, MA_indels_dedup, MUT_indels_dedup, erroneous_indellist, normal_erroneous, duplicates, indels_dedup, filtered_indels)

## 2.2 Filtering substitutions

#Normal panel with some QC
normal_panel = names(CD2Mutant)[grep("N2:1",CD2Mutant[1:ind])][1:6]
normal_erroneous <- vcfs[normal_panel]
normal_erroneous <- sapply(normal_erroneous, function(vcf)
  vcf[genoToCoverage(vcf)[,"NORMAL"] < 150 &     # maximal coverage on control
        genoToCoverage(vcf)[,"NORMAL"] >= 15 &                # minimal coverage in control
        genoToCoverage(vcf)[,"TUMOUR"] >= 15 &                # minimal coverage in test
        altCounts(vcf)[,1] > 0 & altCounts(vcf)[,2] > 0])

#Quality control
vcfs_filter1 <- sapply(vcfs[-match(normal_panel,names(vcfs))], function(vcf) {
  if (is.na(vcf))
    return(NA)
  return(vcf[genoToCoverage(vcf)[,"NORMAL"] < 150 &          # no more than 150 reads in control
               geno(vcf)$PM[,"NORMAL"]==0 &                        # no report of variant in control
               geno(vcf)$PM[,"TUMOUR"]>=0.2 &                      # >= 20% reads in test worm report mutant allele?
               genoToCoverage(vcf)[,"NORMAL"] >= 15 &              # >= 15 reads in test&control
               genoToCoverage(vcf)[,"TUMOUR"] >= 15 &
               altCounts(vcf)[,1] > 0 & altCounts(vcf)[,2] > 0])   # at least 1 read reporting variant in both directions
})
vcfs_filter1[normal_panel] <- normal_erroneous
vcfs_filter1 <- vcfs_filter1[sort(names(vcfs_filter1))]
barplot(as.vector(sapply(vcfs_filter1,length)), main="Number of substitutions after quality filter")

#Indels filter
newindelFilter <- function(vcf, pind) {
  if (length(vcf)==0 || is.na(vcf) || length(pind)==0 || is.na(pind)) return(vcf)
  pind <- pind[info(pind)$PC %in% c("D","I") &
                 abs(width(granges(pind)$REF)-width(unlist(granges(pind)$ALT)))==1 &
                 (width(granges(pind)$REF)==1 | width(unlist(granges(pind)$ALT))==1)]
  hits <- findOverlaps(query = ranges(vcf), subject = ranges(pind))
  if (length(queryHits(hits))>0)
    return(vcf[-unique(queryHits(hits))])
  else return(vcf)
}
vcfs_filter12 <- sapply(names(vcfs_filter1), function(sampleName) newindelFilter(vcfs_filter1[[sampleName]], filtered_indels[[sampleName]]))
dimensions_after_filter12 <- as.vector(sapply(vcfs_filter12,length))
barplot(dimensions_after_filter12, main="Number of variants after indel-specific filter")

# Deduplication for mutagen exposed samples (vs all)
vcfs_filter12 <- vcfs_filter12[!is.na(vcfs_filter12)]
mut.acc <- mut.acc[mut.acc %in% names(vcfs_filter12)]
full_allVcfs <- do.call("rbind",vcfs_filter12)
full_duplicates <- unique(full_allVcfs[duplicated(full_allVcfs)])
# fct for subs deduplication
deduplicate <- function(vcf, ref) {
  if (length(vcf)==0 || is.na(vcf)) return(vcf)
  variants_in_vcf <- paste0(seqnames(granges(vcf)),":",start(granges(vcf)),"_",as.character(granges(vcf)$REF),"/",as.character(unlist(granges(vcf)$ALT)))
  variants_in_ref <- paste0(seqnames(granges(ref)),":",start(granges(ref)),"_",as.character(granges(ref)$REF),"/",as.character(unlist(granges(ref)$ALT)))
  bad <- variants_in_vcf %in% variants_in_ref
  return(vcf[!bad])
}
MUT_vcfs_dedup <- lapply(vcfs_filter12[!(names(vcfs_filter12) %in% mut.acc)], FUN=deduplicate, full_duplicates)
MUT_allVcfs <- do.call("rbind",vcfs_filter12[!is.na(vcfs_filter12)][!(names(vcfs_filter12[!is.na(vcfs_filter12)]) %in% mut.acc)]) # gather the rest for MA filtering

#Mutation accumulation
mut.acc <- setdiff(mut.acc, normal_panel)
all_genotypes <- data$Genotype.new[match(mut.acc, data$Sample)]
all_generations <- data$Generation[match(mut.acc, data$Sample)]
MA_vcfs_dedup <- list()
for (worm in mut.acc) {
  
  worm_genotype <- data$Genotype.new[match(worm, data$Sample)]
  if (all_generations[match(worm,mut.acc)] > 1) {
    worms_to_compare <- mut.acc[all_genotypes != worm_genotype | all_generations < 2]
    if (!(worm %in% worms_to_compare)) worms_to_compare <- c(worm,worms_to_compare)
    allVcfs <- rbind(do.call("rbind",vcfs_filter12[worms_to_compare]),MUT_allVcfs)
    duplicates <- unique(allVcfs[duplicated(allVcfs)])
    MA_vcfs_dedup[[worm]] <- deduplicate(vcf = vcfs_filter12[[worm]],ref = duplicates)
  }
  else {
    MA_vcfs_dedup[[worm]] <- deduplicate(vcf = vcfs_filter12[[worm]],ref = full_duplicates)
  }
  print(c(worm,CD2Mutant[worm],length(MA_vcfs_dedup[[worm]])))
  
}

# Bring them back together
vcfs_dedup <- c(MUT_vcfs_dedup, MA_vcfs_dedup)
dimensions_after_dedup <- as.vector(sapply(vcfs_dedup,length))
barplot(dimensions_after_dedup, main="Number of substitutions after deduplication")

# Cleanup
rm(allVcfs, vcfs_dedup, MUT_allVcfs, duplicates, vcfs_filter12, vcfs_filter1, MUT_vcfs_dedup, MA_vcfs_dedup, filtered_indels)

## 2.3 Filtering SVs

# SVs are filtered and deduplicated in the preparation/classification step.
# run the script SV_filtering.R
load('filtered_SV.RData')

## 3. SPECTRA

# Working on: indels_dedup, vcfs_dedup, delly.filtcounts

# Categories:
#   - Substituions: single base subs in 96 contexts, DNV (NN > NN), MNV;
#   - Indels: Deletions - 1bp in/not in  homopolymers, 2-5bp in/not in repetitive region, 5-50, 50-400; Insertions - 1bp in/not in  homopolymers, 2-5bp in/not in repetitive region, 5-50, 50-400; DeletionsInsertions (comploex indels) - 1-50bp, 50-400bp;
#   - SVs: TD, DEL, INV, COMPLEX, TRSL, INTCHR, FOLDBACK

# adjust chromosome names, remove mitochondrial variants
for (i in 1:length(vcfs_dedup)) {
  if (length(vcfs_dedup[[i]])==0 || is.na(vcfs_dedup[[i]])) next
  seqlevels(vcfs_dedup[[i]]) <- seqnames(get(worm_ref_genome))
  #seqlevels(vcfs_dedup[[i]]) <- c("chrM","chrIV","chrIII","chrX","chrI","chrV","chrII")
}

vcfs_dedup <- sapply(vcfs_dedup, function(vcf) {
  if (length(vcf)==0 || is.na(vcf)) return(vcf)
  return(vcf[as.character(seqnames(vcf)) != "chrM"])
})

types.full <- paste(rep(rep(c('A','C','G','T'), each = 4), 6), '[', rep(c('C','T'), each = 48), '>', rep(c('A','G','T','A','C','G'), each=16), ']', rep(c('A','C','G','T'), 24), sep='')

spectrum <- list()
for (worm in names(CD2Mutant)){
  
  # substitutions
  ind <- match(worm,names(vcfs_dedup))
  sub.counts <- rep(NA,96)
  names(sub.counts) <- types.full
  dnv.counts <- rep(NA,2)
  names(dnv.counts) <- c('DNV','MNV')
  if (length(vcfs_dedup[[ind]])==0) { sub.counts <- vector("numeric",96); dnv.counts <- vector("numeric",2) }
  else if (!is.na(vcfs_dedup[[ind]])) {
    vcf <- vcfs_dedup[[ind]]
    tmp <- isMNV(vcf)
    if (sum(!tmp)==0)
      sub.counts <- vector("numeric",96)
    else {
      type_context = as.character(tncToPyrimidine(getTrinucleotideSubs(vcf[!tmp], worm_ref_genome)))
      counts <- table(type_context)
      sub.counts[names(counts)] <- counts
      sub.counts[is.na(sub.counts)] <- 0
    }
    if (sum(tmp)==0)
      #      dnv.counts <- vector("numeric",79)
      dnv.counts <- vector("numeric",2)
    else {
      vcf <- vcf[tmp]
      d <- diff(start(vcf))
      tmp_dif <- rle(d)
      dnv.counts[1] <- sum(tmp_dif$lengths==1 & tmp_dif$values==1)
      dnv.counts[2] <- length(which(tmp_dif$lengths>1 & tmp_dif$values==1))
    }
  }
  
  
  # repeat will mean repetitive sequence with >3 units of length up to 6
  # indels
  indel.counts <- rep(NA,14)
  ind <- match(worm,names(pindels_final))
  names(indel.counts) <- c("D.1.rep","D.1.nonrep",
                           "D.2.5.rep","D.2.5.nonrep",
                           "D.5.50","D.50.400",
                           "DI.small","DI.large",
                           "I.1.rep","I.1.nonrep",
                           "I.2.5.rep","I.2.5.nonrep",
                           "I.5.50","I.50.400")
  if (length(pindels_final[[ind]])==0) indel.counts <- vector("numeric",14)
  else if (!is.na(pindels_final[[ind]])) {
    vcf <- pindels_final[[ind]]
    tmp <- is.in.repeat(vcf)
    indel.counts["D.1.rep"] <- length(which(info(vcf)$PC=="D" &  info(vcf)$LEN==1 & tmp))
    indel.counts["D.1.nonrep"] <- length(which(info(vcf)$PC=="D" &  info(vcf)$LEN==1 & !tmp))
    indel.counts["D.2.5.rep"] <- length(which(info(vcf)$PC=="D" &  info(vcf)$LEN>1 & info(vcf)$LEN<6 & tmp))
    indel.counts["D.2.5.nonrep"] <- length(which(info(vcf)$PC=="D" &  info(vcf)$LEN>1 & info(vcf)$LEN<6 & !tmp))
    indel.counts["D.5.50"] <- length(which(info(vcf)$PC=="D" &  info(vcf)$LEN>5 & info(vcf)$LEN<51))
    indel.counts["D.50.400"] <- length(which(info(vcf)$PC=="D" &  info(vcf)$LEN>50))
    indel.counts["DI.small"] <- length(which(info(vcf)$PC=="DI" &  info(vcf)$LEN<51))
    indel.counts["DI.large"] <- length(which(info(vcf)$PC=="DI" &  info(vcf)$LEN>50))
    indel.counts["I.1.rep"] <- length(which(info(vcf)$PC=="I" &  info(vcf)$LEN==1 & tmp))
    indel.counts["I.1.nonrep"] <- length(which(info(vcf)$PC=="I" &  info(vcf)$LEN==1 & !tmp))
    indel.counts["I.2.5.rep"] <- length(which(info(vcf)$PC=="I" &  info(vcf)$LEN>1 & info(vcf)$LEN<6 & tmp))
    indel.counts["I.2.5.nonrep"] <- length(which(info(vcf)$PC=="I" &  info(vcf)$LEN>1 & info(vcf)$LEN<6 & !tmp))
    indel.counts["I.5.50"] <- length(which(info(vcf)$PC=="I" &  info(vcf)$LEN>5 & info(vcf)$LEN<51))
    indel.counts["I.50.400"] <- length(which(info(vcf)$PC=="I" &  info(vcf)$LEN>50))
    indel.counts[is.na(indel.counts)] <- 0
  }
  
  # rearrangements
  rearr.stuff = rep(NA,8)
  names(rearr.stuff) <- c("TD","DEL","INV","COMPLEX","TRSL","INTCHR","FOLDBACK","ALL")
  
  ind <- match(worm,names(delly.filtcounts))
  if (length(ind)>1) ind <- ind[1]
  if (!is.na(ind)) rearr.stuff <- c(delly.filtcounts[[ind]][1:4],delly.filtcounts[[ind]][5] + delly.filtcounts[[ind]][8],delly.filtcounts[[ind]][6:7],delly.filtcounts[[ind]][9])
  
  if (is.na(rearr.stuff[1])) print('no SV for this worm!')
  spectrum[[worm]] <- c(sub.counts,dnv.counts,indel.counts,rearr.stuff)
  
  print(worm)
}
spectrum <- do.call("rbind",spectrum)
spectrum <- as.data.frame(spectrum)

colnames(spectrum) <- c(types.full, 'DNV', 'MNV',
                        c("D.1.rep","D.1.nonrep",
                          "D.2.5.rep","D.2.5.nonrep",
                          "D.5.50","D.50.400",
                          "DI.small","DI.large",
                          "I.1.rep","I.1.nonrep",
                          "I.2.5.rep","I.2.5.nonrep",
                          "I.5.50","I.50.400"),
                        "TD","DEL","INV","COMPLEX","TRSL","INTCHR","FOLDBACK","ALL")

write.csv(spectrum, file = 'Spectrum.csv') # saved in Supplementary Table 1, sheet 3.

