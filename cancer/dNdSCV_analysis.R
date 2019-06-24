# dNdS analysis in a nutshell

# Upload the package
library(dndscv)
library(openxlsx)
library(ggplot2)
library(reshape2)
library(beeswarm)
PATHTOCANCERDATA='/path/to/tables/with/mutations/per/cancer/type'
# https://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0204912.s006

# Undercovered genes: from Wang et al 2018
download.file(url = 'https://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0204912.s006', 
              destfile = 'undercovered_genes.xlsx', method, quiet = FALSE, mode = "w",
              cacheOK = TRUE,
              extra = getOption("download.file.extra"))
a1 <- xlsx::read.xlsx('undercovered_genes.xlsx', sheetIndex = 5)
a2 <- xlsx::read.xlsx('undercovered_genes.xlsx', sheetIndex = 6)
a3 <- xlsx::read.xlsx('undercovered_genes.xlsx', sheetIndex = 7)
undercover <- unique(c(as.character(a1$Gene.Name[a1$Fraction.of.bases.with.probe.coverage<0.7]),
                       as.character(a2$Gene.Name[a2$Fraction.of.bases.with.probe.coverage<0.7]),
                       as.character(a3$Gene.Name[a3$Fraction.of.bases.with.probe.coverage<0.7])))

# Get list of all the genes
all_genes <- read.table('nb2/dNdS/Reference_list_genes_Human_CDS_dNdS_new.txt',header=F,sep='\t',stringsAsFactors = F)[,1]
# List of repair genes
repair_genes <- openxlsx::read.xlsx('Supplementary Table 4. TCGA samples with mutations in DNA pathways.xlsx',sheet = 1)

##############################################

# Read in the pancancer mutation set

mutations <- read.delim(paste0(PATHTOCANCERDATA,'/all_TCGA_mutations_cavemanpindel.txt'), header=T, skip=66)
mutations <- mutations[,c('Sample','Chrom','Pos','Alt','Ref','PM.Tum','Type','Effect','Gene')]
mutations <- mutations[mutations$Type == 'Sub',]
indels <- read.delim(paste0(PATHTOCANCERDATA,'/all_TCGA_mutations_pindel.vcf'), header = T, skip = 65)
indels$PM.Tum <- (indels$PU.Tum + indels$NU.Tum) / (indels$PR.Tum + indels$NR.Tum)
indels <- indels[,c('Sample','Chrom','Pos','Alt','Ref','PM.Tum','Type','Effect','Gene')]
mutations <- rbind(mutations, indels)
mutations <- mutations[mutations$Gene %in% all_genes,]
ref_muts <- data.frame(
  sampleID = mutations$Sample,
  chr = mutations$Chrom,
  pos = mutations$Pos,
  ref = mutations$Ref,
  mut = mutations$Alt,
  stringsAsFactors = FALSE
)
distances <- diff(ref_muts$pos)
bad <- sort(c(which(distances == 1), which(distances == 1) + 1))
ref_muts <- ref_muts[-bad,]

###############################################

# Single gene analysis under 1000 coding mutations

res.1000 <- dndscv(ref_muts, max_coding_muts_per_sample = 1000, outmats = T,
                   max_muts_per_gene_per_sample = 5, # restore some mutations in the set
                   constrain_wnon_wspl = F) # analyze splice site mutations and nonsense separately
confintervals <- geneci(res.1000, gene_list = repair_genes$Gene)

df.1000.single.gene <- confintervals

# Missense
mis_sds <- ((log(confintervals$mis_high) - log(confintervals$mis_mle)) / 1.96)
scores <- ((log(confintervals$mis_mle))**2 / (mis_sds)**2)
df.1000.single.gene$pvalue.mis <- 1 - pchisq(scores, df = 1)
df.1000.single.gene$qvalue.mis <- p.adjust(1 - pchisq(scores, df = 1),method='BH')
confintervals$gene[which(p.adjust(1 - pchisq(scores, df = 1),method='BH') < 0.1)]
# Nonsense
non_sds <- ((log(confintervals$non_high) - log(confintervals$non_mle)) / 1.96)
scores <- ((log(confintervals$non_mle))**2 / (non_sds)**2)
df.1000.single.gene$pvalue.non <- 1 - pchisq(scores, df = 1)
df.1000.single.gene$qvalue.non <- p.adjust(1 - pchisq(scores, df = 1),method='BH')
confintervals$gene[which(p.adjust(1 - pchisq(scores, df = 1),method='BH') < 0.1)]

###############################################

# Pathway analysis

dndslist <-list()
for (cancer in cancerlist) {
  for (pw in pwlist) {
    
    mutations <- read.delim(paste0(PATHTOCANCERDATA,'/',cancer,'_cavemanpindel.txt'), header=T, skip=66)
    mutations <- mutations[,c('Sample','Chrom','Pos','Alt','Ref','PM.Tum','Type','Effect','Gene')]
    mutations <- mutations[mutations$Type == 'Sub',]
    indels <- read.delim(paste0(PATHTOCANCERDATA,'/',cancer,'_pindel.vcf'), header = T, skip = 65)
    indels$PM.Tum <- (indels$PU.Tum + indels$NU.Tum) / (indels$PR.Tum + indels$NR.Tum)
    indels <- indels[,c('Sample','Chrom','Pos','Alt','Ref','PM.Tum','Type','Effect','Gene')]
    mutations <- rbind(mutations, indels)
    mutations <- mutations[mutations$Gene %in% all_genes,]
    all_muts <- data.frame(
      sampleID = mutations$Sample,
      chr = mutations$Chrom,
      pos = mutations$Pos,
      ref = mutations$Ref,
      mut = mutations$Alt,
      stringsAsFactors = FALSE
    )
    distances <- diff(all_muts$pos)
    bad <- sort(c(which(distances == 1), which(distances == 1) + 1))
    all_muts <- all_muts[-bad,]
    
    if (nrow(all_muts) > 100) {
      dndslist[[paste(pw,cancer,sep='_')]] <- dnds(all_muts, max_coding_muts_per_sample = 1000,
                                                   gene_list = intersect(all_genes, repair_genes$Gene[repair_genes$Pathway == pw & repair_genes$CORE]))$globaldnds
    }
  }
}

# Missense
missdndslist <- lapply(dndslist, function(l) l[1,,drop = F])
missdndslist <- do.call('rbind',missdndslist)
mis_sds <- (log(missdndslist$cihigh) - log(missdndslist$mle)) / 1.96
mis_scores <- (log(missdndslist$mle) / mis_sds)**2
names(mis_scores) <- rownames(missdndslist)
which(p.adjust(1 - pchisq(mis_scores, df = 1),method='BH') < 0.05)
# Nonsense
nonedndslist <- lapply(dndslist, function(l) l[2,,drop = F])
nonedndslist <- do.call('rbind',nonedndslist)
non_sds <- (log(nonedndslist$cihigh) - log(nonedndslist$mle)) / 1.96
scores <- (log(nonedndslist$mle) / non_sds)**2
names(scores) <- rownames(nonedndslist)
which(p.adjust(1 - pchisq(scores, df = 1),method='BH') < 0.05)

library(reshape2)
df <- data.frame(Name = names(dndslist),
                 wmis = missdndslist$mle,
                 wmis_low = missdndslist$cilow,
                 wmis_high = missdndslist$cihigh,
                 wnon = nonedndslist$mle,
                 wnon_low = nonedndslist$cilow,
                 wnon_high = nonedndslist$cihigh,
                 pvalue.mis = 1 - pchisq(mis_scores, df = 1),
                 qvalue.mis = p.adjust(1 - pchisq(mis_scores, df = 1),method='BH'),
                 pvalue.non = 1 - pchisq(scores, df = 1),
                 qvalue.non = p.adjust(1 - pchisq(scores, df = 1),method='BH')
                 )

###########################################################3
# Gather together for plotting

df2 <- melt(df[,c(1,2,5)])
df2$qvalue <- NA
df2$qvalue[df2$variable == 'wmis'] <- p.adjust(1 - pchisq(mis_scores, df = 1),method='BH')[match(df2$Name[df2$variable == 'wmis'], names(dndslist))]
df2$qvalue[df2$variable == 'wnon'] <- p.adjust(1 - pchisq(scores, df = 1),method='BH')[match(df2$Name[df2$variable == 'wnon'], names(dndslist))]
df2$value <- log10(df2$value)
df2$Name <- as.character(df2$Name)
df2$variable <- as.character(df2$variable)

# Add per-gene results
df3 <- melt(df.1000.single.gene[,c(1:3)], id.vars = c(1))
df3$qvalue[df3$variable == 'mis_mle'] <- df.1000.single.gene$qvalue.mis[match(df.1000.single.gene$gene,df3$gene[df3$variable == 'mis_mle'])]
df3$qvalue[df3$variable == 'non_mle'] <- df.1000.single.gene$qvalue.non[match(df.1000.single.gene$gene,df3$gene[df3$variable == 'non_mle'])]
df3$gene <- as.character(df3$gene)
df3$variable <- as.character(df3$variable)
colnames(df3)[1] <- 'Name'
df3$value <- log10(df3$value)
df2 <- rbind(df2,df3)
df2$variable <- factor(df2$variable)


# Plot the results
par(mar = c(8,4,2,4))
library(beeswarm)
fff <- beeswarm(value ~ variable, data = df2,
                corral = 'wrap', col = 'lightgrey',
                corralWidth = 0.5, cex = 0.7, xlim = c(0.5,5),
                las = 2, lty = 0,
                xlab = '', ylab = '', ylim = c(-1,2), 
                pch =  16, yaxt = 'n', xaxt = 'n')
axis(side = 1, at = c(1:4), labels = c('missense','nonsense','missense \nper pathway','nonsense \nper pathway'),
     col = 'white',cex.axis = 1.2, las=2)
axis(side = 2, las = 2, at = c(log10(0.1*c(2:9)),log10(c(2:9)),log10(10*c(2:9))), labels = rep('',24),  cex.axis = 1.2, tck = -0.01)
axis(side = 2, las = 2, at = c(-1,0,1,2), labels = c(0.1, 1, 10,100),  cex.axis = 1.2)
abline(h = 0, lty = 2, col = 'red')

mis_plotting <- fff[grep('mis_mle',fff$x.orig),]
mis_plotting <- mis_plotting[mis_plotting$'y.orig' %in% df2$value[df2$qvalue < 0.1 & df2$variable == 'mis_mle' ],]#& df2$value > 0],]
points(x = mis_plotting$x, y = mis_plotting$y, col = 'gray20', pch = 16)
text(x = 1.5, y = seq(2,0.05,length.out = sum(df2$qvalue < 0.1 & df2$variable == 'mis_mle' & df2$value > 0, na.rm = T)), cex = 0.7,font=3,
     labels = df2[df2$qvalue < 0.1 & df2$variable == 'mis_mle' & df2$value > 0,'Name'][order(df2[df2$qvalue < 0.1 & 
                                                                                                   df2$variable == 'mis_mle' & df2$value > 0,'value'], decreasing = T)])
#text(x = 1.5, y = seq(-0.2,-0.5,length.out = sum(df2$qvalue < 0.1 & df2$variable == 'wmis' & df2$value < 0, na.rm = T)),
#     labels = df2[df2$qvalue < 0.1 & df2$variable == 'wmis' & df2$value < 0,'Name'][order(df2[df2$qvalue < 0.1 & df2$variable == 'wmis' & df2$value < 0,'value'], decreasing = T)])
# nothing on the negative side

mis_plotting <- fff[grep('non_mle',fff$x.orig),]
mis_plotting <- mis_plotting[mis_plotting$'y.orig' %in% df2$value[df2$qvalue < 0.1 & df2$variable == 'non_mle'],]
points(x = mis_plotting$x, y = mis_plotting$y, col = 'gray20', pch = 16)
text(x = 2.5, y = seq(2,0.05,length.out = sum(df2$qvalue < 0.1 & df2$variable == 'non_mle' & df2$value > 0, na.rm = T)),cex = 0.7,font=3,
     labels = df2[df2$qvalue < 0.1 & df2$variable == 'non_mle' & df2$value > 0,'Name'][order(df2[df2$qvalue < 0.1 & df2$variable == 'non_mle' & df2$value > 0,'value'], decreasing = T)])
#text(x = 2.5, y = seq(-0.2,-0.05,length.out = sum(df2$qvalue < 0.1 & df2$variable == 'non_mle' & df2$value < 0, na.rm = T)),cex = 0.7,
#     labels = df2[df2$qvalue < 0.1 & df2$variable == 'non_mle' & df2$value < 0,'Name'][order(df2[df2$qvalue < 0.1 & df2$variable == 'non_mle' & df2$value < 0,'value'], decreasing = T)])
# nothing on the negative side

mis_plotting <- fff[grep('wmis',fff$x.orig),]
mis_plotting <- mis_plotting[mis_plotting$'y.orig' %in% df2$value[df2$qvalue < 0.1 & df2$variable == 'wmis' ],]#& df2$value > 0],]
points(x = mis_plotting$x, y = mis_plotting$y, col = 'gray20', pch = 16)
text(x = 3.5, y = seq(2,0.05,length.out = sum(df2$qvalue < 0.1 & df2$variable == 'wmis' & df2$value > 0, na.rm = T)), cex = 0.7,
     labels = df2[df2$qvalue < 0.1 & df2$variable == 'wmis' & df2$value > 0,'Name'][order(df2[df2$qvalue < 0.1 & df2$variable == 'wmis' & df2$value > 0,'value'], decreasing = T)])
#text(x = 3.5, y = seq(-0.2,-0.5,length.out = sum(df2$qvalue < 0.1 & df2$variable == 'wmis' & df2$value < 0, na.rm = T)),
#     labels = df2[df2$qvalue < 0.1 & df2$variable == 'wmis' & df2$value < 0,'Name'][order(df2[df2$qvalue < 0.1 & df2$variable == 'wmis' & df2$value < 0,'value'], decreasing = T)])
# nothing on the negative side

mis_plotting <- fff[grep('wnon',fff$x.orig),]
mis_plotting <- mis_plotting[mis_plotting$'y.orig' %in% df2$value[df2$qvalue < 0.1 & df2$variable == 'wnon'],]
points(x = mis_plotting$x, y = mis_plotting$y, col = 'gray20', pch = 16)
text(x = 4.5, y = seq(2,0.05,length.out = sum(df2$qvalue < 0.1 & df2$variable == 'wnon' & df2$value > 0, na.rm = T)),cex = 0.7,
     labels = df2[df2$qvalue < 0.1 & df2$variable == 'wnon' & df2$value > 0,'Name'][order(df2[df2$qvalue < 0.1 & df2$variable == 'wnon' & df2$value > 0,'value'], decreasing = T)])
#text(x = 4.5, y = seq(-0.2,-0.05,length.out = sum(df2$qvalue < 0.1 & df2$variable == 'wnon' & df2$value < 0, na.rm = T)),cex = 0.7,
#     labels = df2[df2$qvalue < 0.1 & df2$variable == 'wnon' & df2$value < 0,'Name'][order(df2[df2$qvalue < 0.1 & df2$variable == 'wnon' & df2$value < 0,'value'], decreasing = T)])
# nothing on the negative side

# save df2 to a sheet of SupTable5

openxlsx::write.xlsx(list(df.1000.single.gene, df), file = 'Supplementary Table 5. dNdS values for DNA repair genes and pathways.xlsx',
                     colNames = T, rowNames = F, sheetName = c('Per-gene dNdS', 'Per pathway dNdS'))


