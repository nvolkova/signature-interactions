########################################################
## MCMC sampling model for interaction effects in UV  ##
## N. Volkova, EBI, November 2018                     ##
########################################################
library(VariantAnnotation)
library(MASS)
source('../useful_functions.R')
library(ggplot2)
library(reshape2)
source('../plotting_functions.R')
library(greta)
library(openxlsx)


types.full <- paste(rep(rep(c('A','C','G','T'), each = 4), 6), '[', rep(c('C','T'), each = 48), '>', rep(c('A','G','T','A','C','G'), each=16), ']', rep(c('A','C','G','T'), 24), sep='')
indels <- c("D.1.5","D.5.50","D.50.400","DI.small","DI.large","I.1.5","I.5.50","I.50.400")

# Get the new signature set from here: https://www.synapse.org/#!Synapse:syn11967914
new_cancer_signatures <- read.csv('sigProfiler_exome_SBS_signatures.csv')


# Data
bigmat <- read.table('TCGA.caveman.matrix.dat',sep='\t')
rownames(bigmat) <- paste0('TCGA',substr(rownames(bigmat),5,nchar(rownames(bigmat))))
metadata <- read.table('TCGA_caveman_patients_to_files.dat', sep='\t', header=T)
skcm <- as.character(metadata$tumour[metadata$project=='SKCM'])

skcm <- skcm[which(sapply(skcm, function(x) cosine(as.numeric(uv[x,1:96]),new_cancer_signatures[,'SBS7a'] + new_cancer_signatures[,'SBS7b']))>0.8)]
uv <- bigmat[skcm,]
colnames(uv) <- c(types.full,indels)

# Expression - FPKM values from TCGA
metadata <- metadata[match(skcm, metadata$tumour),]
for (i in 6:ncol(metadata)) {
  metadata[,i] <- metadata[,i] / median(metadata[,i], na.rm = T)
}
thr <- 0.2

NER.core <- c('XPC','XPA','CUL5','ERCC1','ERCC2','ERCC4','ERCC5','ERCC6','CUL3')
TLS.core <- c('POLN','POLQ','REV3L','REV1','SHPRH','POLK','POLH','MAD2L2','POLI')

# Germline
PATHTOSNP='/path/to/germline/SNP/tables/per/cancertype'
# Germline mutations in TLS and NER: take a list of all TCGA germline mutations in these pathways, put them through Ensembl VEP
SNPs <- list()
tmp_NER <- read.table('germline_SKCM_NER_38.txt', sep = '\t', header = T)
tmp_TLS <- read.table('germline_TLS_mutations.txt', sep = '\t', header = T)
for (x in c(NER.core,TLS.core)) {
  if (x %in% NER.core)
    SNPs[[x]] <- unique(as.character(tmp_NER$Uploaded_variation[tmp_NER$IMPACT %in% c('HIGH') & grepl(x, tmp_NER$SYMBOL)]))
  else
    SNPs[[x]] <- unique(as.character(tmp_TLS$Uploaded_variation[tmp_TLS$IMPACT %in% c('HIGH') & grepl(x, tmp_TLS$SYMBOL)]))
}
SNPs <- SNPs[sapply(SNPs,length)>0]
cancername = 'SKCM'
snps_per_sample <- matrix(0,nrow = length(skcm), ncol = length(SNPs), dimnames = list(skcm, names(SNPs)))
for (genename in names(SNPs)) {
  f <- list.files(paste0(PATHTOSNP,'/',genename,'/',cancername))
  for (i in 1:length(skcm)) {
    file <- f[grep(skcm[i],f)]
    if (length(file)==0) next
    if (file.info(paste0(PATHTOSNP,'/',genename,'/',cancername,'/',file))$size==0) next
    vcf <- read.delim(paste0(PATHTOSNP,'/',genename,'/',cancername,'/',file), header = F)
    identifier <- paste0(vcf$V1,':',vcf$V2,'_',vcf$V4,'/',vcf$V5) 
    het <- substr(vcf$V10,1,3)[identifier %in% SNPs[[genename]]]
    if (length(het)>0) {
      snps_per_sample[skcm[i], genename] <- length(intersect(identifier, SNPs[[genename]]))
      if (sum(het == 1|1) >0)
        print(paste(c(skcm[i],genename,'homozygous SNP'),collapse=' '))
    }
  }
  print(genename)
}
colSums(snps_per_sample>0)


# Check defects per gene
X <- data.frame(NER = as.numeric(substr(skcm,1,12) %in% unlist(samples[paste0(NER.core,'_hom')]) |
                     rowSums(sapply(intersect(NER.core,colnames(metadata)), function(y) skcm %in% metadata$tumour[metadata[,y] < thr]))>0 |
                     snps_per_sample[,1]>0),
                TLS = as.numeric(substr(skcm,1,12) %in% unlist(samples[paste0(TLS.core,'_hom')]) |
                      rowSums(sapply(intersect(TLS.core,colnames(metadata)), function(y) skcm %in% metadata$tumour[metadata[,y] < thr]))>0 |
                     rowSums(snps_per_sample[,2:ncol(snps_per_sample)])>0))


# Data for sampling
skcm.data <- list(
  N = nrow(uv),
  R = 104,
  S = 1,
  K = ncol(X),
  X = as.matrix(X),
  y = as.matrix(uv),
  M = max(rowSums(uv))
)

# Sampling model
S <- variable(lower = 0, upper = 1, dim = c(skcm.data[['S']],skcm.data[['R']]))
S <- S / (greta::rowSums(S) %*% matrix(1,nrow = 1, ncol = skcm.data[['R']]))
E <- variable(lower = 0, upper = skcm.data[['M']], dim = c(skcm.data[['N']],skcm.data[['S']]))
sigma = 0.5
beta = normal(mean = 0, sd = sigma, dim = c(skcm.data[['K']],skcm.data[['R']]))
size = 50
mu1 = (E %*% S) * exp(skcm.data[['X']] %*% beta) / 
    ((exp(skcm.data[['X']] %*% beta) %*% t(S)) %*% matrix(1,nrow = 1, ncol=104))
prob = size/(size + mu1)
distribution(skcm.data[['y']]) = negative_binomial(size = size * matrix(1,nrow=skcm.data[['N']],ncol=skcm.data[['R']]),prob = prob)
                                                    
m <- model(S,E,beta)

# Sampling
draws <- greta::mcmc(m, n_samples = 500, warmup = 500, chains = 4)

# Visualize the draws
library(bayesplot)
mcmc_trace(draws[,grep('S',colnames(draws[[1]]))[sample(skcm.data[['S']]*c(1:104),12)]])
mcmc_trace(draws[,grep('beta',colnames(draws[[1]]))[sample(skcm.data[['K']]*c(1:104),12)]])

draws_all <- do.call('rbind',draws)

S_est <- matrix(colMeans(draws_all[,grep('S',colnames(draws_all), fixed = T)]), nrow = skcm.data[['S']], ncol = skcm.data[['R']])
S_low <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = skcm.data[['S']], ncol = skcm.data[['R']])
S_high <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = skcm.data[['S']], ncol = skcm.data[['R']])
rownames(S_est) = rownames(S_low) = rownames(S_high) = paste0('Sig',c(1:skcm.data[['S']]))
E_est <- matrix(colMeans(draws_all[,grep('E',colnames(draws_all))]), ncol = skcm.data[['S']])
beta_est <- matrix(colMeans(draws_all[,grep('beta',colnames(draws_all), fixed = T)]), nrow = skcm.data[['K']], ncol = skcm.data[['R']])
beta_low <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = skcm.data[['K']], ncol = skcm.data[['R']])
beta_high <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = skcm.data[['K']], ncol = skcm.data[['R']])
beta_var <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,var), nrow = skcm.data[['K']], ncol = skcm.data[['R']])

#plot_subindel_wb(t(S_est), CI = T, low = t(S_low), high = t(S_high))
plot_all_changes(plot_subindel_wb, S_est, beta_est)
boxplot(log10(E_est[,1]) ~ skcm.data[['X']][,1])
boxplot(log10(E_est[,1]) ~ skcm.data[['X']][,2])
interaction_effect_plot_human(beta_est[1,], lwd = 2, CI = T, low = beta_low[1,], high = beta_high[1,], at = c(-1,0,1))
interaction_effect_plot_human(beta_est[2,], lwd = 2, CI = T, low = beta_low[2,], high = beta_high[2,], at = c(-1,0,1))


pdf('UV_real_predicted.pdf',7,7)
mu_uv <- E_est %*% S_est * exp(skcm.data[['X']] %*% beta_est) / 
  ((exp(skcm.data[['X']] %*% beta_est) %*% t(S_est)) %*% matrix(1,nrow = 1, ncol=105))
plot(as.vector(skcm.data[['y']]),
     as.vector(mu_uv), pch = 16, xlab = 'Observed', ylab = 'Predicted', bty = 'n')
abline(a=0,b=1,col='red',lty=2)
boxplot(log10(E_est[,1]) ~ skcm.data[['X']][,1],bty = 'n', xaxt = 'n', frame = F, yaxt = 'n')
axis(side = 1, at = c(1,2), labels = c('Normal','NER def.'), tick = T, col = 'white', col.ticks = 'black')
axis(side = 2, at = c(2,3,4), labels = c(100,1000,10000), las = 2, tick = T, col = 'white', col.ticks = 'black')
dev.off()

pdf('UV_LFC.pdf',12,4)
interaction_effect_plot_human(beta_est[1,], lwd = 2, CI = T,
                              low = beta_low[1,],
                              high = beta_high[1,], at = c(-1,0,1),
                              labels = c('<0.1',1,10),
                              plot_main = 'UV change w.r.t. NER mutations')
dev.off()

S_var <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,var), nrow = skcm.data[['S']], ncol = skcm.data[['R']])
tmp <- data.frame(S_est[1,],
                  S_est[1,] * exp(beta_est[1,]))
var_tmp <- exp(2*beta_est[1,]) * S_var[1,] + S_est[1,]**2 * exp(2*beta_est[1,]) * beta_var[1,] # variance of S_est[1,] * exp(beta_est[2,])
tmp_low <- data.frame(S_low[1,],
                      S_est[1,] * exp(beta_est[1,]) - 1.96 * sqrt(var_tmp))
tmp_up <- data.frame(S_high[1,],
                     S_est[1,] * exp(beta_est[1,]) + 1.96 * sqrt(var_tmp))
colnames(tmp) = colnames(tmp_low) = colnames(tmp_up) <- c('UV','UV + muts')
plot_subindel_wb(tmp,CI = T,low = tmp_low, high = tmp_up, ymax = max(tmp_up), norm = F) + theme(panel.grid = element_blank())


AGE = ids$donor_age_at_diagnosis[match(substr(rownames(skcm.data[['y']]),1,12), ids$submitted_donor_id)]
no_muts_average <- mean(rowSums(skcm.data[['y']])[skcm.data[['X']][,1]==0] / AGE[skcm.data[['X']][,1]==0], na.rm = T)
tmp$UV <- tmp$UV * no_muts_average
tmp_low$UV <- tmp_low$UV * no_muts_average
tmp_up$UV <- tmp_up$UV * no_muts_average

muts_average <- mean(rowSums(skcm.data[['y']])[skcm.data[['X']][,1]==1] / AGE[skcm.data[['X']][,1]==1], na.rm = T)
tmp_low$`UV + muts` <- tmp_low$`UV + muts` * muts_average / sum(tmp$`UV + muts`)
tmp_up$`UV + muts` <- tmp_up$`UV + muts` * muts_average / sum(tmp$`UV + muts`)
tmp$`UV + muts` <- tmp$`UV + muts` * muts_average / sum(tmp$`UV + muts`)

interaction_effect_plot_human(log10(tmp$`UV + muts` / tmp$UV),
                              CI = F, log = T, at = c(-1,0,1))

pdf('UV_signatures.pdf',12,5)
plot_subindel_wb(tmp,CI = T,low = tmp_low, high = tmp_up, ymax = max(tmp_up), norm = F) + theme(panel.grid = element_blank())
dev.off()

sd_uv_lfc <- sqrt(1/tmp$`UV + muts`**2 * var_tmp * muts_average**2 + 
                    1/tmp$UV**2  * S_var[1,] * no_muts_average**2)

UV_LFC <- mean(tmp$`UV + muts` / tmp$UV) # 1.814657
UV_LFC_sd <- sqrt(sum(sd_uv_lfc**2) / 104**2) # 0.04520428

pdf('UV_LFC.pdf',12,4)
interaction_effect_plot_human(log(tmp$`UV + muts` / tmp$UV), lwd = 2, CI = T,
                              low = log(tmp$`UV + muts` / tmp$UV) - 1.96*sd_uv_lfc,
                              high = log(tmp$`UV + muts` / tmp$UV) + 1.96*sd_uv_lfc, at = c(-1,0,log10(2),1),
                              labels = c('<0.1',1,2,10),
                              plot_main = 'UV change w.r.t. NER mutations')
dev.off()





