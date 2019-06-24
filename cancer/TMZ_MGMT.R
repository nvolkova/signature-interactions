########################################################
## MCMC sampling model for interaction effects in GBM ##
## N. Volkova, EBI, November 2018                     ##
########################################################

library(greta)
library(openxlsx)
source('..//useful_functions.R')

cosine <- function(x,y) {
  x %*% y / sqrt(sum(x**2)) / sqrt(sum(y**2))
}
types.full <- paste(rep(rep(c('A','C','G','T'), each = 4), 6), '[', rep(c('C','T'), each = 48), '>', rep(c('A','G','T','A','C','G'), each=16), ']', rep(c('A','C','G','T'), 24), sep='')
indels <- c("D.1.5","D.5.50","D.50.400","DI.small","DI.large","I.1.5","I.5.50","I.50.400")

# From Kim et al 2015
# https://genome.cshlp.org/content/early/2015/02/03/gr.180612.114
# Download supplements
maf <- openxlsx::read.xlsx('supp_gr.180612.114_Supplemental_Tables_S1-S6.xlsx', sheet = 6, colNames = TRUE)
info <- openxlsx::read.xlsx('supp_gr.180612.114_Supplemental_Tables_S1-S6.xlsx', sheet = 1, colNames = TRUE)
info <- info[-c(1,25),]
head(maf)
colnames(maf)[c(7,8,13,14)] <- c('Chrom','Pos','Ref','Alt')

ref_genome='BSgenome.Hsapiens.UCSC.hg19'
library(BSgenome.Hsapiens.UCSC.hg19)

TMZmat <- lapply(unique(maf$Patient), function(x) {
  inds <- grep(x,maf$Patient)
  tab <- maf[inds,][maf$Tumor_Type[inds]=='recurrence',]
  mutations = numeric(104)
  names(mutations)=c(types.full, indels)
  type_context = as.character(tncToPyrimidine(getTrinucleotideSubs_human_table(tab)))
  counts <- table(type_context)
  mutations[names(counts)] <- counts
  return(t(mutations))
})
names(TMZmat) <- unique(maf$Patient)
TMZmat <- do.call('rbind',TMZmat)
rownames(TMZmat) <- unique(maf$Patient)
colnames(TMZmat)[1:96] <- types.full

plot_subindel_wb(t(TMZmat), norm=F)

non_TMZ <- c("TCGA-06-0152","TCGA-06-0210","TCGA-06-0211","TCGA-06-0221","TCGA-14-0736","TCGA-14-1034")
TMZmat <- TMZmat[-match(non_TMZ,rownames(TMZmat)),]
mgmt_status <- sapply(rownames(TMZmat), function(x) info$MGMT[grep(x,info$Patient)])
to.plot <- data.frame(MGMT = ifelse(grepl('MET',mgmt_status),'MET','non-MET'))
to.plot$number <- rowSums(TMZmat)
to.plot$similarity <- sapply(1:nrow(TMZmat), function(i) cosine(TMZmat[i,1:96], (beta_M_greta_full[,'EMS'] * exp(beta_I_greta_full[,'agt.1.EMS']))[1:96]))

pdf('~/plot_for_TMZ.pdf',10,7)
boxplot(similarity ~ MGMT, data=to.plot, xlab = 'MGMT Methylation Status',
        ylab = 'Cosine similarity to agt-1:EMS profile',
        main = 'TMZ treated GBM samples')
boxplot(number ~ MGMT, data=to.plot, xlab = 'MGMT Methylation Status',
        ylab = 'Number of mutations', main = 'TMZ treated GBM samples')
barplot(c(rowSums(TMZmat[names(mgmt_status)[grepl('MET',mgmt_status)],]), rowSums(TMZmat[names(mgmt_status)[!grepl('MET',mgmt_status)],])), col=c(rep('darkred',5),rep('darkblue',11)),
        xlab = 'Sample', main='Number of mutations')
dev.off()

fit <- glm(number ~ MGMT, data = to.plot[c(2:nrow(to.plot),1),], family=stats::poisson(link = 'identity'))
summary(fit)

library(beeswarm)
to.plot$number <- log10(to.plot$number)
pdf('TMZ_boxplot.pdf',3,4)
f <- beeswarm(to.plot$number ~ to.plot$MGMT, pch = 16, col = 'gray74', 
              bty = 'n', ylab = 'Total burden', xaxt = 'n', xlab = 'MGMT status',
              corral = 'wrap',
              method='hex', yaxt = 'n',
              corralWidth = 0.3)
axis(side = 1, at = c(1,2), labels = c('MET (6)','non-MET (11)'), col = 'white')
axis(side = 2, at = c(1,2,3,4), labels = c(10,100,1000,10000), las = 2)
axis(side = 2, at = c(log10(c(20,30,40,50,60,70,80,90)),
                      log10(c(200,300,400,50,600,700,800,900)),
                      log10(c(2000,3000,4000,5000,6000,7000,8000,9000))), labels = rep('',24), las = 2)
boxplot(to.plot$number ~ to.plot$MGMT, frame = F, add = T,
        col=rgb(1,1,1,alpha=0.2), outline = F, xaxt = 'n', yaxt = 'n', staplewex=0, lty = 1,boxwex=.5)
axis(side=3, tcl=0.25, at=1:2, labels=NA)
mtext(side=3, paste0("P = ",signif(summary(fit)$coef[2,4],1)))
dev.off()

tmz.data <- list(
  y = TMZmat,
  N = nrow(TMZmat),
  R = ncol(TMZmat),
  K = 1,
  X = matrix(as.numeric(to.plot$MGMT=='MET')),
  M = max(rowSums(TMZmat)),
  S = 2
)

S <- variable(lower = 0, upper = 1, dim = c(tmz.data[['S']],tmz.data[['R']]))
S <- S / (greta::rowSums(S) %*% matrix(1,nrow = 1, ncol = tmz.data[['R']]))
E <- variable(lower = 0, upper = tmz.data[['M']], dim = c(tmz.data[['N']],tmz.data[['S']]))

sigma = variable(lower = 0, upper = 5)
beta = normal(mean = 0, sd = sigma, dim = c(tmz.data[['K']],tmz.data[['R']]))

size = 50

mu1 = E[,-(tmz.data[['S']]),drop=F] %*% S[-(tmz.data[['S']]),,drop=F] + 
  (E[,tmz.data[['S']],drop=F] %*% S[tmz.data[['S']],,drop = F]) * exp(tmz.data[['X']] %*% beta) / 
  ((exp(tmz.data[['X']] %*% beta) %*% t(S[tmz.data[['S']],,drop=F])) %*% matrix(1,nrow = 1, ncol=104))
prob = size/(size + mu1)
distribution(tmz.data[['y']]) = negative_binomial(size = matrix(1, nrow = tmz.data[['N']], ncol = tmz.data[['R']]) * size, prob = prob)
# defining the model
m <- model(S,E,beta,sigma)

# sampling
draws <- mcmc(m, n_samples = 1000, warmup = 1500, verbose = F, chains = 4)

# Visualize the draws
library(bayesplot)
mcmc_trace(draws[,grep('S',colnames(draws[[1]]))[sample(1:(2*104),9)]])
mcmc_trace(draws[,grep('beta',colnames(draws[[1]]))[1:9]])
mcmc_trace(draws[,grep('sigma',colnames(draws[[1]])),drop=F])

draws_all <- do.call('rbind',draws)

S_est <- matrix(colMeans(draws_all[,grep('S',colnames(draws_all), fixed = T)]), nrow = tmz.data[['S']], ncol = tmz.data[['R']])
S_low <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = tmz.data[['S']], ncol = tmz.data[['R']])
S_high <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = tmz.data[['S']], ncol = tmz.data[['R']])
rownames(S_est) = rownames(S_low) = rownames(S_high) = c('Sig1','Sig2')
E_est <- matrix(colMeans(draws_all[,grep('E',colnames(draws_all))]), ncol = tmz.data[['S']])
beta_est <- matrix(colMeans(draws_all[,grep('beta',colnames(draws_all), fixed = T)]), nrow = tmz.data[['K']], ncol = tmz.data[['R']])
beta_low <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = tmz.data[['K']], ncol = tmz.data[['R']])
beta_high <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = tmz.data[['K']], ncol = tmz.data[['R']])
beta_var <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,var), nrow = tmz.data[['K']], ncol = tmz.data[['R']])

plot_subindel_wb(t(S_est), CI = T, low = t(S_low), high = t(S_high), ymax = 0.2)
plot_subindel_wb(cbind(t(S_est), S_est[2,] * exp(beta_est[1,])))

pdf('temozolomide_lfc.pdf',12,4)
interaction_effect_plot_human(c(beta_est[1,],0), lwd = 2, CI = T, low = c(beta_low[1,],0), high = c(beta_high[1,],0), at = c(-1,0,1,2,3))
dev.off()

pdf('tmz_quality.pdf',7,7)
mu_tmz <- E_est[,-(tmz.data[['S']]),drop=F] %*% S_est[-(tmz.data[['S']]),,drop=F] + 
  (E_est[,tmz.data[['S']],drop=F] %*% S_est[tmz.data[['S']],,drop = F]) * exp(tmz.data[['X']] %*% beta_est) / 
  ((exp(tmz.data[['X']] %*% beta_est) %*% t(S_est[tmz.data[['S']],,drop=F])) %*% matrix(1,nrow = 1, ncol=104))

plot(as.vector(tmz.data[['y']]),
     as.vector(mu_tmz), pch = 16)
abline(a=0,b=1,col='red',lty=2)
dev.off()

pdf('temozolomide_all_signatures.pdf',12,5)
tmp <- data.frame(t(S_est), S_est[2,] * exp(beta_est[1,]))
colnames(tmp) <- c('Sig1','Sig2','Sig2 + MGMT')
plot_subindel_wb(tmp)
dev.off()
pdf('temozolomide_some_signatures.pdf',12,5)
S_var <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,var), nrow = tmz.data[['S']], ncol = tmz.data[['R']])
var_tmp <- exp(2*beta_est[1,]) * S_var[2,] + S_est[2,]**2 * exp(2*beta_est[1,]) * beta_var[1,]
tmp_low <- data.frame(S_low[2,], S_est[2,] * exp(beta_est[1,]) - 1.96 * sqrt(var_tmp) / sqrt(1000))
tmp_up <- data.frame(S_high[2,], S_est[2,] * exp(beta_est[1,]) + 1.96 * sqrt(var_tmp) / sqrt(1000))
tmp <- data.frame(S_est[2,], S_est[2,] * exp(beta_est[1,]))
colnames(tmp) = colnames(tmp_up) = colnames(tmp_low) <- c('Sig2','Sig2 + MGMT')
plot_subindel_wb(tmp,CI = T,low = tmp_low, high = tmp_up, ymax = 0.25) + theme(panel.grid = element_blank())
dev.off()
boxplot(E_est[,1] ~ tmz.data[['X']][,1])
boxplot(log10(E_est[,2]) ~ tmz.data[['X']][,1])

S_var <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,var), nrow = tmz.data[['S']], ncol = tmz.data[['R']])
var_tmp <- exp(2*beta_est[1,]) * S_var[2,] + S_est[2,]**2 * exp(2*beta_est[1,]) * beta_var[1,]
tmp <- data.frame(S_est[2,],
                  S_est[2,] * exp(beta_est[1,]) / sum(S_est[2,] * exp(beta_est[1,]) ))
tmp_low <- data.frame(S_low[2,],
                      tmp[,2] - 1.96 * sqrt(var_tmp) / sum(S_est[2,] * exp(beta_est[1,]) ))
tmp_up <- data.frame(S_high[2,],
                     tmp[,2] + 1.96 * sqrt(var_tmp) / sum(S_est[2,] * exp(beta_est[1,]) ))
tmp_low[tmp_low<0] <- 0
colnames(tmp) = colnames(tmp_low) = colnames(tmp_up) <- c('TMZ','TMZ + MGMT')
plot_subindel_wb(tmp, CI = T, low = tmp_low, high = tmp_up)

no_muts_average <- mean(E_est[,2][tmz.data[['X']]==0])
tmp$TMZ <- tmp$TMZ * no_muts_average
tmp_low$TMZ <- tmp_low$TMZ * no_muts_average
tmp_up$TMZ <- tmp_up$TMZ * no_muts_average

muts_average <- mean(E_est[,2][tmz.data[['X']]==1])
tmp_low$`TMZ + MGMT` <- tmp_low$`TMZ + MGMT` * muts_average 
tmp_up$`TMZ + MGMT` <- tmp_up$`TMZ + MGMT` * muts_average 
tmp$`TMZ + MGMT` <- tmp$`TMZ + MGMT` * muts_average 

interaction_effect_plot_human(c(log(tmp$`TMZ + MGMT` / tmp$TMZ),0),
                              CI = F, log = T, at = c(-1,0,1,2,3,4))

pdf('temozolomide_some_signatures.pdf',12,5)
plot_subindel_wb(tmp,CI = T,low = tmp_low, high = tmp_up, ymax =0.2) + theme(panel.grid = element_blank())
dev.off()

sd_tmz_lfc <- sqrt(1/tmp$`TMZ + MGMT`**2 * var_tmp + 
                    1/tmp$TMZ**2 * S_var[2,])

TMZ_LFC <- mean(tmp$`TMZ + MGMT` / tmp$TMZ) # 51.60243
TMZ_LFC_sd <- sqrt(sum(sd_tmz_lfc**2) / 104**2) # 0.9135637

pdf('temozolomide_lfc.pdf',12,4)
interaction_effect_plot_human(c(log(tmp$`TMZ + MGMT` / tmp$TMZ),0), lwd = 2, CI = T,
                              low = c(log(tmp$`TMZ + MGMT` / tmp$TMZ) - 1.96*sd_tmz_lfc,0),
                              high = c(log(tmp$`TMZ + MGMT` / tmp$TMZ) + 1.96*sd_tmz_lfc,0), at = c(-1,0,1,2,3,4),
                              #labels = c('<0.1',0.5,1,10),
                              plot_main = 'AA change w.r.t. TC-NER mutations')
dev.off()

