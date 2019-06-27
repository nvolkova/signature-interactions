#############################################################################
## MCMC sampling model for interaction effects in MMR-def cancers          ##
## N. Volkova, EBI, November 2018                                          ##
#############################################################################

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

bigmat <- read.table('TCGA.caveman.matrix.dat',sep='\t')
metadata <- read.table('TCGA_caveman_patients_to_files.dat', sep='\t', header=T)
rownames(bigmat) <- paste0('TCGA', substr(rownames(bigmat), 5,nchar(rownames(bigmat))))

# MSI status for TCGA from TCGA Clinical Explorer
files <- list.files('~/Downloads/MSI_TCGA/')
cancer_types <- as.character(sapply(files, function(x) unlist(strsplit(x, split='[_]'))[1]))
status <- lapply(files, function(x) read.table(paste0('~/Downloads/MSI_TCGA/',x), sep='\t', header=T))
msi_status = data.frame(sample = do.call('c', lapply(status, function(x) as.character(x$SampleCode))),
                        cancer = rep(cancer_types, times = sapply(status, function(x) nrow(x))),
                        MSI = do.call('c', lapply(status, function(x) as.character(x$MSIstatus))))
# Only COADREAD, ESCA, STAD, UCEC :(

mmr_data <- list()
mmr_data[['y']] <- do.call('rbind', lapply(msi_status$sample[msi_status$MSI=='MSI-H'],
                                           function(x) bigmat[grep(x,rownames(bigmat)),]))
mmr_data[['N']] <- nrow(mmr_data[['y']])
mmr_data[['R']] <- ncol(mmr_data[['y']])
mmr_data[['K']] <- 3
mmr_data[['X']] <- matrix(0, nrow = mmr_data[['N']], ncol = 3)
mmr_data[['X']][,1] <- as.numeric(metadata$project[match(rownames(mmr_data[['y']]), metadata$tumour)] %in% c('COAD','READ'))
mmr_data[['X']][,2] <- as.numeric(metadata$project[match(rownames(mmr_data[['y']]), metadata$tumour)] == 'STAD')
mmr_data[['X']][,3] <- as.numeric(metadata$project[match(rownames(mmr_data[['y']]), metadata$tumour)] == 'UCEC')
#mmr_data[['X']][,4] <- as.numeric(metadata$project[match(rownames(mmr_data[['y']]), metadata$tumour)] == 'ESCA')
#mmr_data[['X']][,5] <- as.numeric(metadata$project[match(rownames(mmr_data[['y']]), metadata$tumour)] == 'UCS')
mmr_data[['S']] <- 1
mmr_data[['M']] <- max(rowSums(mmr_data[['y']]))

# Greta model

S <- variable(lower = 0, upper = 1, dim = c(mmr_data[['S']],mmr_data[['R']]))
S <- S / (greta::rowSums(S) %*% matrix(1,nrow = 1, ncol = mmr_data[['R']]))
E <- variable(lower = 0, upper = mmr_data[['M']], dim = c(mmr_data[['N']],mmr_data[['S']]))

sigma = variable(lower = 0, upper = 1)
beta = normal(mean = 0, sd = sigma, dim = c(mmr_data[['K']],mmr_data[['R']]))

size = 50

mu1 = E %*% S * exp(mmr_data[['X']] %*% beta) / ((exp(mmr_data[['X']] %*% beta) %*% t(S[mmr_data[['S']],,drop=F])) %*% matrix(1,nrow = 1, ncol=104))
#mu1 = E[,-(mmr_data[['S']]),drop=F] %*% S[-(mmr_data[['S']]),,drop=F] + 
#  (E[,mmr_data[['S']],drop=F] %*% S[mmr_data[['S']],,drop = F]) * exp(mmr_data[['X']] %*% beta) / 
#((exp(mmr_data[['X']] %*% beta) %*% t(S[mmr_data[['S']],,drop=F])) %*% matrix(1,nrow = 1, ncol=104))

prob = size/(size + mu1)
distribution(mmr_data[['y']]) = negative_binomial(size = matrix(1, nrow = mmr_data[['N']], ncol = mmr_data[['R']]) * size, prob = prob)

#m <- model(S,E)
m <- model(S,E,beta,sigma)

# sampling
draws <- mcmc(m,n_samples = 1000, warmup = 500,verbose=F)
draws_2 <- mcmc(m, n_samples = 1000, warmup = 500, verbose = F)
draws_3 <- mcmc(m, n_samples = 1000, warmup = 500, verbose = F)
draws_4 <- mcmc(m, n_samples = 1000, warmup = 500, verbose = F)

draws[[2]] <- draws_2[[1]]
draws[[3]] <- draws_3[[1]]
draws[[4]] <- draws_4[[1]]

# Visualize the draws
library(bayesplot)
mcmc_trace(draws[,grep('S',colnames(draws[[1]]))[sample(1:(4*104),9)]])
mcmc_trace(draws[,grep('beta',colnames(draws[[1]]))[1:9]])
mcmc_trace(draws[,grep('sigma',colnames(draws[[1]])),drop=F])

draws_all <- draws[[2]]
S_est <- matrix(colMeans(draws_all[,grep('S',colnames(draws_all), fixed = T)]), nrow = mmr_data[['S']], ncol = mmr_data[['R']])
S_low <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = mmr_data[['S']], ncol = mmr_data[['R']])
S_high <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = mmr_data[['S']], ncol = mmr_data[['R']])
rownames(S_est) = rownames(S_low) = rownames(S_high) = paste0('Sig',1:mmr_data[['S']])
E_est <- matrix(colMeans(draws_all[,grep('E',colnames(draws_all))]), ncol = mmr_data[['S']])
beta_est <- matrix(colMeans(draws_all[,grep('beta',colnames(draws_all), fixed = T)]), nrow = mmr_data[['K']], ncol = mmr_data[['R']])
beta_low <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = mmr_data[['K']], ncol = mmr_data[['R']])
beta_high <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = mmr_data[['K']], ncol = mmr_data[['R']])
beta_var <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,var), nrow = mmr_data[['K']], ncol = mmr_data[['R']])

plot_subindel_wb(t(S_est), CI = T, low = t(S_low), high = t(S_high))
plot_subindel_wb(cbind(t(S_est), S_est[mmr_data[['S']],] * exp(beta_est[1,]),S_est[mmr_data[['S']],] * exp(beta_est[2,]),S_est[mmr_data[['S']],] * exp(beta_est[3,])))
interaction_effect_plot_human(beta_est[1,], lwd = 2, CI = T, low = beta_low[1,], high = beta_high[1,], at = c(-1,0,1))
interaction_effect_plot_human(beta_est[2,], lwd = 2, CI = T, low = beta_low[2,], high = beta_high[2,], at = c(-1,0,1))
interaction_effect_plot_human(beta_est[3,], lwd = 2, CI = T, low = beta_low[3,], high = beta_high[3,], at = c(-1,0,1), lwd = 5)

mu <- E_est[,-(mmr_data[['S']]),drop=F] %*% S_est[-(mmr_data[['S']]),,drop=F] + 
  (E_est[,mmr_data[['S']],drop=F] %*% S_est[mmr_data[['S']],,drop = F]) * exp(mmr_data[['X']] %*% beta_est)

plot(rowSums(mmr_data[['y']]),rowSums(mu), pch = 16)
abline(a=0,b=1,col='red',lty=2)

# Make a tsne plot
library(tsne)
mut_mat <- t(mmr_data[['y']])
cosdist <- function(x,y) {
  x0 <- x/sum(x)
  y0 <- y/sum(y)
  x0 %*% y0 / sqrt(x0%*%x0)/sqrt(y0%*%y0)
}
D <- as.dist(sapply(1:ncol(mut_mat), function(i) sapply(1:ncol(mut_mat), function(j) 1-cosdist(mut_mat[,i],mut_mat[,j]))))
set.seed(123)
t <- tsne(D)
rownames(t) <- colnames(mut_mat)

# Get the MMR and POLE status
coad <- mmr_data[['X']][,1]>0
stad <- mmr_data[['X']][,2]>0
ucec <- mmr_data[['X']][,3]>0

decomposition <- E_est
colnames(decomposition) <- c('Sig26','sig14-15','brca','sig20+ind=6')
for (i in 1:nrow(decomposition))
  decomposition[i,] <- decomposition[i,] / sum(decomposition[i,])

# Visualize it
library(RColorBrewer)
source('~/MMR/plotting functions/scatterpie.R')
library(mg14)
o1 <- order(colSums(mut_mat[,coad]),decreasing = T)
o2 <- order(colSums(mut_mat[,stad]),decreasing = T)
o3 <- order(colSums(mut_mat[,ucec]),decreasing = T)
#pdf('~/MMR_tsne.pdf', 15, 10)
#par(bty="n", mar=c(0,0,0,0))
plot(NA,NA, xlab="", ylab="", xlim=c(-25,25), ylim=c(-25,25), xaxt="n", yaxt="n", bty='n')
corr_scatterpie(t[ucec,1][o3], t[ucec,2][o3], p=decomposition[ucec,][o3,],
                r=sqrt(colSums(mut_mat)[ucec][o3])/75, labels=NA, col=brewer.pal(4,'Set3'),
                lty=0, circles=TRUE, lwd.circle=rep(2.5,sum(ucec)),
                lty.circle=rep(1,sum(ucec)), add=TRUE, col.circle = 'blue')
corr_scatterpie(t[coad,1][o1], t[coad,2][o1],
                p=decomposition[coad,][o1,],
                r=sqrt(colSums(mut_mat)[coad][o1])/75,
                labels=NA, col=brewer.pal(4,'Set3'), lty=0, circles=TRUE,
                lwd.circle=rep(0.5,sum(coad)),lty.circle=rep(1,sum(coad)),
                add=TRUE, col.circle = 'white')
corr_scatterpie(t[stad,1][o2], t[stad,2][o2], p=decomposition[stad,][o2,],
                r=sqrt(colSums(mut_mat)[stad][o2])/75, labels=NA, col=brewer.pal(4,'Set3'),
                lty=0, circles=TRUE, lwd.circle=rep(2.5,sum(stad)),
                lty.circle=rep(1,sum(stad)), add=TRUE, col.circle = 'black')
mg14:::.pie(x0=-20, y0=20, x=matrix(rep(1,4), nrow=1), r=sqrt(10000)/50, labels=colnames(decomposition),
            col=brewer.pal(4,'Set3'), lty=0, circles=TRUE, add=TRUE)
us <- par("usr")
pr <- (us[2]-us[1])/(us[4]-us[3])
fr <- par("pin")[1]/par("pin")[2]
for(i in c(1,10,100,1000,10000)){
  polygon(-20 + cos(seq(0,2*pi, l=100)) * sqrt(i)/75, 15+(1+sin(seq(0,2*pi, l=100))) * sqrt(i)/75 / pr * fr, col=NA)
  if (i>10) text(-20, 15 + 2*sqrt(i)/75 / pr * fr + 0.3,labels = as.character(i),cex=0.8)
}
polygon(-10 + cos(seq(0,2*pi, l=100)) * sqrt(i)/75, 23 + (1+sin(seq(0,2*pi, l=100))) * sqrt(i)/75 / pr * fr,
        lwd=0.5, col=NA, border = 'black')
polygon(-10 + cos(seq(0,2*pi, l=100)) * sqrt(i)/75, 19 + (1+sin(seq(0,2*pi, l=100))) * sqrt(i)/75 / pr * fr,
        lwd=1, col=NA, border = 'black')
polygon(-10 + cos(seq(0,2*pi, l=100)) * sqrt(i)/75, 15 +(1+sin(seq(0,2*pi, l=100))) * sqrt(i)/75 / pr * fr,
        lwd=1, col=NA, border = 'blue')
text(x = -20,y=13.8,labels = "Number of mutations")
text(x = -10,y=21.8,labels = "COAD")
text(x = -10,y=17.8,labels = "STAD")
text(x = -10,y=13.8,labels = "UCEC")
dev.off()

# PCA for the sake of sanity
p <- prcomp(x = mut_mat)
plot(p)
coords <- p$rotation[,1:2]
plot(NA,NA, xlab="", ylab="", xlim=c(0,0.2), ylim=c(-0.8,0.4), bty='n')
corr_scatterpie(coords[coad,1][o1], coords[coad,2][o1],
                p=decomposition[coad,][o1,],
                r=sqrt(colSums(mut_mat)[coad][o1])/10000,
                labels=NA, col=brewer.pal(4,'Set3'), lty=0, circles=TRUE,
                lwd.circle=rep(0.5,sum(coad)),lty.circle=rep(1,sum(coad)),
                add=TRUE, col.circle = 'white')
corr_scatterpie(coords[stad,1][o2], coords[stad,2][o2], p=decomposition[stad,][o2,],
                r=sqrt(colSums(mut_mat)[stad][o2])/10000, labels=NA, col=brewer.pal(4,'Set3'),
                lty=0, circles=TRUE, lwd.circle=rep(2.5,sum(stad)),
                lty.circle=rep(1,sum(stad)), add=TRUE, col.circle = 'black')
corr_scatterpie(coords[ucec,1][o3], coords[ucec,2][o3], p=decomposition[ucec,][o3,],
                r=sqrt(colSums(mut_mat)[ucec][o3])/10000, labels=NA, col=brewer.pal(4,'Set3'),
                lty=0, circles=TRUE, lwd.circle=rep(2.5,sum(ucec)),
                lty.circle=rep(1,sum(ucec)), add=TRUE, col.circle = 'blue')
mg14:::.pie(x0=-20, y0=20, x=matrix(rep(1,4), nrow=1), r=sqrt(10000)/50, labels=colnames(decomposition),
            col=brewer.pal(4,'Set3'), lty=0, circles=TRUE, add=TRUE)
us <- par("usr")
pr <- (us[2]-us[1])/(us[4]-us[3])
fr <- par("pin")[1]/par("pin")[2]
for(i in c(1,10,100,1000,10000)){
  polygon(-20 + cos(seq(0,2*pi, l=100)) * sqrt(i)/75, 15+(1+sin(seq(0,2*pi, l=100))) * sqrt(i)/75 / pr * fr, col=NA)
  if (i>10) text(-20, 15 + 2*sqrt(i)/75 / pr * fr + 0.3,labels = as.character(i),cex=0.8)
}
polygon(-10 + cos(seq(0,2*pi, l=100)) * sqrt(i)/75, 23 + (1+sin(seq(0,2*pi, l=100))) * sqrt(i)/75 / pr * fr,
        lwd=0.5, col=NA, border = 'black')
polygon(-10 + cos(seq(0,2*pi, l=100)) * sqrt(i)/75, 19 + (1+sin(seq(0,2*pi, l=100))) * sqrt(i)/75 / pr * fr,
        lwd=1, col=NA, border = 'black')
polygon(-10 + cos(seq(0,2*pi, l=100)) * sqrt(i)/75, 15 +(1+sin(seq(0,2*pi, l=100))) * sqrt(i)/75 / pr * fr,
        lwd=1, col=NA, border = 'blue')
text(x = -20,y=13.8,labels = "Number of mutations")
text(x = -10,y=21.8,labels = "COAD")
text(x = -10,y=17.8,labels = "STAD")
text(x = -10,y=13.8,labels = "UCEC")

#####################################
# Take all the MMR samples

mmr_info <- read.table('~/Downloads/ncomms15180-s16.xls.txt', sep ='\t', header = T)
mmr_info <- mmr_info[mmr_info$MSI.H_0.75>0 & mmr_info$MSS_0.75==0,]

mmr_data <- list()
mmr_data[['y']] <- do.call('rbind', lapply(unique(c(as.character(msi_status$sample[msi_status$MSI=='MSI-H']),
                                              as.character(mmr_info$Barcode))),
                                           function(x) bigmat[grep(x,rownames(bigmat)),]))
# BRCA, COAD, LIHC, OV, STAD, UCEC

mmr_data[['y']] <- mmr_data[['y']][rowSums(mmr_data[['y']][,1:96])>100,]
mmr_data[['y']] <- mmr_data[['y']][(sapply(1:nrow(mmr_data[['y']]), function(x) cosine(as.numeric(mmr_data[['y']][x,1:96]), cancer_signatures[,2])) < 0.75) &
                                     (sapply(1:nrow(mmr_data[['y']]), function(x) cosine(as.numeric(mmr_data[['y']][x,1:96]), cancer_signatures[,11])) < 0.75),]
mmr_data[['X']] <- model.matrix( ~ ., data = data.frame(v1 = as.character(metadata$project[match(rownames(mmr_data[['y']]),metadata$tumour)])))[,-1]
mmr_data[['X']] <- mmr_data[['X']][,colSums(mmr_data[['X']])>4]
mmr_data[['y']] <- mmr_data[['y']][rowSums(mmr_data[['X']])>0,]
mmr_data[['X']] <- mmr_data[['X']][rowSums(mmr_data[['X']])>0,] 
mmr_data[['N']] <- nrow(mmr_data[['y']])
mmr_data[['R']] <- ncol(mmr_data[['y']])
mmr_data[['K']] <- ncol(mmr_data[['X']])
mmr_data[['S']] <- 1
mmr_data[['M']] <- max(rowSums(mmr_data[['y']]))

# Greta model

S <- variable(lower = 0, upper = 1, dim = c(mmr_data[['S']],mmr_data[['R']]))
S <- S / (greta::rowSums(S) %*% matrix(1,nrow = 1, ncol = mmr_data[['R']]))
E <- variable(lower = 0, upper = mmr_data[['M']], dim = c(mmr_data[['N']],mmr_data[['S']]))

sigma = 0.1
beta = normal(mean = 0, sd = sigma, dim = c(mmr_data[['K']],mmr_data[['R']]))

size = 50

mu1 = E %*% S
#mu1 = E[,-(mmr_data[['S']]),drop=F] %*% S[-(mmr_data[['S']]),,drop=F] + 
#  (E[,mmr_data[['S']],drop=F] %*% S[mmr_data[['S']],,drop = F]) * exp(mmr_data[['X']] %*% beta) #/ 
#((exp(mmr_data[['X']] %*% beta) %*% t(S[mmr_data[['S']],,drop=F])) %*% matrix(1,nrow = 1, ncol=104))

prob = size/(size + mu1)
distribution(mmr_data[['y']]) = negative_binomial(size = matrix(1, nrow = mmr_data[['N']], ncol = mmr_data[['R']]) * size, prob = prob)

#m <- model(S,E)
m <- model(S,E,beta)

# sampling
draws <- mcmc(m,n_samples = 1000, warmup = 500)
draws_2 <- mcmc(m, n_samples = 1000, warmup = 500, verbose = F)
draws_3 <- mcmc(m, n_samples = 1000, warmup = 500, verbose = F)
draws_4 <- mcmc(m, n_samples = 1000, warmup = 500, verbose = F)

draws[[2]] <- draws_2[[1]]
draws[[3]] <- draws_3[[1]]
draws[[4]] <- draws_4[[1]]

# Visualize the draws
library(bayesplot)
mcmc_trace(draws[,grep('S',colnames(draws[[1]]))[sample(1:(mmr_data[['S']]*104),9)]])
mcmc_trace(draws[,grep('beta',colnames(draws[[1]]))[sample(1:(mmr_data[['K']]*104),9)]])


draws_all <- draws[[1]]
S_est <- matrix(colMeans(draws_all[,grep('S',colnames(draws_all), fixed = T)]), nrow = mmr_data[['S']], ncol = mmr_data[['R']])
S_low <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = mmr_data[['S']], ncol = mmr_data[['R']])
S_high <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = mmr_data[['S']], ncol = mmr_data[['R']])
rownames(S_est) = rownames(S_low) = rownames(S_high) = paste0('Sig',1:mmr_data[['S']])
E_est <- matrix(colMeans(draws_all[,grep('E',colnames(draws_all))]), ncol = mmr_data[['S']])
beta_est <- matrix(colMeans(draws_all[,grep('beta',colnames(draws_all), fixed = T)]), nrow = mmr_data[['K']], ncol = mmr_data[['R']])
beta_low <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = mmr_data[['K']], ncol = mmr_data[['R']])
beta_high <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = mmr_data[['K']], ncol = mmr_data[['R']])
beta_var <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,var), nrow = mmr_data[['K']], ncol = mmr_data[['R']])

#plot_subindel_wb(t(S_est), CI = T, low = t(S_low), high = t(S_high))
plot_subindel_wb(cbind(t(S_est[1:mmr_data[['S']],]),
                       S_est[mmr_data[['S']],] * exp(beta_est[1,]),
                       S_est[mmr_data[['S']],] * exp(beta_est[2,]),
                       S_est[mmr_data[['S']],] * exp(beta_est[3,]),
                       S_est[mmr_data[['S']],] * exp(beta_est[4,]),
                       S_est[mmr_data[['S']],] * exp(beta_est[5,])))
                       #S_est[mmr_data[['S']],] * exp(beta_est[6,]),
                       #S_est[mmr_data[['S']],] * exp(beta_est[7,]),S_est[mmr_data[['S']],] * exp(beta_est[8,]),S_est[mmr_data[['S']],] * exp(beta_est[9,]),
                       #S_est[mmr_data[['S']],] * exp(beta_est[10,]),S_est[mmr_data[['S']],] * exp(beta_est[11,])))
interaction_effect_plot_human(beta_est[1,], lwd = 2, CI = T, low = beta_low[1,], high = beta_high[1,], at = c(-1,0,1,2)) # BRCA
interaction_effect_plot_human(beta_est[2,], lwd = 2, CI = T, low = beta_low[2,], high = beta_high[2,], at = c(-1,0,1,2)) # ? CESC
interaction_effect_plot_human(beta_est[3,], lwd = 2, CI = T, low = beta_low[3,], high = beta_high[3,], at = c(-1,0,1,2)) # COAD
interaction_effect_plot_human(beta_est[4,], lwd = 2, CI = T, low = beta_low[4,], high = beta_high[4,], at = c(-1,0,1,2)) 
interaction_effect_plot_human(beta_est[5,], lwd = 2, CI = T, low = beta_low[5,], high = beta_high[5,], at = c(-1,0,1,2))
interaction_effect_plot_human(beta_est[6,], lwd = 2, CI = T, low = beta_low[6,], high = beta_high[6,], at = c(-1,0,1,2))
interaction_effect_plot_human(beta_est[7,], lwd = 2, CI = T, low = beta_low[7,], high = beta_high[7,], at = c(-1,0,1,2)) # LUSC
interaction_effect_plot_human(beta_est[8,], lwd = 2, CI = T, low = beta_low[8,], high = beta_high[8,], at = c(-1,0,1,2)) # OV
interaction_effect_plot_human(beta_est[9,], lwd = 2, CI = T, low = beta_low[9,], high = beta_high[9,], at = c(-1,0,1,2)) # READ (CG>TG)
interaction_effect_plot_human(beta_est[10,], lwd = 2, CI = T, low = beta_low[10,], high = beta_high[10,], at = c(-1,0,1,2)) # STAD
interaction_effect_plot_human(beta_est[11,], lwd = 2, CI = T, low = beta_low[11,], high = beta_high[11,], at = c(-1,0,1,2)) # UCEC

mu <- E_est[,-(mmr_data[['S']]),drop=F] %*% S_est[-(mmr_data[['S']]),,drop=F] + 
  (E_est[,mmr_data[['S']],drop=F] %*% S_est[mmr_data[['S']],,drop = F]) * exp(mmr_data[['X']] %*% beta_est)

plot(rowSums(mmr_data[['y']]),rowSums(mu), pch = 16)
abline(a=0,b=1,col='red',lty=2)

# Make a tsne plot
library(tsne)
mut_mat <- t(mmr_data[['y']])
cosdist <- function(x,y) {
  x0 <- x/sum(x)
  y0 <- y/sum(y)
  x0 %*% y0 / sqrt(x0%*%x0)/sqrt(y0%*%y0)
}
D <- as.dist(sapply(1:ncol(mut_mat), function(i) sapply(1:ncol(mut_mat), function(j) 1-cosdist(mut_mat[,i],mut_mat[,j]))))
set.seed(123)
t <- tsne(D, perplexity = 20)
rownames(t) <- colnames(mut_mat)

# Get the MMR and POLE status

decomposition <- E_est
colnames(decomposition) <- c('mmr1','mmr2','mmr3','brca')
for (i in 1:nrow(decomposition))
  decomposition[i,] <- decomposition[i,] / sum(decomposition[i,])

# Visualize it
library(RColorBrewer)
source('~/MMR/plotting functions/scatterpie.R')
library(mg14)
#pdf('~/MMR_tsne.pdf', 15, 10)
#par(bty="n", mar=c(0,0,0,0))
plot(NA,NA, xlab="", ylab="", xlim=c(-50,36), ylim=c(-37,42), xaxt="n", yaxt="n", bty='n')
for (j in 1:mmr_data[['K']]) {
  sel <- (mmr_data[['X']][,j]>0)
  o1 <- order(rowSums(mmr_data[['y']])[sel],decreasing = T)
  corr_scatterpie(t[sel,1][o1], t[sel,2][o1], p=decomposition[sel,][o1,],
                  r=sqrt(colSums(mut_mat)[sel][o1])/75, labels=NA, col=brewer.pal(4,'Set1'),
                  lty=0, circles=TRUE, lwd.circle=rep(2.5,sum(sel)),
                  lty.circle=rep(1,sum(sel)), add=TRUE, col.circle = brewer.pal(11,'Set3')[j])
}
mg14:::.pie(x0=-40, y0=-30, x=matrix(rep(1,4), nrow=1), r=sqrt(10000)/50, labels=colnames(decomposition),
            col=brewer.pal(4,'Set1'), lty=0, circles=TRUE, add=TRUE)
us <- par("usr")
pr <- (us[2]-us[1])/(us[4]-us[3])
fr <- par("pin")[1]/par("pin")[2]
#for(i in c(1,10,100,1000,10000)){
#  polygon(-20 + cos(seq(0,2*pi, l=100)) * sqrt(i)/75, 15+(1+sin(seq(0,2*pi, l=100))) * sqrt(i)/75 / pr * fr, col=NA)
#  if (i>10) text(-20, 15 + 2*sqrt(i)/75 / pr * fr + 0.3,labels = as.character(i),cex=0.8)
#}
for (j in 1:mmr_data[['K']]) {
  print(polygon(-45 + cos(seq(0,2*pi, l=100)) , 40-j*2 + (1+sin(seq(0,2*pi, l=100))) / pr * fr,
          lwd=0.5, col=brewer.pal(11,'Set3')[j]), border = 'black')
  print(text(x = -40,y=40-j*2+1,labels = colnames(mmr_data[['X']])[j]))
}
#dev.off()

####################################################################

# Take all the MMR samples and label them by MMR gene mutations (or methylation?)
# methylation = silencing = low expression
# hence missense mutations or low expression should be enough!

# Labels
mmr.genes <- c('MLH1','PMS2','MSH2','MSH3','MSH6')
#mmr.genes <- MMR.core[-6]

# Expression
metadata <- read.table('TCGA_caveman_patients_to_files.dat', sep='\t', header=T)
metadata1 <- metadata
check.exp <- function(z, tr, m = metadata1) {
  return(which(m[,match(z,colnames(m))] < tr))
}
for (i in 6:ncol(metadata1)) {
  exp_meds <- sapply(unique(metadata1$project), function(x) median(metadata1[metadata1$project==x,i], na.rm = T))
  names(exp_meds) <- unique(metadata1$project)
  metadata1[,i] <- metadata1[,i] / exp_meds[metadata1$project]
  #  metadata[,i] <- rank(metadata[,i], na.last='keep') / length(metadata[,i])
}
exp.mmr <- lapply(mmr.genes, function(x) metadata1$tumour[sort(check.exp(x,0.2))])

#samples <- list()
#for (x in c(mmr.genes)) {
#  tmp <- read.delim(paste0('~/TCGAmutations/results/',x,'.txt'), sep='\t', header=T)
#  tmp <- tmp[tmp$IMPACT %in% c('HIGH') & grepl(x, tmp$SYMBOL),]
#  if (nrow(tmp)==0) next
#  tmp2 <- sapply(as.character(tmp$Uploaded_variation), function(x) unlist(strsplit(x,split='[:]'))[1])
#  tmp2 <- paste0('TCGA',substr(as.character(tmp2),5,nchar(tmp2)-2))
#  samples[[x]] <- unique(tmp2)
#}

# Germline

SNPs <- list()
tmp <- read.table('~/TCGAmutations/results/germline_MMR_38.txt', sep = '\t', header = T)
for (x in mmr.genes) {
  SNPs[[x]] <- unique(as.character(tmp$Uploaded_variation[tmp$IMPACT %in% c('HIGH') & grepl(x, tmp$SYMBOL)]))
}
SNPs <- SNPs[sapply(SNPs,length)>0]

#ddb2_hom <- NULL
intsam <- list(); z = 0
for (cancername in unlist(strsplit(as.character(unique(metadata$project)),split = ' '))[-c(1,5)]) {
  coad <- as.character(metadata$tumour[metadata$project==cancername])
#snps_per_sample <- matrix(0,nrow = length(coad), ncol = length(SNPs), dimnames = list(as.character(coad), names(SNPs)))
  for (genename in names(SNPs)) {
    f <- list.files(paste0('~/SNP/',genename,'/',cancername))
    for (i in 1:length(coad)) {
      file <- f[grep(coad[i],f)]
      if (length(file)==0) next
      if (file.info(paste0('~/SNP/',genename,'/',cancername,'/',file))$size==0) next
      vcf <- read.delim(paste0('~/SNP/',genename,'/',cancername,'/',file), header = F)
      identifier <- paste0(vcf$V1,':',vcf$V2,'_',vcf$V4,'/',vcf$V5) 
      het <- substr(vcf$V10,1,3)[identifier %in% SNPs[[genename]]]
      if (length(het)>0) {
        z = z + 1
        #if (genename=='DDB2' & het=='1|1') ddb2_hom <- c(ddb2_hom,as.character(lung[i]))
        print(c(as.character(coad[i]),het,genename))
        intsam[[z]] <- c(as.character(coad[i]),het,genename)
      }
      #snps_per_sample[as.character(lung[i]), genename] <- length(intersect(identifier, SNPs[[genename]]))
    }
    print(genename)
  }
}
sapply(intsam, function(x) x[1])

# Damaging mutations

#copynumber <- read.delim('filtered.combined.segments.txt', sep = '\t', header = T)
#copynumber <- GRanges(copynumber)
#ids <- read.table('~/Downloads/ICGCtoTCGA.tsv', sep = '\t', header = T)
#samples <- list()
#for (x in mmr.genes) {
#  tmp <- read.delim(paste0('~/TCGAmutations/results/',x,'.txt'), sep='\t', header=T)
#  tmp <- tmp[tmp$IMPACT %in% c('HIGH') & grepl(x, tmp$SYMBOL),]
#  if (nrow(tmp)==0) next
#  mutations <- read.delim(paste0('/Volumes/nobackup/MMR/',x,'.dat'), sep='\t', header = T)
#  identifier <- paste0(mutations$Sample,'-',mutations$Chrom,':',mutations$Pos,'_',mutations$Ref,'/',mutations$Alt)
#  a <- sapply(as.character(tmp$Uploaded_variation), function(y) substr(unlist(strsplit(y,split='[:]'))[1],1,12))
#  vaf <- mutations$PM.Tum[match(as.character(tmp$Uploaded_variation), identifier)]
#  samples[[x]] <- paste0('TCGA', substr(a,5,12))[vaf>0.4]
#  samples[[x]] <- unlist(sapply(samples[[x]], function(y) metadata$tumour[grep(y,metadata$patient_id)]))
#  cn <- copynumber[queryHits(findOverlaps(copynumber, genes_of_interest[genes_of_interest$name==x]))]
#  cn <- cn[cn$tissue %in% names(samples[[x]])]
#  samples[[paste0(x,'LOH')]] <- samples[[x]][as.character(cn$tissue[cn$minor==0])]
#}

#mmr_info <- read.table('~/Downloads/ncomms15180-s16.xls.txt', sep ='\t', header = T)
#mmr_info <- mmr_info[mmr_info$MSI.H_0.75>0 & mmr_info$MSS_0.75==0,]

#mmr_data <- list()

#mmr_data[['y']] <- do.call('rbind', lapply(unique(c(as.character(msi_status$sample[msi_status$MSI=='MSI-H']),
#                                                    as.character(mmr_info$Barcode))),
#                                           function(x) bigmat[grep(x,rownames(bigmat)),]))

#new.mmr <- unique(c(as.character(unlist(exp.mmr)),
#                    as.character(unlist(samples[paste0(MMR.core,'_hom')])),
                    #as.character(unlist(samples[paste0(MMR.core[-2],'_het')]))))
#                    unlist(sapply(unique(c(as.character(msi_status$sample[msi_status$MSI=='MSI-H']),
#                                    as.character(mmr_info$Barcode))),
#                           function(x) rownames(bigmat)[grep(x,rownames(bigmat))]))))
#new.mmr <- unique(substr(new.mmr,1,12))

#################################################################################################################

bigmat <- read.table('TCGA.caveman.matrix.dat',sep='\t')
bigmat <- bigmat[rowSums(bigmat) > 1000 & rowSums(bigmat) < 20000,]
alternative_rownames_bigmat <- paste0('TCGA',substr(rownames(bigmat),5,nchar(rownames(bigmat))))

mmr_data <- list()
new.mmr <- unique(c(as.character(unlist(samples[paste0(mmr.genes,'_hom')])),
             as.character(unlist(samples[paste0(mmr.genes,'_het')])),
             as.character(sapply(intsam, function(x) x[1])),
             substr(unlist(exp.mmr),1,12),
             msi_status$sample[msi_status$MSI == 'MSI-H']))
             #substr(rownames(mmr_data[['y']]),1,12)))
new.mmr <- alternative_rownames_bigmat[substr(alternative_rownames_bigmat,1,12) %in% new.mmr]

new.mmr <- new.mmr[which(sapply(new.mmr, function(x) cosine(as.numeric(bigmat[match(x,alternative_rownames_bigmat),1:96]), cancer_signatures[,2])) < 0.75 &
                           sapply(new.mmr, function(x) cosine(as.numeric(bigmat[match(x,alternative_rownames_bigmat),1:96]), cancer_signatures[,13])) < 0.75 &
                           sapply(new.mmr, function(x) cosine(as.numeric(bigmat[match(x,alternative_rownames_bigmat),1:96]), cancer_signatures[,7])) < 0.75 &
                           sapply(new.mmr, function(x) cosine(as.numeric(bigmat[match(x,alternative_rownames_bigmat),1:96]), cancer_signatures[,10])) < 0.75 &
                           sapply(new.mmr, function(x) cosine(as.numeric(bigmat[match(x,alternative_rownames_bigmat),1:96]), cancer_signatures[,4])) < 0.75)] # 370
table(metadata$project[match(new.mmr, metadata$tumour)])


mmr_data <- list()
mmr_data[['y']] <- bigmat[match(new.mmr,alternative_rownames_bigmat),]
mmr_data[['X']] <- model.matrix( ~ .,
    data = data.frame(v1 = as.character(metadata$project[match(paste0('TCGA',substr(rownames(mmr_data[['y']]),5,nchar(rownames(mmr_data[['y']])))),
                        metadata$tumour)])))[,-1]
colSums(mmr_data[['X']])
mmr_data[['X']] <- cbind(mmr_data[['X']], v1ACC = as.numeric(rowSums(mmr_data[['X']])==0))
# COAD - DEFAULT
mmr_data[['X']] <- cbind(mmr_data[['X']], ADRENO = as.numeric(mmr_data[['X']][,'v1ACC']>0),
                                          #BRAIN = as.numeric((mmr_data[['X']][,'v1GBM'] + mmr_data[['X']][,'v1LGG']) > 0),
                                          BREAST = as.numeric(mmr_data[['X']][,'v1BRCA']>0),
                                          #BLADDER = as.numeric(mmr_data[['X']][,'v1BLCA'] > 0),
                                          CERVIX = as.numeric(mmr_data[['X']][,'v1CESC']>0),
                                          HNSC = as.numeric(mmr_data[['X']][,'v1HNSC']>0),
                                          STOMACH = as.numeric(mmr_data[['X']][,'v1STAD']>0),
                                          COLORECT = as.numeric((mmr_data[['X']][,'v1COAD'] + mmr_data[['X']][,'v1READ']) > 0),
                                          #KIRC = as.numeric((mmr_data[['X']][,'v1KIRC']) > 0),
                                          LIVER = as.numeric(mmr_data[['X']][,'v1LIHC'] > 0),
                                          LUNG = as.numeric((mmr_data[['X']][,'v1LUSC']+mmr_data[['X']][,'v1LUAD']) > 0),
                                          #OVARY = as.numeric(mmr_data[['X']][,'v1OV']>0),
                                          PROSTATE = as.numeric(mmr_data[['X']][,'v1PRAD']>0),
                                          UTERUS = as.numeric((mmr_data[['X']][,'v1UCEC'] + mmr_data[['X']][,'v1UCS'])>0))
mmr_data[['X']] <- mmr_data[['X']][,-c(1:22)]
#mmr_data[['X']] <- mmr_data[['X']][,colSums(mmr_data[['X']]) > 10]
ind <- which(rowSums(mmr_data[['X']])>0)
mmr_data[['X']] <- mmr_data[['X']][ind,]
mmr_data[['y']] <- mmr_data[['y']][ind,]

mmr_data[['X']] <- mmr_data[['X']][,-match('COLORECT',colnames(mmr_data[['X']]))]
#mmr_data[['X']] <- cbind(mmr_data[['X']],AGE = ids$donor_age_at_diagnosis[
#  match(metadata$patient_id[match(rownames(mmr_data[['y']]),
#                                  metadata$tumour)], 
#        ids$submitted_donor_id)])
#mmr_data[['y']] <- mmr_data[['y']][!is.na(mmr_data[['X']][,'AGE']),] 
#mmr_data[['X']] <- mmr_data[['X']][!is.na(mmr_data[['X']][,'AGE']),] 
#mmr_data[['X']] <- mmr_data[['X']][,colSums(mmr_data[['X']])>10]
#mmr_data[['y']] <- mmr_data[['y']][rowSums(mmr_data[['X']])>0,]
#mmr_data[['X']] <- mmr_data[['X']][rowSums(mmr_data[['X']])>0,]
#mmr_data[['X']][,'AGE'] <- scale(mmr_data[['X']][,'AGE'], center = T)

mmr_data[['N']] <- nrow(mmr_data[['y']]) 
mmr_data[['R']] <- ncol(mmr_data[['y']]) 
mmr_data[['K']] <- ncol(mmr_data[['X']])
mmr_data[['S']] <- 3
mmr_data[['M']] <- max(rowSums(mmr_data[['y']]))

# Greta model

save(mmr_data, file = 'MMR_data.RData')

S <- variable(lower = 0, upper = 1, dim = c(mmr_data[['S']],mmr_data[['R']]))
S <- S / (greta::rowSums(S) %*% matrix(1,nrow = 1, ncol = mmr_data[['R']]))
E <- variable(lower = 0, upper = mmr_data[['M']], dim = c(mmr_data[['N']],mmr_data[['S']]))

sigma = 0.5
beta = normal(mean = 0, sd = sigma, dim = c(mmr_data[['K']],mmr_data[['R']]))

size = 50

#mu1 = E %*% S
mu1 = E[,-(mmr_data[['S']]),drop=F] %*% S[-(mmr_data[['S']]),,drop=F] + 
  (E[,mmr_data[['S']],drop=F] %*% S[mmr_data[['S']],,drop = F]) * exp(mmr_data[['X']] %*% beta) / 
((exp(mmr_data[['X']] %*% beta) %*% t(S[mmr_data[['S']],,drop=F])) %*% matrix(1,nrow = 1, ncol=104))

mu1 = E %*% S * exp(mmr_data[['X']] %*% beta) / ((exp(mmr_data[['X']] %*% beta) %*% t(S[mmr_data[['S']],,drop=F])) %*% matrix(1,nrow = 1, ncol=104))

prob = size/(size + mu1)
distribution(mmr_data[['y']]) = negative_binomial(size = matrix(1, nrow = mmr_data[['N']], ncol = mmr_data[['R']]) * size, prob = prob)

#m <- model(S,E)
m <- model(S,E,beta)

# sampling
draws <- mcmc(m,n_samples = 500, warmup = 500, chains = 4)
#load('/Volumes/yoda1/MMR_sampling_4_sigs_no_normalization.RData')

# Visualize the draws
library(bayesplot)
mcmc_trace(draws[,grep('S',colnames(draws[[1]]))[sample(1:(mmr_data[['S']]*104),9)]])
mcmc_trace(draws[,grep('beta',colnames(draws[[1]]))[sample(1:(mmr_data[['K']]*104),9)]])

save(mmr_data, draws, file = '~/MMR_draws_norm_3_sigs.RData')

draws_all <- draws[[1]] # 1,3
S_est <- matrix(colMeans(draws_all[,grep('S',colnames(draws_all), fixed = T)]), nrow = mmr_data[['S']], ncol = mmr_data[['R']])
S_var <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,var), nrow = mmr_data[['S']], ncol = mmr_data[['R']])
S_low <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = mmr_data[['S']], ncol = mmr_data[['R']])
S_high <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = mmr_data[['S']], ncol = mmr_data[['R']])
rownames(S_est) = rownames(S_low) = rownames(S_high) = paste0('Sig',1:mmr_data[['S']])
E_est <- matrix(colMeans(draws_all[,grep('E',colnames(draws_all))]), ncol = mmr_data[['S']])
beta_est <- matrix(colMeans(draws_all[,grep('beta',colnames(draws_all), fixed = T)]), nrow = mmr_data[['K']], ncol = mmr_data[['R']])
beta_low <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = mmr_data[['K']], ncol = mmr_data[['R']])
beta_high <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = mmr_data[['K']], ncol = mmr_data[['R']])
beta_var <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,var), nrow = mmr_data[['K']], ncol = mmr_data[['R']])

#plot_subindel_wb(t(S_est), CI = T, low = t(S_low), high = t(S_high))
rownames(beta_est) <- colnames(mmr_data[['X']])
plot_subindel_wb(t(S_est), CI = T, low = t(S_low), high = t(S_high))
plot_all_changes(plot_subindel_wb, S_est, beta_est)

pdf('~/MMR_effects.pdf',10,4)
interaction_effect_plot_human(c(beta_est[1,],0), lwd = 2, CI = T, low = beta_low[1,], high = beta_high[1,], at = c(-1,0,1), plot_main = colnames(mmr_data[['X']])[1]) 
interaction_effect_plot_human(c(beta_est[2,],0), lwd = 2, CI = T, low = beta_low[2,], high = beta_high[2,], at = c(-1,0,1), plot_main = colnames(mmr_data[['X']])[2]) 
interaction_effect_plot_human(c(beta_est[3,],0), lwd = 2, CI = T, low = beta_low[3,], high = beta_high[3,], at = c(-1,0,1), plot_main = colnames(mmr_data[['X']])[3]) 
interaction_effect_plot_human(c(beta_est[4,],0), lwd = 2, CI = T, low = beta_low[4,], high = beta_high[4,], at = c(-1,0,1), plot_main = colnames(mmr_data[['X']])[4]) 
interaction_effect_plot_human(c(beta_est[5,],0), lwd = 2, CI = T, low = beta_low[5,], high = beta_high[5,], at = c(-1,0,1), plot_main = colnames(mmr_data[['X']])[5])
#interaction_effect_plot_human(c(beta_est[6,],0), lwd = 2, CI = T, low = beta_low[6,], high = beta_high[6,], at = c(-1,0,1), plot_main = colnames(mmr_data[['X']])[6])
#interaction_effect_plot_human(c(beta_est[7,],0), lwd = 2, CI = T, low = beta_low[7,], high = beta_high[7,], at = c(-1,0,1), plot_main = colnames(mmr_data[['X']])[7]) 
#interaction_effect_plot_human(c(beta_est[8,],0), lwd = 2, CI = T, low = beta_low[8,], high = beta_high[8,], at = c(-1,0,1), plot_main = colnames(mmr_data[['X']])[8]) 

dev.off()

mu <- E_est[,-(mmr_data[['S']]),drop=F] %*% S_est[-(mmr_data[['S']]),,drop=F] + 
  (E_est[,mmr_data[['S']],drop=F] %*% S_est[mmr_data[['S']],,drop = F]) * exp(mmr_data[['X']] %*% beta_est)

mu <- E_est %*% S_est * exp(mmr_data[['X']] %*% beta_est) / ((exp(mmr_data[['X']] %*% beta_est) %*% t(S_est[mmr_data[['S']],,drop=F])) %*% matrix(1,nrow = 1, ncol=104))

pdf('~/MMR_observed_predicted.pdf',5,5)
plot(rowSums(mmr_data[['y']]),rowSums(mu), pch = 16,
     xlab = 'Observed', ylab= 'Predicted')
abline(a=0,b=1,col='red',lty=2)
dev.off()

# Make a tsne plot
library(tsne)
mut_mat <- t(mmr_data[['y']])
cosdist <- function(x,y) {
  x0 <- x/sum(x)
  y0 <- y/sum(y)
  x0 %*% y0 / sqrt(x0%*%x0)/sqrt(y0%*%y0)
}
D <- as.dist(sapply(1:ncol(mut_mat), function(i) sapply(1:ncol(mut_mat), function(j) 1-cosdist(mut_mat[,i],mut_mat[,j]))))
set.seed(123)
t <- tsne(D, perplexity = 20)
rownames(t) <- colnames(mut_mat)

# Get the MMR and POLE status

decomposition <- E_est
colnames(decomposition) <- c('sig6','other')
for (i in 1:nrow(decomposition))
  decomposition[i,] <- decomposition[i,] / sum(decomposition[i,])

# Visualize it
library(RColorBrewer)
source('~/MMR/plotting functions/scatterpie.R')
library(mg14)
pdf('~/MMR_tsne_bethesda_and_expression_and_mutations_4_sigs.pdf', 10, 10)
#par(bty="n", mar=c(0,0,0,0))
cancers <- colnames(mmr_data[['X']])
plot(NA,NA, xlab="", ylab="", xlim=c(-50,35), ylim=c(-33,34), xaxt="n", yaxt="n", bty='n')
for (j in 1:length(cancers)) {
  sel <- (mmr_data[['X']][,cancers[j]]>0)
  o1 <- order(rowSums(mmr_data[['y']])[sel],decreasing = T)
  corr_scatterpie(t[sel,1][o1], t[sel,2][o1], p=decomposition[sel,][o1,],
                  r=sqrt(colSums(mut_mat)[sel][o1])/30, labels=NA, col=brewer.pal(4,'Set1'),
                  lty=0, circles=TRUE, lwd.circle=rep(2.5,sum(sel)),
                  lty.circle=rep(1,sum(sel)), add=TRUE, col.circle = brewer.pal(11,'Set3')[j])
}
sel <- (rowSums(mmr_data[['X']])==0)
o1 <- order(rowSums(mmr_data[['y']])[sel],decreasing = T)
corr_scatterpie(t[sel,1][o1], t[sel,2][o1], p=decomposition[sel,][o1,],
                r=sqrt(colSums(mut_mat)[sel][o1])/30, labels=NA, col=brewer.pal(4,'Set1'),
                lty=0, circles=TRUE, lwd.circle=rep(2.5,sum(sel)),
                lty.circle=rep(1,sum(sel)), add=TRUE, col.circle = 'black')#brewer.pal(11,'Set3')[j])
mg14:::.pie(x0=-40, y0=32, x=matrix(rep(1,2), nrow=1), r=sqrt(10000)/20, labels=colnames(decomposition),
            col=brewer.pal(2,'Set1'), lty=0, circles=TRUE, add=TRUE, cex = 0.8)
us <- par("usr")
pr <- (us[2]-us[1])/(us[4]-us[3])
fr <- par("pin")[1]/par("pin")[2]
#for(i in c(1,10,100,1000,10000)){
#  polygon(-20 + cos(seq(0,2*pi, l=100)) * sqrt(i)/75, 15+(1+sin(seq(0,2*pi, l=100))) * sqrt(i)/75 / pr * fr, col=NA)
#  if (i>10) text(-20, 15 + 2*sqrt(i)/75 / pr * fr + 0.3,labels = as.character(i),cex=0.8)
#}
for (j in 1:length(cancers)) {
  print(polygon(-30 + cos(seq(0,2*pi, l=100)) , -10-j*3 + (1+sin(seq(0,2*pi, l=100))) / pr * fr,
                lwd=0.5, col=brewer.pal(11,'Set3')[j]), border = 'black')
  print(text(x = -20,y=-10-j*3+1,labels = cancers[j]))
}
dev.off()

pdf('MMR_signature.pdf',10,4)
tmp <- data.frame(t(S_est[1,,drop = F]))
plot_subindel_wb(tmp, CI = T, low = t(S_low[1,,drop=F]), high = t(S_high[1,,drop=F])) + theme(panel.grid = element_blank())
dev.off()
pdf('MMR_signatures.pdf',12,8)
tmp <- data.frame(t(S_est))
tmp_low <- data.frame(t(S_low))
tmp_high <- data.frame(t(S_high))
colnames(tmp) = colnames(tmp_high) = colnames(tmp_low) <- c('Sig26', 'BRCA','POLE','MMR')
plot_subindel_wb(t(S_est), CI = T, low = tmp_low, high = tmp_high) + theme(panel.grid = element_blank())
dev.off()
pdf('MMR_alterations.pdf',12,20)
tmp <- data.frame(cbind(S_est[mmr_data[['S']],],
                      S_est[mmr_data[['S']],]* exp(beta_est[1,]),
                       S_est[mmr_data[['S']],] * exp(beta_est[2,]),
                       S_est[mmr_data[['S']],] * exp(beta_est[3,]),
                       #S_est[mmr_data[['S']],] * exp(beta_est[4,]),
                       S_est[mmr_data[['S']],] * exp(beta_est[5,])))
tmp_low <- data.frame(t(S_low))
tmp_high <- data.frame(t(S_high))
colnames(tmp) = colnames(tmp_high) = colnames(tmp_low) <- c('Sig26', 'BRCA','POLE','MMR')
plot_subindel_wb(t(S_est), CI = T, low = tmp_low, high = tmp_high) + theme(panel.grid = element_blank())
dev.off()

rownames(beta_est) <- colnames(mmr_data[['X']])
altered_signatures <- data.frame(S_est[mmr_data[['S']],],
                                 S_est[mmr_data[['S']],] * exp(beta_est[1,]),
                                 S_est[mmr_data[['S']],] * exp(beta_est[2,]),
                                 S_est[mmr_data[['S']],] * exp(beta_est[3,]),
                                 S_est[mmr_data[['S']],] * exp(beta_est[4,]),
                                 S_est[mmr_data[['S']],] * exp(beta_est[5,]),
                                 S_est[mmr_data[['S']],] * exp(beta_est[6,]),
                                 S_est[mmr_data[['S']],] * exp(beta_est[7,]),
                                 S_est[mmr_data[['S']],] * exp(beta_est[8,]))
colnames(altered_signatures) <- c('COLORECT',rownames(beta_est))
var_tmp <- t(rbind(S_var[mmr_data[['S']],],exp(2*beta_est) * S_var[mmr_data[['S']],] + S_est[mmr_data[['S']],]**2 * exp(2*beta_est) * beta_var)) # variance of S_est[1,] * exp(beta_est[1,])

altered_low <- altered_signatures - 1.96*sqrt(var_tmp)
altered_high <- altered_signatures + 1.96*sqrt(var_tmp)

plot_subindel_wb(altered_signatures, CI = T, low = altered_low, high = altered_high)

pdf('MMR_signature.pdf',12,7)
plot_subindel_wb(altered_signatures[,c(1,5,7,9)], CI = T, low = altered_low[,c(1,5,7,9)], high = altered_high[,c(1,5,7,9)]) + theme(panel.grid = element_blank())
dev.off()

altered_sums <- colSums(altered_signatures)
altered_sums_var <- colSums(var_tmp)

new_var_tmp <- sqrt(1/altered_sums**2 * var_tmp + altered_signatures**2 / altered_sums**4 * altered_sums_var)
new_altered_signatures <- apply(altered_signatures,2,function(x) x/sum(x))

plot_subindel_wb(new_altered_signatures, CI = T, low = new_altered_signatures - 1.96*new_var_tmp, high = new_altered_signatures + 1.96 * new_var_tmp)



bigmat <- read.table('TCGA.caveman.matrix.dat',sep='\t')
bigmat <- bigmat[rowSums(bigmat) < 20000 & rowSums(bigmat) > 50,]
alternative_rownames_bigmat <- paste0('TCGA',substr(rownames(bigmat),5,nchar(rownames(bigmat))))

mean_c_to_t <- list()
mean_c_to_t_sd <- list()
mean_c_to_t_low <- list()
mean_c_to_t_high <- list()
mmr_mean_c_to_t <- list()
mmr_mean_c_to_t_sd <- list()
mmr_mean_c_to_t_low <- list()
mmr_mean_c_to_t_high <- list()
mean_del <- list()
mean_del_sd <- list()
mean_del_low <- list()
mean_del_high <- list()
mmr_mean_del <- list()
mmr_mean_del_sd <- list()
mmr_mean_del_low <- list()
mmr_mean_del_high <- list()
mean_ins <- list()
mean_ins_sd <- list()
mean_ins_low <- list()
mean_ins_high <- list()
mmr_mean_ins <- list()
mmr_mean_ins_sd <- list()
mmr_mean_ins_low <- list()
mmr_mean_ins_high <- list()
age <- ids$donor_age_at_diagnosis[match(substr(alternative_rownames_bigmat,1,12),ids$submitted_donor_id)]
ctype <- sapply(alternative_rownames_bigmat, function(x) metadata$project[match(x,metadata$tumour)])
for (proj in c('COAD','BRCA','CESC','HNSC','STAD','LIHC','LUSC','PRAD','UCEC')) {
  
  rest_of_samples <- setdiff(alternative_rownames_bigmat[ctype == proj],new.mmr)
  
  #mean_prof[[proj]] <- rowMeans(apply(bigmat[match(rest_of_samples,alternative_rownames_bigmat),],1,function(x) x/sum(x)))
  
  #mean_sd[[proj]] <- apply(apply(bigmat[match(rest_of_samples,alternative_rownames_bigmat),],1,function(x) x/sum(x)),1,sd) / length(rest_of_samples)
  
  mean_c_to_t[[proj]] <- mean(rowSums(bigmat[match(rest_of_samples,alternative_rownames_bigmat),c(35,39,43,47)]) / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T)
  mean_c_to_t_sd[[proj]] <- sd(rowSums(bigmat[match(rest_of_samples,alternative_rownames_bigmat),c(35,39,43,47)]) / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T) #/ 
    #sum(!is.na(rowSums(bigmat[match(rest_of_samples,alternative_rownames_bigmat),c(35,39,43,47)]) / age[match(rest_of_samples,alternative_rownames_bigmat)]))
  #mean_c_to_t_low[[proj]] <- quantile(rowSums(bigmat[match(rest_of_samples,alternative_rownames_bigmat),c(35,39,43,47)]) / age[match(rest_of_samples,alternative_rownames_bigmat)],
  #                                    0.05, na.rm = T)
  #mean_c_to_t_high[[proj]] <- quantile(rowSums(bigmat[match(rest_of_samples,alternative_rownames_bigmat),c(35,39,43,47)]) / age[match(rest_of_samples,alternative_rownames_bigmat)],
  #                                 0.95, na.rm = T)
  
  mean_del[[proj]] <- mean(bigmat[match(rest_of_samples,alternative_rownames_bigmat),97] / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T)
  mean_del_sd[[proj]] <- sd(bigmat[match(rest_of_samples,alternative_rownames_bigmat),97] / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T) #/ 
    #sum(!is.na(bigmat[match(rest_of_samples,alternative_rownames_bigmat),97] / age[match(rest_of_samples,alternative_rownames_bigmat)]))
#  mean_del_low[[proj]] <- quantile(bigmat[match(rest_of_samples,alternative_rownames_bigmat),97] / age[match(rest_of_samples,alternative_rownames_bigmat)],
 #                                     0.05, na.rm = T)
#  mean_del_high[[proj]] <- quantile(bigmat[match(rest_of_samples,alternative_rownames_bigmat),97] / age[match(rest_of_samples,alternative_rownames_bigmat)],
  #                                 0.95, na.rm = T)
  
  
  mean_ins[[proj]] <- mean((bigmat[match(rest_of_samples,alternative_rownames_bigmat),102]) / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T)
  mean_ins_sd[[proj]] <- sd((bigmat[match(rest_of_samples,alternative_rownames_bigmat),102]) / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T) #/ 
#    sum(!is.na((bigmat[match(rest_of_samples,alternative_rownames_bigmat),102]) / age[match(rest_of_samples,alternative_rownames_bigmat)]))
##  mean_ins_low[[proj]] <- quantile(bigmat[match(rest_of_samples,alternative_rownames_bigmat),102] / age[match(rest_of_samples,alternative_rownames_bigmat)],
#                                      0.05, na.rm = T)
#  mean_ins_high[[proj]] <- quantile(bigmat[match(rest_of_samples,alternative_rownames_bigmat),102] / age[match(rest_of_samples,alternative_rownames_bigmat)],
#                                   0.95, na.rm = T)
  
  rest_of_samples <- intersect(alternative_rownames_bigmat[ctype == proj],new.mmr)
  mmr_mean_c_to_t[[proj]] <- mean(rowSums(bigmat[match(rest_of_samples,alternative_rownames_bigmat),c(35,39,43,47)]) / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T)
  mmr_mean_c_to_t_sd[[proj]] <- sd(rowSums(bigmat[match(rest_of_samples,alternative_rownames_bigmat),c(35,39,43,47)]) / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T) #/ 
  #  sum(!is.na(rowSums(bigmat[match(rest_of_samples,alternative_rownames_bigmat),c(35,39,43,47)]) / age[match(rest_of_samples,alternative_rownames_bigmat)]))
  #mmr_mean_c_to_t_low[[proj]] <- quantile(rowSums(bigmat[match(rest_of_samples,alternative_rownames_bigmat),c(35,39,43,47)]) / age[match(rest_of_samples,alternative_rownames_bigmat)],
  #                                    0.05, na.rm = T)
  #mmr_mean_c_to_t_high[[proj]] <- quantile(rowSums(bigmat[match(rest_of_samples,alternative_rownames_bigmat),c(35,39,43,47)]) / age[match(rest_of_samples,alternative_rownames_bigmat)],
  #                                 0.95, na.rm = T)
  
  mmr_mean_del[[proj]] <- mean(bigmat[match(rest_of_samples,alternative_rownames_bigmat),97] / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T)
  mmr_mean_del_sd[[proj]] <- sd(bigmat[match(rest_of_samples,alternative_rownames_bigmat),97] / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T) #/ 
  #  sum(!is.na(bigmat[match(rest_of_samples,alternative_rownames_bigmat),97] / age[match(rest_of_samples,alternative_rownames_bigmat)]))
  #mmr_mean_del_low[[proj]] <- quantile(bigmat[match(rest_of_samples,alternative_rownames_bigmat),97] / age[match(rest_of_samples,alternative_rownames_bigmat)],
  #                                    0.05, na.rm = T)
  #mmr_mean_del_high[[proj]] <- quantile(bigmat[match(rest_of_samples,alternative_rownames_bigmat),97] / age[match(rest_of_samples,alternative_rownames_bigmat)],
  #                                 0.95, na.rm = T)
  
  mmr_mean_ins[[proj]] <- mean((bigmat[match(rest_of_samples,alternative_rownames_bigmat),102]) / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T)
  mmr_mean_ins_sd[[proj]] <- sd((bigmat[match(rest_of_samples,alternative_rownames_bigmat),102]) / age[match(rest_of_samples,alternative_rownames_bigmat)], na.rm = T) #/ 
  #  sum(!is.na((bigmat[match(rest_of_samples,alternative_rownames_bigmat),102]) / age[match(rest_of_samples,alternative_rownames_bigmat)]))
  #mmr_mean_ins_low[[proj]] <- quantile(bigmat[match(rest_of_samples,alternative_rownames_bigmat),102] / age[match(rest_of_samples,alternative_rownames_bigmat)],
  #                                    0.05, na.rm = T)
  #mmr_mean_ins_high[[proj]] <- quantile(bigmat[match(rest_of_samples,alternative_rownames_bigmat),102] / age[match(rest_of_samples,alternative_rownames_bigmat)],
  #                                 0.95, na.rm = T)
}

pdf('MMR_C_T_scatterplot.pdf',4,4)
plot(y = unlist(mmr_mean_c_to_t), 
     x = unlist(mean_c_to_t),pch = 16,col='red',
     xlab = 'Non-MMR', ylab = 'MMR', bty = 'n', main = 'CpG > TpG')
for (j in 1:length(mean_c_to_t)) {
  lines(y = c(mmr_mean_c_to_t[[j]] - 1.96*mmr_mean_c_to_t_sd[[j]],mmr_mean_c_to_t[[j]] + 1.96*mmr_mean_c_to_t_sd[[j]]),
        x = c(mean_c_to_t[[j]],mean_c_to_t[[j]]), col = 'red', lwd = 1)
  lines(x = c(mean_c_to_t[[j]] - 1.96*mean_c_to_t_sd[[j]],mean_c_to_t[[j]] + 1.96*mean_c_to_t_sd[[j]]),
        y = c(mmr_mean_c_to_t[[j]],mmr_mean_c_to_t[[j]]), col = 'red', lwd = 1)
}
#lm(unlist(mmr_mean_c_to_t) ~ unlist(mean_c_to_t))
#abline(a = 0, b = 3.2, lty = 2)
#text(x = 6, y = 12, expression(y == -0.76 + 3.2 * x))
text(x = unlist(mean_c_to_t), y = unlist(mmr_mean_c_to_t), labels = names(unlist(mean_c_to_t)))
dev.off()

pdf('MMR_DEL_scatterplot.pdf',4,4)
plot(y = unlist(mmr_mean_del), 
     x = unlist(mean_del), pch = 16,col='orange',
     xlab = 'Non-MMR', ylab = 'MMR', bty = 'n', main = '1bp deletions')
for (j in 1:length(mean_del)) {
  lines(y = c(mmr_mean_del[[j]] - 1.96*mmr_mean_del_sd[[j]],mmr_mean_del[[j]] + 1.96*mmr_mean_del_sd[[j]]),
        x = c(mean_del[[j]],mean_del[[j]]), col = 'orange', lwd = 2)
  lines(x = c(mean_del[[j]] - 1.96*mean_del_sd[[j]],mean_del[[j]] + 1.96*mean_del_sd[[j]]),
        y = c(mmr_mean_del[[j]],mmr_mean_del[[j]]), col = 'orange', lwd = 1.5)
}
#lm(unlist(mmr_mean_del) ~ 0 + unlist(mean_del))
#abline(a = 0, b = 5.5, lty = 2)
text(x = unlist(mean_del), y = unlist(mmr_mean_del), labels = names(unlist(mean_del)))
#text(x = 1.2, y = 4, expression(y == -0.5 + 7.3 * x))
dev.off()

pdf('MMR_INS_scatterplot.pdf',4,4)
plot(y = unlist(mmr_mean_ins), 
     x = unlist(mean_ins), pch = 16,col='purple',
     xlab = 'Non-MMR', ylab = 'MMR', bty = 'n', main = '1 bp insertions')
for (j in 1:length(mean_ins)) {
  lines(y = c(mmr_mean_ins[[j]] - 1.96*mmr_mean_ins_sd[[j]],mmr_mean_ins[[j]] + 1.96*mmr_mean_ins_sd[[j]]),
        x = c(mean_ins[[j]],mean_ins[[j]]), col = 'purple', lwd = 2)
  lines(x = c(mean_ins[[j]] - 1.96*mean_ins_sd[[j]],mean_ins[[j]] + 1.96*mean_ins_sd[[j]]),
        y = c(mmr_mean_ins[[j]],mmr_mean_ins[[j]]), col = 'purple', lwd = 1.5)
}
#lm(unlist(mmr_mean_ins) ~ 0 + unlist(mean_ins))
#abline(a = 0, b = 4, lty = 2)
text(x = unlist(mean_ins), y = unlist(mmr_mean_ins), labels = names(unlist(mean_ins)))
#text(x = 0.4, y = 1, expression(y == 4 * x))
dev.off()




abline(a = 0, b = 1)
abline(a = 0, b = 10,col='red')

mean_prof <- do.call('rbind',mean_prof)
mean_sd <- do.call('rbind',mean_sd)

plot_subindel_wb(t(mean_prof))

to.show.mean.sd <- sqrt(rbind(rowSums(mean_sd[,c(35,39,43,47)]**2),rbind(t(mean_sd[,c(97,102)]**2))))

to.show.new_var_tmp <- sqrt(rbind(colSums(new_var_tmp[c(35,39,43,47),]**2),rbind(new_var_tmp[c(97,102),]**2)))

to.show.mean <- rbind(rowSums(mean_prof[,c(35,39,43,47)]),rbind(t(mean_prof[,c(97,102)])))

to.show.altered_sig <- rbind(colSums(new_altered_signatures[c(35,39,43,47),]),rbind(new_altered_signatures[c(97,102),]))

show.sd <- sqrt(1 / to.show.mean**2 * to.show.new_var_tmp**2 + to.show.altered_sig**2 / to.show.mean**4 * to.show.mean.sd**2)

to.show <- rbind(colSums(new_altered_signatures[c(35,39,43,47),]),rbind(new_altered_signatures[c(97,102),])) / 
  rbind(rowSums(mean_prof[,c(35,39,43,47)]),rbind(t(mean_prof[,c(97,102)])))

to.show.low <- to.show - 1.96 * show.sd
to.show.high <- to.show + 1.96 * show.sd

colnames(to.show) = colnames(to.show.low) = colnames(to.show.high) <- colnames(altered_signatures)



df.mean = data.frame(name = colnames(to.show),
                     mean = as.numeric(to.show[3,]),
                     low = as.numeric(to.show.low[3,]),
                     high = as.numeric(to.show.high[3,]))
X_axis <- c(1:nrow(df.mean))
pdf('mmr_plot_ins.pdf',3,4)
par(mar = c(8,4,4,4))
plot(x =  c(0.5,1:(nrow(df.mean)+1)), xaxt = 'n', y = rep(NA,length(X_axis)+2), ylim = c(-1,5),
     yaxt = 'n', xlab = '', ylab = '', bty='n', main = 'Foldchange in average number of C>T mutations')
for (j in 1:length(X_axis)) {
  cur.col <- 'gray88'; cur.lwd = 1
  if (df.mean$low[j] > 1) {
    cur.col <- 'gray48'
    cur.lwd <- 2
  }
  polygon(x = c(X_axis[j] - 0.5,rep(X_axis[j] + 0.5,2),rep(X_axis[j] - 0.5,2)),
          y = c(rep(df.mean$high[j],2),rep(df.mean$low[j],2),df.mean$high[j]),
          col = cur.col,
          border = 0.01)
  lines(x = c(X_axis[j] - 0.5, X_axis[j] + 0.5), y = rep(df.mean$mean[j],2), lwd = cur.lwd)
}
abline(h = 0, lty = 2)
#abline(h = quantile(log10(av.changes), 0.05), col = 'red', lty = 2)
#abline(h = quantile(log10(av.changes), 0.95), col = 'red', lty = 2)

axis(side = 2, labels = c(0:5), at = c(0:5),las=2)
axis(side = 1, 
     labels = df.mean$name,
     at = X_axis,
     las=2, col ='white')
dev.off()



pdf('MMR_change_2904.pdf',10,5)
par(mar = c(6,3,2,2))
fff <- barplot(as.matrix(to.show),ylim = c(0,10), ylab = 'Fold-change in contribution compared to non-MMR',
               beside = T, las = 2, col = c('red',"orange","purple"), border = NA)
legend('top',legend = c('C>T at CpG sites','Deletions','Insertions'), fill = c('red',"orange","purple"), bty = 'n')
arrows(x0 = as.vector(fff),
       y0 = as.vector(as.matrix(to.show.low)),
       y1 = as.vector(as.matrix(to.show.high)),
       col = 'darkgrey', length = 0.05, angle = 90, lwd = 2)
arrows(x0 = as.vector(fff),
       y1 = as.vector(as.matrix(to.show.low)),
       y0 = as.vector(as.matrix(to.show.high)),
       col = 'darkgrey', length = 0.05, angle = 90, lwd = 2)
abline(h = 1, lty = 2)
dev.off()


pdf('MMR_signature.pdf',10,4)
tmp <- data.frame(mmr = rowSums(new_altered_signatures))
tmp.sd <- data.frame(mmr = rowSums(new_var_tmp))
tmp.low <- tmp - 1.96*tmp.sd
tmp.low[tmp.low<0] <- 0
plot_subindel_wb(tmp, CI = T, low = tmp.low, high = tmp + 1.96*tmp.sd) + theme(panel.grid = element_blank())
dev.off()



#####################################################

# check differences per gene?


# Labels
mmr.genes <- c('MLH1','MLH3','PMS2','MSH2','MSH3','MSH6')

# Expression
metadata <- read.table('TCGA_caveman_patients_to_files.dat', sep='\t', header=T)
check.exp <- function(z, tr, m = metadata) {
  return(which(m[,match(z,colnames(m))] < tr))
}
for (i in 6:ncol(metadata)) {
  metadata[,i] <- metadata[,i] / median(metadata[,i], na.rm = T)
  #  metadata[,i] <- rank(metadata[,i], na.last='keep') / length(metadata[,i])
}
exp.mmr <- lapply(mmr.genes, function(x) as.character(metadata$tumour[sort(check.exp(x,0.2))]))
names(exp.mmr) <- mmr.genes

# Germline

SNPs <- list()
tmp <- read.table('~/TCGAmutations/results/germline_MMR_38.txt', sep = '\t', header = T)
for (x in mmr.genes) {
  SNPs[[x]] <- unique(as.character(tmp$Uploaded_variation[tmp$IMPACT %in% c('HIGH') & grepl(x, tmp$SYMBOL)]))
}
SNPs <- SNPs[sapply(SNPs,length)>0]

intsam <- list(); z = 0
for (cancername in unlist(strsplit(as.character(unique(metadata$project)),split = ' '))[-c(1,5)]) {
  coad <- as.character(metadata$tumour[metadata$project==cancername])
  for (genename in names(SNPs)) {
    f <- list.files(paste0('/Volumes/SNP1/',genename,'/',cancername))
    for (i in 1:length(coad)) {
      file <- f[grep(coad[i],f)]
      if (length(file)==0) next
      if (file.info(paste0('/Volumes/SNP1/',genename,'/',cancername,'/',file))$size==0) next
      vcf <- read.delim(paste0('/Volumes/SNP1/',genename,'/',cancername,'/',file), header = F)
      identifier <- paste0(vcf$V1,':',vcf$V2,'_',vcf$V4,'/',vcf$V5) 
      het <- substr(vcf$V10,1,3)[identifier %in% SNPs[[genename]]]
      if (length(het)>0) {
        z = z + 1
        print(c(as.character(coad[i]),het,genename))
        intsam[[z]] <- c(as.character(coad[i]),het,genename)
      }
    }
    print(genename)
  }
}
sapply(intsam, function(x) x[1])

# Damaging mutations
# Somatic

copynumber <- read.delim('filtered.combined.segments.txt', sep = '\t', header = T)
copynumber <- GRanges(copynumber)
ids <- read.table('~/Downloads/ICGCtoTCGA.tsv', sep = '\t', header = T)
samples <- list()
for (x in mmr.genes) {
  tmp <- read.delim(paste0('~/TCGAmutations/results/',x,'.txt'), sep='\t', header=T)
  tmp <- tmp[tmp$IMPACT %in% c('HIGH') & grepl(x, tmp$SYMBOL),]
  if (nrow(tmp)==0) next
  mutations <- read.delim(paste0('/Volumes/nobackup/MMR/',x,'.dat'), sep='\t', header = T)
  identifier <- paste0(mutations$Sample,'-',mutations$Chrom,':',mutations$Pos,'_',mutations$Ref,'/',mutations$Alt)
  a <- sapply(as.character(tmp$Uploaded_variation), function(y) substr(unlist(strsplit(y,split='[:]'))[1],1,12))
  vaf <- mutations$PM.Tum[match(as.character(tmp$Uploaded_variation), identifier)]
  samples[[x]] <- paste0('TCGA', substr(a,5,12))[vaf>0.3]
  samples[[x]] <- unlist(sapply(samples[[x]], function(y) metadata$tumour[grep(y,metadata$patient_id)]))
  cn <- copynumber[queryHits(findOverlaps(copynumber, genes_of_interest[genes_of_interest$name==x]))]
  cn <- cn[cn$tissue %in% names(samples[[x]])]
  samples[[paste0(x,'LOH')]] <- samples[[x]][as.character(cn$tissue[cn$minor==0])]
}

save(samples, intsam, exp.mmr, file = 'nb/MMR/data_on_TCGA_MMR_samples_intsam_exp.mmr.RData')
samples -> samples.mmr 
new.mmr.2 <- unique(c(as.character(unlist(exp.mmr)),
                    as.character(unlist(samples.mmr)),
                    as.character(sapply(intsam, function(x) x[1]))))


new.mmr.2 <- unique(c(as.character(unlist(exp.mmr)),
                      rownames(bigmat)[which(substr(alternative_rownames_bigmat,1,12) %in% as.character(unlist(samples[paste0(MMR.core,'_hom')])))],
                      as.character(sapply(intsam, function(x) x[1]))))

length(intersect(substr(new.mmr.2,1,12), unlist(samples[paste0(MMR.core,'_hom')])))
length(intersect(substr(new.mmr.2,1,12), unlist(samples[paste0(MMR.core,'_het')])))

new.mmr.2 <- new.mmr.2[new.mmr.2 %in% alternative_rownames_bigmat]
new.mmr.2 <- new.mmr.2[rowSums(bigmat[match(new.mmr.2,alternative_rownames_bigmat),1:96])>100] # 263

per.gene <- list()
for (y in mmr.genes) {
  per.gene[[y]] <- unique(c(as.character(unlist(exp.mmr[[y]])),
                   as.character(unlist(samples[[y]])),
                   as.character(unlist(samples[[paste0(y,'LOH')]])),
                   as.character(sapply(intsam[sapply(intsam, function(x) x[length(x)])==y], function(x) x[1]))))
}


table(metadata$project[match(new.mmr.2, metadata$tumour)])

plot_subindel_wb(t(bigmat[new.mmr.2[1:20],]))

# include tissue, gene, age (too many?)

new.mmr.2 <- new.mmr.2[sapply(new.mmr.2, function(x) cosine(as.numeric(bigmat[x,1:96]), cancer_signatures[,13])) < 0.75 &
            sapply(new.mmr.2, function(x) cosine(as.numeric(bigmat[x,1:96]), cancer_signatures[,2])) < 0.75] # 548


#tmp <- read.delim('/Volumes/r/nadia/tcga_uterus/POLE.dat', sep='\t', header = T)
#tmp <- tmp[as.character(tmp$Effect)=='missense' & tmp$Gene=='POLE' & nchar(as.character(tmp$Protein))==7,]
#tmp <- tmp[as.numeric(substr(as.character(tmp$Protein),4,6)) < 472 & 
#             as.numeric(substr(as.character(tmp$Protein),4,6)) > 267,]
#data$POLE <- gsub(ucec,pattern = 'TCGA',replacement = 'H_LR') %in% as.character(tmp$Sample)
#tmp <- tmp[as.character(tmp$Sample) %in% gsub(ucec,pattern = 'TCGA',replacement = 'H_LR'),]
#data$whichPOLE <- tmp$Protein[match(gsub(ucec,pattern = 'TCGA',replacement = 'H_LR'), as.character(tmp$Sample))]


mmr_data_2 <- list()
mmr_data_2[['y']] <- bigmat[match(new.mmr.2,alternative_rownames_bigmat),]
mmr_data_2[['X']] <- model.matrix( ~ .,data = data.frame(v1 = as.character(metadata$project[match(new.mmr.2,metadata$tumour)])))[,-1]
mmr_data_2[['X']] <- cbind(mmr_data_2[['X']], BRAIN = as.numeric((mmr_data_2[['X']][,'v1GBM'] + mmr_data_2[['X']][,'v1LGG']) > 0),
                           BREAST = as.numeric(mmr_data_2[['X']][,'v1BRCA'] > 0),
                         STES = as.numeric((mmr_data_2[['X']][,'v1ESCA'] + mmr_data_2[['X']][,'v1STAD']) > 0),
                         COLORECT = as.numeric((mmr_data_2[['X']][,'v1COAD'] + mmr_data_2[['X']][,'v1READ']) > 0),
                         KIPC = as.numeric((mmr_data_2[['X']][,'v1KIRC'] + mmr_data_2[['X']][,'v1KIRP'] + mmr_data_2[['X']][,'v1KICH']) > 0),
                         OV = as.numeric(mmr_data_2[['X']][,'v1OV'] > 0),
                         LUNG = as.numeric((mmr_data_2[['X']][,'v1LUSC'] + mmr_data_2[['X']][,'v1LUAD']) > 0),
                         UTERUS = as.numeric(mmr_data_2[['X']][,'v1UCEC'] > 0))
mmr_data_2[['X']] <- mmr_data_2[['X']][,-c(1:25)]

mmr_data_2[['X']] <- mmr_data_2[['X']][,colSums(mmr_data_2[['X']])>10]
mmr_data_2[['y']] <- mmr_data_2[['y']][rowSums(mmr_data_2[['X']])>0,]
mmr_data_2[['X']] <- mmr_data_2[['X']][rowSums(mmr_data_2[['X']])>0,]

#mmr_data_2[['X']] <- cbind(MLH1 = as.numeric(rownames(mmr_data_2[['y']]) %in% per.gene[['MLH1']]),
#                           MLH3 = as.numeric(rownames(mmr_data_2[['y']]) %in% per.gene[['MLH3']]),
#                           PMS2 = as.numeric(rownames(mmr_data_2[['y']]) %in% per.gene[['PMS2']]),
#                           MSH2 = as.numeric(rownames(mmr_data_2[['y']]) %in% per.gene[['MSH2']]),
#                           MSH3 = as.numeric(rownames(mmr_data_2[['y']]) %in% per.gene[['MSH3']]),
#                           MSH6 = as.numeric(rownames(mmr_data_2[['y']]) %in% per.gene[['MSH6']]))

#mmr_data_2[['X']] <- mmr_data_2[['X']][,colSums(mmr_data_2[['X']])>10]
#mmr_data_2[['y']] <- mmr_data_2[['y']][rowSums(mmr_data_2[['X']])>0,]
#mmr_data_2[['X']] <- mmr_data_2[['X']][rowSums(mmr_data_2[['X']])>0,]

#colnames(mmr_data_2[['X']])[1:5] <- substr(colnames(mmr_data_2[['X']])[1:5],3,nchar(colnames(mmr_data_2[['X']])[1:5]))
#mmr_data_2[['X']] <- mmr_data_2[['X']][,-match('COLORECT', colnames(mmr_data_2[['X']]))]


#mmr_data_2[['X']] <- cbind(mmr_data_2[['X']],AGE = ids$donor_age_at_diagnosis[
#  match(metadata$patient_id[match(rownames(mmr_data_2[['y']]),
#                                  metadata$tumour)], 
#        ids$submitted_donor_id)])
#mmr_data_2[['y']] <- mmr_data_2[['y']][!is.na(mmr_data_2[['X']][,'AGE']),] 
#mmr_data_2[['X']] <- mmr_data_2[['X']][!is.na(mmr_data_2[['X']][,'AGE']),] 
#mmr_data_2[['X']] <- mmr_data_2[['X']][,colSums(mmr_data_2[['X']])>10]
#mmr_data_2[['y']] <- mmr_data_2[['y']][rowSums(mmr_data_2[['X']])>0,]
#mmr_data_2[['X']] <- mmr_data_2[['X']][rowSums(mmr_data_2[['X']])>0,]
#mmr_data_2[['X']][,'AGE'] <- scale(mmr_data_2[['X']][,'AGE'], center = T)


mmr_data_2[['N']] <- nrow(mmr_data_2[['y']]) # 482
mmr_data_2[['R']] <- ncol(mmr_data_2[['y']]) # 104
mmr_data_2[['K']] <- ncol(mmr_data_2[['X']]) # 10
mmr_data_2[['S']] <- 3
mmr_data_2[['M']] <- max(rowSums(mmr_data_2[['y']]))

library(greta)
S <- variable(lower = 0, upper = 1, dim = c(mmr_data_2[['S']],mmr_data_2[['R']]))
S <- S / (greta::rowSums(S) %*% matrix(1,nrow = 1, ncol = mmr_data_2[['R']]))
E <- variable(lower = 0, upper = mmr_data_2[['M']], dim = c(mmr_data_2[['N']],mmr_data_2[['S']]))

sigma = 0.1
beta = normal(mean = 0, sd = sigma, dim = c(mmr_data_2[['K']],mmr_data_2[['R']]))

size = 50

#mu1 = E %*% S
mu1 = E[,-(mmr_data_2[['S']]),drop=F] %*% S[-(mmr_data_2[['S']]),,drop=F] + 
  (E[,mmr_data_2[['S']],drop=F] %*% S[mmr_data_2[['S']],,drop = F]) * exp(mmr_data_2[['X']] %*% beta) / 
((exp(mmr_data_2[['X']] %*% beta) %*% t(S[mmr_data_2[['S']],,drop=F])) %*% matrix(1,nrow = 1, ncol=104))

prob = size/(size + mu1)
distribution(mmr_data_2[['y']]) = negative_binomial(size = matrix(1, nrow = mmr_data_2[['N']], ncol = mmr_data_2[['R']]) * size, prob = prob)

#m <- model(S,E)
m <- model(S,E,beta)

# sampling
draws <- mcmc(m,n_samples = 1000, warmup = 1000,verbose=F)
draws_2 <- mcmc(m, n_samples = 1000, warmup = 1000, verbose = F)
draws_3 <- mcmc(m, n_samples = 1000, warmup = 1000, verbose = F)
draws_4 <- mcmc(m, n_samples = 1000, warmup = 1000, verbose = F)

draws[[2]] <- draws_2[[1]]
draws[[3]] <- draws_3[[1]]
draws[[4]] <- draws_4[[1]]

save(draws,'/Volumes/nobackup/MMR/mmr_2_draws_23012019.RData')

mcmc_trace(draws[,grep('S',colnames(draws[[1]]))[sample(1:(mmr_data_2[['S']]*104),9)]])
mcmc_trace(draws[,grep('beta',colnames(draws[[1]]))[sample(1:(mmr_data_2[['K']]*104),9)]])

draws_all <- draws[[1]] # 3
S_est <- matrix(colMeans(draws_all[,grep('S',colnames(draws_all), fixed = T)]), nrow = mmr_data_2[['S']], ncol = mmr_data_2[['R']])
S_low <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = mmr_data_2[['S']], ncol = mmr_data_2[['R']])
S_high <- matrix(apply(draws_all[,grep('S',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = mmr_data_2[['S']], ncol = mmr_data_2[['R']])
rownames(S_est) = rownames(S_low) = rownames(S_high) = paste0('Sig',1:mmr_data_2[['S']])
E_est <- matrix(colMeans(draws_all[,grep('E',colnames(draws_all))]), ncol = mmr_data_2[['S']])
beta_est <- matrix(colMeans(draws_all[,grep('beta',colnames(draws_all), fixed = T)]), nrow = mmr_data_2[['K']], ncol = mmr_data_2[['R']])
beta_low <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.025), nrow = mmr_data_2[['K']], ncol = mmr_data_2[['R']])
beta_high <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,quantile,0.975), nrow = mmr_data_2[['K']], ncol = mmr_data_2[['R']])
beta_var <- matrix(apply(draws_all[,grep('beta',colnames(draws_all), fixed = T)],2,var), nrow = mmr_data_2[['K']], ncol = mmr_data_2[['R']])

#plot_subindel_wb(t(S_est), CI = T, low = t(S_low), high = t(S_high))
plot_subindel_wb(t(S_est), CI = T, low = t(S_low), high = t(S_high))
plot_subindel_wb(data.frame(t(S_est),
                       S_est[mmr_data_2[['S']],] * exp(beta_est[1,]),
                       S_est[mmr_data_2[['S']],] * exp(beta_est[2,]),
                       S_est[mmr_data_2[['S']],] * exp(beta_est[3,]),
                       S_est[mmr_data_2[['S']],] * exp(beta_est[4,]),
                       S_est[mmr_data_2[['S']],] * exp(beta_est[5,]),
                       S_est[mmr_data_2[['S']],] * exp(beta_est[6,]),
                       S_est[mmr_data_2[['S']],] * exp(beta_est[7,]),
                       S_est[mmr_data_2[['S']],] * exp(beta_est[8,])))#, 
#                       STES = S_est[mmr_data_2[['S']],] * exp(beta_est[7,]), 
#                       KIDNEY = S_est[mmr_data_2[['S']],] * exp(beta_est[8,]), 
#                       LUNG = S_est[mmr_data_2[['S']],] * exp(beta_est[9,])))

par(mfrow = c(5,2))
for (j in 1:10)
  interaction_effect_plot_human(beta_est[j,], lwd = 2, CI = T,
                                low = beta_low[j,],
                                high = beta_high[j,],
                                at = c(-1,0,1),
                                labels = c('<0.1',1,10), plot_main = colnames(mmr_data_2[['X']])[j])

mu <- E_est[,-(mmr_data_2[['S']]),drop=F] %*% S_est[-(mmr_data_2[['S']]),,drop=F] + 
  (E_est[,mmr_data_2[['S']],drop=F] %*% S_est[mmr_data_2[['S']],,drop = F]) * exp(mmr_data_2[['X']] %*% beta_est)

plot(rowSums(mmr_data_2[['y']]),rowSums(mu), pch = 16,
     xlab = 'Observed', ylab= 'Predicted')
abline(a=0,b=1,col='red',lty=2)

#save(draws, file = '/Volumes/nobackup/MMR_draws_061218_mutated_samples_cancer_types_only.RData')
save(draws, file = 'nb/MMR_draws_071218_mutated_samples_cancer_types_and_age_only.RData')


pdf('~/MMR_observed_predicted.pdf',5,5)
plot(rowSums(mmr_data_2[['y']]),rowSums(mu), pch = 16,
     xlab = 'Observed', ylab= 'Predicted')
abline(a=0,b=1,col='red',lty=2)
dev.off()

# Make a tsne plot
library(tsne)
mut_mat <- t(mmr_data_2[['y']])
cosdist <- function(x,y) {
  x0 <- x/sum(x)
  y0 <- y/sum(y)
  x0 %*% y0 / sqrt(x0%*%x0)/sqrt(y0%*%y0)
}
D <- as.dist(sapply(1:ncol(mut_mat), function(i) sapply(1:ncol(mut_mat), function(j) 1-cosdist(mut_mat[,i],mut_mat[,j]))))
set.seed(123)
t <- tsne(D, perplexity = 20)
rownames(t) <- colnames(mut_mat)

# Get the MMR and POLE status

decomposition <- E_est
colnames(decomposition) <- c('POLE','BRCA','26','MMR')
for (i in 1:nrow(decomposition))
  decomposition[i,] <- decomposition[i,] / sum(decomposition[i,])

# Visualize it
library(RColorBrewer)
source('~/MMR/plotting functions/scatterpie.R')
library(mg14)
pdf('~/MMR_tsne.pdf', 10, 15)
#par(bty="n", mar=c(0,0,0,0))
cancers <- colnames(mmr_data_2[['X']])[-11]
plot(NA,NA, xlab="", ylab="", xlim=c(-33,33), ylim=c(-42,55), xaxt="n", yaxt="n", bty='n')
for (j in 1:length(cancers)) {
  sel <- (mmr_data_2[['X']][,cancers[j]]>0)
  o1 <- order(rowSums(mmr_data_2[['y']])[sel],decreasing = T)
  corr_scatterpie(t[sel,1][o1], t[sel,2][o1], p=decomposition[sel,][o1,],
                  r=sqrt(colSums(mut_mat)[sel][o1])/30, labels=NA, col=brewer.pal(4,'Set1'),
                  lty=0, circles=TRUE, lwd.circle=rep(2.5,sum(sel)),
                  lty.circle=rep(1,sum(sel)), add=TRUE, col.circle = brewer.pal(11,'Set3')[j])
}
mg14:::.pie(x0=-30, y0=50, x=matrix(rep(1,4), nrow=1), r=sqrt(10000)/20, labels=colnames(decomposition),
            col=brewer.pal(4,'Set1'), lty=0, circles=TRUE, add=TRUE, cex = 0.8)
us <- par("usr")
pr <- (us[2]-us[1])/(us[4]-us[3])
fr <- par("pin")[1]/par("pin")[2]
#for(i in c(1,10,100,1000,10000)){
#  polygon(-20 + cos(seq(0,2*pi, l=100)) * sqrt(i)/75, 15+(1+sin(seq(0,2*pi, l=100))) * sqrt(i)/75 / pr * fr, col=NA)
#  if (i>10) text(-20, 15 + 2*sqrt(i)/75 / pr * fr + 0.3,labels = as.character(i),cex=0.8)
#}
for (j in 1:length(cancers)) {
  print(polygon(0 + cos(seq(0,2*pi, l=100)) , 50-j*3 + (1+sin(seq(0,2*pi, l=100))) / pr * fr,
                lwd=0.5, col=brewer.pal(11,'Set3')[j]), border = 'black')
  print(text(x = 15,y=50-j*3+1,labels = cancers[j]))
}
dev.off()

pdf('MMR_signature.pdf',12,4)
tmp <- data.frame(t(S_est[4,,drop = F]))
plot_subindel_wb(tmp, CI = T, low = t(S_low[4,,drop=F]), high = t(S_high[4,,drop=F])) + theme(panel.grid = element_blank())
dev.off()
pdf('MMR_signatures.pdf',12,8)
rownames(S_est) = rownames(S_high) = rownames(S_low) <- c('Sig26', 'POLEMMR','BRCA','MMR')
plot_subindel_wb(t(S_est), CI = T, low = t(S_low), high = t(S_high)) + theme(panel.grid = element_blank())
dev.off()
pdf('MMR_alterations.pdf',12,20)
tmp <- data.frame(cbind(S_est[mmr_data_2[['S']],] * exp(beta_est[1,]),
                        S_est[mmr_data_2[['S']],] * exp(beta_est[2,]),
                        S_est[mmr_data_2[['S']],] * exp(beta_est[3,]),
                        S_est[mmr_data_2[['S']],] * exp(beta_est[4,]),
                        S_est[mmr_data_2[['S']],] * exp(beta_est[5,]),
                        S_est[mmr_data_2[['S']],] * exp(beta_est[6,]),
                        S_est[mmr_data_2[['S']],] * exp(beta_est[7,]),
                        S_est[mmr_data_2[['S']],] * exp(beta_est[8,]),
                        S_est[mmr_data_2[['S']],] * exp(beta_est[9,]),
                        S_est[mmr_data_2[['S']],] * exp(beta_est[10,]),
                        S_est[mmr_data_2[['S']],] * exp(beta_est[11,])))
tmp_low <- data.frame(t(S_low))
tmp_high <- data.frame(t(S_high))
colnames(tmp) = colnames(tmp_high) = colnames(tmp_low) <- c('Sig26', 'BRCA','POLE','MMR')
plot_subindel_wb(t(S_est), CI = T, low = tmp_low, high = tmp_high) + theme(panel.grid = element_blank())
dev.off()

pdf('MMR_change.pdf',12,5)
rownames(beta_est) <- colnames(mmr_data_2[['X']])
fff <- barplot(t(exp(beta_est[-12,c(97,102)])), beside = T, las = 2, legend.text = c('Deletions','Insertions'), col = c("orange","purple"), border = NA)
arrows(fff, y0 = t(exp(beta_low[-12,c(97,102)])), y1 = t(exp(beta_high[-12,c(97,102)])), col = 'darkgrey', length = 0.01, angle = 90, lwd = 2)
abline(h = 1, lty = 2)
dev.off()

pdf('MMR_effects.pdf', 12,12)
par(mfrow = c(5,2))
for (i in 1:10)
  interaction_effect_plot_human(beta_est[i,], CI = T, low = beta_low[i,], high = beta_high[i,],
                                plot_main = cancers[i], at = c(-1,0,1), labels = c('<0.1', 1, 10))
dev.off()

########################################################################

pcawg.msi <- read.xlsx('~/Downloads/media-2.xlsx', sheetIndex = 3, startRow = 2)
pcawg.msi <- pcawg.msi[,1:6]
pcawg.msi$proportion <- pcawg.msi$Number.of.mutations.in.informative.microsatellite / pcawg.msi$Total.number.of.analyzed.informative.microsatellite
mmr.pcawg <- as.character(pcawg.msi$ID)[pcawg.msi$proportion >= 0.03]

########################################################################


