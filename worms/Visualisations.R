# Visualizing all the WORM signatures

source('../useful_functions.R')
source('../plotting_functions.R')
library(reshape2)
library(ggplot2)
library(ggpubr)
clrs <- c("#2EBAED","#000000","#DE1C14","#D4D2D2","#ADCC54","#F0D0CE","brown","#8DD3C7","#FFFFB3","#BEBADA","darkmagenta")

################################ Supplementary figures 1 & 2 ###################################

# Run the Worm_model.R first and have all the draws and coefficient ready

# Genotype signatures
no_germline <- c('exo.1','agt.1','rad.51')
q <- list()
j = 1
for (gene in colnames(G1)) {
  ym = max(beta_GH_greta_full_high[,gene,drop=F])
  if (ym < 0.5) ym = 0.5
  proper_name <- data$Genotype[match(rownames(G1)[G1[,gene]>0],data$Sample)][1]
  max_generation <- max(as.numeric(data$Generation)[match(rownames(G1)[G1[,gene]>0],data$Sample)])
  title1 <- paste0('Effects for ',proper_name, '; ', round(sum(beta_GH_greta_full[,gene,drop=F]),1), ' (',
                   round(sum(beta_GH_greta_full_low[,gene,drop=F]),1),'-', round(sum(beta_GH_greta_full_high[,gene,drop=F]),1),
                   ') het. mut-s per gen., max ', max_generation, ' gen. with ',
                   round(mean(rowSums(Y)[names(CD2Mutant)[CD2Mutant==paste0(proper_name,':',max_generation)]])), ' (',
                   min(rowSums(Y)[names(CD2Mutant)[CD2Mutant==paste0(proper_name,':',max_generation)]]), '-',
                   max(rowSums(Y)[names(CD2Mutant)[CD2Mutant==paste0(proper_name,':',max_generation)]]), ') mut-s')
  if (gene %in% no_germline)
    title1 <- paste0('Estimated effects for ',proper_name,'; ', round(sum(beta_GH_greta_full[,gene,drop=F]),1), ' (',
                     round(sum(beta_GH_greta_full_low[,gene,drop=F]),1),'-', round(sum(beta_GH_greta_full_high[,gene,drop=F]),1),
                     ') het. mut-s per gen., max 1 generation')
  if (j <52)
    q[[j]] <- plot_dnvhugesig_wb(beta_GH_greta_full[,gene,drop=F], CI = T,
                                 low = beta_GH_greta_full_low[,gene,drop=F], ymax = ym,
                                 high = beta_GH_greta_full_high[,gene,drop=F], norm = F,
                                 colors = c("#2EBAED","#000000","#DE1C14","#D4D2D2","#ADCC54","#F0D0CE",
                                            "brown","#8DD3C7","#FFFFB3","#BEBADA","darkmagenta"),
                                 rownames = F) + 
      theme(panel.grid = element_blank(), strip.text.y = element_blank()) + ggtitle(title1)
  else q[[j]] <- plot_dnvhugesig_wb(beta_GH_greta_full[,gene,drop=F], CI = T,
                                    low = beta_GH_greta_full_low[,gene,drop=F], ymax = ym,
                                    high = beta_GH_greta_full_high[,gene,drop=F], norm = F,
                                    colors = c("#2EBAED","#000000","#DE1C14","#D4D2D2","#ADCC54","#F0D0CE",
                                               "brown","#8DD3C7","#FFFFB3","#BEBADA","darkmagenta")) + 
    theme(panel.grid = element_blank(), strip.text.y = element_blank()) + ggtitle(title1)
  j = j + 1
}
pdf('SupplementaryFigure1.pdf',35,50)
ggarrange(plotlist = q, ncol = 3, nrow = 18)
dev.off()


# Mutagen signatures
unit <- c(1, 10, 10, 100, 0.1, 100, 1, 10, 1, 100, 10, 1)
#  1 microM, 10 microM, 10 microM, 100 Gray, 100 microM, 100 Gray, 1 mM, 10 milliM, 1 mM, 100 Joule, 10 microM, 1 mM
unit_name <- c('muM','muM','muM','Gy', 'muM', 'Gy', 'mM', 'mM', 'mM','J/m2', 'muM','mM')
names(unit_name) = names(unit) = names(avdose) <- colnames(Mall)
library(ggpubr)
pdf('SupplementaryFigure2.pdf',30,20)
q <- list()
for (i in 1:r) {
  ym = max(beta_M_greta_full_high[,i,drop=F])
  if (ym < 0.5) ym = 0.5
  q[[i]] <- plot_dnvhugesig_wb(beta_M_greta_full[,i,drop=F] / avdose[i] * unit[i],
                               CI = T, low = beta_M_greta_full_low[,i,drop=F] / avdose[i] * unit[i],
                               high = beta_M_greta_full_high[,i,drop=F] / avdose[i] * unit[i],
                               norm = F, ymax = max(beta_M_greta_full_high[,i,drop=F] / avdose[i] * unit[i]),
                               colors = c("#2EBAED","#000000","#DE1C14","#D4D2D2","#ADCC54","#F0D0CE","brown","#8DD3C7","#FFFFB3","#BEBADA","darkmagenta"),
                               rownames = F) + 
    theme(panel.grid = element_blank(), strip.text.y = element_blank()) + 
    ggtitle(paste0('Effects for ',colnames(beta_M_greta_full)[i], ', ',
                   round(sum(beta_M_greta_full[,i]/ avdose[i] * unit[i]),1), ' mutations on average per ',unit[i],' ',unit_name[i]))
}
ggarrange(plotlist = q, ncol = 3, nrow = 4)
dev.off()

# More summary plots for Figure 2

plot_dnvhugesig_wb(mut_matrix = sapply(colnames(beta_M_greta_full)[c(1,11,3,8,7,4,10)], function(i) 
                                          beta_M_greta_full[,i] / avdose[i] * unit[i]), 
                   colors = clrs, CI = T,
                   low = sapply(colnames(beta_M_greta_full)[c(1,11,3,8,7,4,10)], function(i) 
                     beta_M_greta_full_low[,i] / avdose[i] * unit[i]),
                   high = sapply(colnames(beta_M_greta_full)[c(1,11,3,8,7,4,10)], function(i) 
                     beta_M_greta_full_high[,i] / avdose[i] * unit[i]), 
                   norm = F, diff_scale = T, diff_limits = c(1,2,0.5,30,20,3,1)) +
  theme(panel.grid = element_blank())

inds <- c('agt.1','mlh.1','polh.1','rev.3','smc.6','xpc.1')
g <- plot_dnvhugesig_wb(beta_GH_greta_full[,inds], colors = clrs,CI = T,
                   low = (beta_GH_greta_full_low[,inds]), 
                   high = (beta_GH_greta_full_high[,inds]),
                   diff_scale = T, norm = F, diff_limits = c(0.5,50,1,1,0.5,0.5)) + 
  theme(panel.grid = element_blank())


#########################################################################################################

# Interaction effects #

for (z in c('agt.1.MMS','polk.1.MMS','agt.1.EMS','polk.1.EMS','xpc.1.UV','xpf.1.AristolochicAcid','rev.3.UV','polh.1.EMS')) {
  pdf(paste0(z,'.beta.I.greta.pdf'),12,3)
  par(mar = c(2, 4, 4, 2) + 0.1)
  current_at <- c(-1,0,1,2)
  if (z == 'xpc.1.UV' || x == 'xpf.1.UV') current_at <- c(-1,0,1,2,3)
  interaction_effect_plot(beta_I_greta_full[,z], CI = T, at = current_at,
                          low = beta_I_greta_full_low[,z],
                          high = beta_I_greta_full_high[,z], plot_main = z, cex = 2, lwd = 2, lwd.means = 4)
  dev.off()
}

# Interaction effects including the dose-dependent amplification of genotype effect #

alternative_I <- beta_I_greta_full
alternative_I_low <- beta_I_greta_full_low
alternative_I_high <- beta_I_greta_full_high
alternative_I_var <- beta_I_greta_full_var
for (z in colnames(W)) {
  back.g <- colnames(G1)[which(G1[rownames(W)[which(W[,z]>0)][1],]>0)]
  back.m <- colnames(Mall)[which(Mall[rownames(W)[W[,z]>0 & doses>0][1],]>0)]
  prof_0 <- t((beta_GH_greta_full[,back.g,drop = F])) + t((beta_M_greta_full[,back.m]))
  prof_1 <- t((beta_GH_greta_full[,back.g,drop = F]))*(1+alpha_GM_greta_var[match(z,colnames(W))]) + 
    t((beta_M_greta_full[,back.m])) * t(exp(beta_I_greta_full[,z]))
  alternative_I[,z] <- log(prof_1 / prof_0)
  
  prof_var <- t((beta_GH_greta_full)[,back.g,drop = F])**2 * (alpha_GM_greta_var[match(z,colnames(W))]) + 
    t((beta_M_greta_full[,back.m]))**2 * t(exp(beta_I_greta_full[,z])) * t(beta_I_greta_full_var[,z])
  
  alternative_I_var[,grep(z,colnames(beta_I_greta_full))] <- prof_var / prof_1**2
  
  alternative_I_low[,z] <- alternative_I[,z] - 1.96 * sqrt(alternative_I_var[,z]) 
  alternative_I_high[,z] <- alternative_I[,z] + 1.96 * sqrt(alternative_I_var[,z]) 

}  

# Inidvidual examples
for (z in c('agt.1.MMS','polk.1.MMS','agt.1.EMS','polk.1.EMS','xpc.1.UV','xpf.1.AristolochicAcid','rev.3.UV','polh.1.EMS')) {
  pdf(paste0(z,'.alternative.I.pdf'),12,3)
  par(mar = c(2, 4, 4, 2) + 0.1)
  current_at <- c(-1,0,1,2)
  #current_labels <- c('<0.1',1,10,100)
  if (z %in% c('rev.3.UV')) {
    current_at <-  c(-1,0,1)
    current_labels <- c('<0.1',1,10)
  }
  interaction_effect_plot(alternative_I[,z], CI = T, low = alternative_I_low[,z], high = alternative_I_high[,z],
                          plot_main = z, at = current_at)
  dev.off()
}


# Visualizing actual signatures

unit <- c(1, 10, 10, 100, 0.1, 100, 1, 10, 1, 100, 10, 1)
#  1 microM, 10 microM, 10 microM, 100 Gray, 100 microM, 100 Gray, 1 mM, 10 milliM, 1 mM, 100 Joule, 10 microM, 1 mM 
names(unit) = names(avdose) <- colnames(Mall)
for (z in c('agt.1.MMS','polk.1.MMS','agt.1.EMS','polk.1.EMS','xpc.1.UV','xpf.1.AristolochicAcid','rev.3.UV','polh.1.EMS')) {
  zg = paste(unlist(strsplit(z,split='[.]'))[1:2],collapse='.')
  if (substr(zg,1,3) == 'bub') {
    zg <- substr(z,1,13)
  }
  if (substr(zg,1,5) == 'rad.T') {
    zg <- substr(z,1,15)
  }
  zm = unlist(strsplit(z,split='[.]'))[length(unlist(strsplit(z,split='[.]')))]
  if (zm == 'B1') zm <- 'Aflatoxin.B1'
  tmp <- data.frame(v1 = beta_M_greta_full[,zm] / avdose[zm] * unit[zm],
                    v2 = beta_GH_greta_full[,zg,drop = F]*(alpha_GM_greta[match(z,colnames(W))] / avdose[zm] * unit[zm]) + 
                      beta_M_greta_full[,zm] * exp(beta_I_greta_full[,z]) * unit[zm] / avdose[zm])
  colnames(tmp) <- c(zm,z)
  tmp_var <- data.frame(v1 = beta_M_greta_full_var[,zm] * (unit[zm] / avdose[zm])**2,
                        v2 = beta_GH_greta_full[,zg,drop = F]**2 * (unit[zm] / avdose[zm])**2 * (alpha_GM_greta_var[match(z,colnames(W))]) + 
                          (unit[zm] / avdose[zm])**2 * beta_M_greta_full[,zm]**2 * exp(beta_I_greta_full[,z]) * beta_I_greta_full_var[,z])
  tmp_low <- tmp - 1.96*sqrt(tmp_var)
  tmp_high <- tmp + 1.96*sqrt(tmp_var)
  tmp_low[tmp_low<0] <- 0
  q2 <- plot_dnvhugesig_wb(tmp, CI = T,
                           low = tmp_low,
                           high = tmp_high,
                           ymax = max(tmp_high),
                           norm = F,
                           colors = c("#2EBAED","#000000","#DE1C14","#D4D2D2","#ADCC54","#F0D0CE","brown","#8DD3C7","goldenrod","#BEBADA","darkmagenta"),
                           diff_scale = T, diff_limits = c(20,20)) +
    theme(panel.grid = element_blank())
  q2
#                           diff_scale = T, diff_limits = c(round(max(tmp_high[,1])+1),round(max(tmp_high[,2])+1))) +
#    theme(panel.grid = element_blank())
  print(q2)
}


############################################################################################################################

# Contributions from different factor groups

mu <- ((as.matrix(G1) %*% t((beta_GH_greta_full))) * ((g + as.matrix(W2) %*% t(t(alpha_G_greta)) + 
                                                 doses * (as.matrix(W) %*% t(t(alpha_GM_greta)))) %*% matrix(1,nrow = 1,ncol=119)) + 
  (as.matrix(Mall) %*% t((beta_M_greta_full))) * exp(as.matrix(W) %*% t(beta_I_greta_full)))

genetic_cont <- ((as.matrix(G1) %*% t((beta_GH_greta_full))) * ((g + as.matrix(W2) %*% t(t(alpha_G_greta))) %*% matrix(1,nrow = 1,ncol=119)))
genmut_cont <- ((as.matrix(G1) %*% t((beta_GH_greta_full))) * ((doses * (as.matrix(W) %*% t(t(alpha_GM_greta)))) %*% matrix(1,nrow = 1,ncol=119)))
matrix_of_dif <- (as.matrix(Mall) %*% t((beta_M_greta_full)) * exp(as.matrix(W) %*% t(beta_I_greta_full)) - as.matrix(Mall) %*% t((beta_M_greta_full)) + genmut_cont)
pdf('factor_group_contributions.pdf',width = 8,height = 6)
par(mar = c(8,5,5,5))
f <- barplot(c(sum(rowSums(genetic_cont)[rowSums(W)>0]), 
               sum(rowSums(as.matrix(Mall) %*% t((beta_M_greta_full)))[rowSums(W)>0]), 
               sum(matrix_of_dif[matrix_of_dif>0]),
               sum(matrix_of_dif[matrix_of_dif<0])),
             main = 'Mutations attributed to different factors',
             col = c('palegreen3', 'skyblue', 'lightpink', 'darksalmon'),
             names.arg = c('Endogenous\n processes', 'Exogenous\n mutagens', 'Interactions+','Interactions-'),
             ylim = c(-20000,100000),
             las = 2, yaxt = 'n', border = NA)

axis(side = 2, at = c(-10000,0,20000,50000,100000), labels = c(-10,0,20,50,100), las = 2)
text(x = as.vector(f), y = c(sum(rowSums(genetic_cont)[rowSums(W)>0]), 
                             sum(rowSums(as.matrix(Mall) %*% t((beta_M_greta_full)))[rowSums(W)>0]), 
                             sum(matrix_of_dif[matrix_of_dif>0]),0) + 1000,
     labels = c(paste0(100*round(sum(rowSums(genetic_cont)[rowSums(W)>0]) / 
                                   sum(rowSums(mu)[rowSums(W)>0]),2), ' %'),
                paste0(100*round(sum(rowSums(as.matrix(Mall) %*% t((beta_M_greta_full)))[rowSums(W)>0]) / sum(rowSums(mu)[rowSums(W)>0]),2), ' %'),
                paste0(100*round(sum(matrix_of_dif[matrix_of_dif>0]) / sum(rowSums(mu)[rowSums(W)>0]),2), ' %'),
                paste0(100*round(sum(matrix_of_dif[matrix_of_dif<0]) / sum(rowSums(mu)[rowSums(W)>0]),2), ' %')),
     pos = 3, col = "black")
dev.off()

############################################################################################################################

# Relative change in the total number of base subs


# Calculate distributions for changes in total numbers of muts, and in subs and in the rest separately
changes_all <- list()
changes_subs <- list()
similarities <- list()
k <- 1
p <- 54
Mall <- X[,(p+2):(p+r+1)]
avdose <- NULL
for (j in 1:ncol(Mall)) {
  avdose <- c(avdose, mean(Mall[Mall[,j]>0,j]))
  Mall[,j] = Mall[,j] / mean(Mall[Mall[,j]>0,j])
}
unit <- c(1, 10, 10, 100, 0.1, 100, 1, 10, 1, 100, 10, 1)
#  1 microM, 10 microM, 10 microM, 100 Gray, 100 microM, 100 Gray, 1 mM, 10 milliM, 1 mM, 100 Joule, 10 microM, 1 mM
unit_name <- c('muM','muM','muM','Gy', 'muM', 'Gy', 'mM', 'mM', 'mM','J/m2', 'muM','mM')
names(unit_name) = names(unit) = names(avdose) <- colnames(Mall)
cosine <- function(x,y) {
  return(sum(x * y) / sqrt(sum(x**2)) / sqrt(sum(y**2)))
}
for (j in 1:2000) {
  beta_GH_greta_full_tmp <- matrix(as.matrix(draws[j,grep('beta_G',colnames(draws))]), nrow = m, ncol = p)
  beta_M_greta_full_tmp <- matrix(as.matrix(draws[j,grep('beta_M',colnames(draws))]), nrow = m, ncol = r)
  alpha_GM_greta_tmp <- as.matrix(draws[j,grep('alpha_G_M',colnames(draws))])
  beta_I_greta_full_tmp <- matrix(as.matrix(draws[j,grep('beta_I',colnames(draws))]), nrow = m, ncol = s)
  changes_all[[k]] <- rep(NA,ncol(W))
  names(changes_all[[k]]) <- colnames(W)
  changes_subs[[k]] <- rep(NA,ncol(W))
  names(changes_subs[[k]]) <- colnames(W)
  similarities[[k]] <- rep(NA,ncol(W))
  names(similarities[[k]]) <- colnames(W)
  for (z in colnames(W)) {
    zg = paste(unlist(strsplit(z,split='[.]'))[1:2],collapse='.')
    if (substr(zg,1,3) == 'bub') {
      zg <- substr(z,1,13)
    }
    if (substr(zg,1,5) == 'rad.T') {
      zg <- substr(z,1,15)
    }
    zm = unlist(strsplit(z,split='[.]'))[length(unlist(strsplit(z,split='[.]')))]
    if (zm == 'B1') zm <- 'Aflatoxin.B1'
    mu_0 <- beta_GH_greta_full_tmp[,match(zg,colnames(beta_GH_greta_full))] + 
      beta_M_greta_full_tmp[,match(zm,colnames(beta_M_greta_full))] / avdose[zm] * unit[zm]
    mu_1 <- beta_GH_greta_full_tmp[,match(zg,colnames(beta_GH_greta_full))] * 
      (1 + alpha_GM_greta_tmp[match(z,colnames(W)),] / avdose[zm] * unit[zm]) +
      beta_M_greta_full_tmp[,match(zm,colnames(beta_M_greta_full))] * 
      exp(beta_I_greta_full_tmp[,match(z,colnames(beta_I_greta_full))]) / avdose[zm] * unit[zm]
    changes_subs[[k]][z] <- (sum(mu_1[1:96]) / sum(mu_0[1:96]))
    changes_all[[k]][z] <- (sum(mu_1) / sum(mu_0))
    similarities[[k]][z] <- cosine(mu_0, mu_1)
  }
  k <- k + 1
  print(k)
}

changes_subs <- do.call('cbind',changes_subs)
changes_all <- do.call('cbind',changes_all)

means_all <- apply(changes_all,1,mean)
sds_all <- apply(changes_all,1,function(x) sd(log(x)))

sapply(names(means_all), function(x) {
  return(1 - pchisq(log(means_all[x])**2 / sds_all[x]**2, df = 1))
}) -> pv
pv_all <- p.adjust(pv,method='BH')

low <- apply(changes_subs,1,quantile,0.025)
high <- apply(changes_subs,1,quantile,0.975)
means <- apply(changes_subs,1,mean)
sds <- apply(changes_subs,1,function(x) sd(log(x)))

o <- order(means)
means <- means[o]
sds <- sds[o]
low <- low[o]
high <- high[o]
X_axis <- c(1:length(means))

sapply(names(means), function(x) {
  return(1 - pchisq(log(means[x])**2 / sds[x]**2, df = 1))
}) -> pv
pv_subs <- p.adjust(pv,method='BH')

pdf('base_subs_fold_changes.pdf', width = 6, height = 8)
par(mar = c(10,4,4,4))
plot(x = X_axis, xaxt = 'n', y = rep(NA,length(X_axis)), ylim = c(log10(0.5),log10(40)), 
     yaxt = 'n', xlab = '', ylab = '', bty='n', main = 'Foldchange in total number of mutations')
for (j in 1:length(X_axis)) {
  cur.col <- 'gray88'; cur.lwd = 1
  if (pv_subs[j] < 0.1) {
    cur.col <- 'gray48'
    cur.lwd <- 2
  }
  polygon(x = c(X_axis[j] - 0.5,rep(X_axis[j] + 0.5,2),rep(X_axis[j] - 0.5,2)),
          y = c(rep(log10(high[j]),2),rep(log10(low[j]),2),log10(high[j])),
          col = cur.col,
          border = 0.01)
  lines(x = c(X_axis[j] - 0.5, X_axis[j] + 0.5), y = rep(log10(means[j]),2), lwd = cur.lwd)
}
abline(h = log10(1), lty = 2)

axis(side = 2, labels = c(0.5,1,2,5,10,20,30,40), at = log10(c(0.5,1,2,5,10,20,30,40)),las = 2)
axis(side = 2, labels = rep('',6), at = log10(c(0.6,0.7,0.8,0.9,3,4)),las = 2, tck = -0.01)
axis(side = 1, 
     labels = names(means),
     at = X_axis,
     las=2, cex.axis = 0.5)
dev.off()

#######################################################################################

# Similarity of profiles
similarities <- do.call('cbind',similarities)
sds <- apply(similarities,1,sd)
low <- apply(1-similarities,1,quantile,0.025)
high <- apply(1-similarities,1,quantile,0.975)
means <- apply(1-similarities,1,mean)

o <- order(means)
low <- low[o]
high <- high[o]
means <- means[o]

X_axis <- c(1:length(means))

pdf('similarities.pdf',width=6,height=6)
par(mar = c(10,4,4,4))
plot(x = X_axis, xaxt = 'n', y = rep(NA,length(X_axis)),
     yaxt = 'n', ylim = c(0,1), xlab = 'Experiment', ylab = 'Distance',
     bty='n', main = 'Distances to the profile without interactions')

for (j in 1:length(X_axis)) {
  cur.col <- 'lightsteelblue1'; cur.lwd = 1
  if (means[j] > 0.2) {
    cur.col <- 'skyblue4'
    cur.lwd <- 2
  }
  polygon(x = c(X_axis[j] - 0.5,rep(X_axis[j] + 0.5,2),rep(X_axis[j] - 0.5,2)),
          y = c(rep(high[j],2),rep(low[j],2),high[j]),
          col = cur.col,
          border = 0.01)
  lines(x = c(X_axis[j] - 0.5, X_axis[j] + 0.5), y = rep(means[j],2), lwd = cur.lwd)
}

abline(h = (0.2), lty = 2)
abline(h = 0, lty = 2)

axis(side = 2, labels = c(0,0.2,0.4,0.6), at = c(0,0.2,0.4,0.6),las = 2)
axis(side = 1, 
     labels = names(means),
     at = X_axis,
     las=2, cex.axis = 0.5)
dev.off()
################################

# Individual interactions


ints <- sort(unique(c(names(means)[pv_subs < 0.1],
                    names(means_s)[means_s > 0.2])))

pdf('SupplementaryFigure5.pdf',30,60)
par(mfrow = c(21,5))
line_X_axis <- cumsum(c(rep(1,101),2,4,3,4,2,rep(3,6),2,4,3,4,2,rep(3,6),2,rep(2.5,7)))
for (w in ints) {
  interaction_effect_plot(alternative_I[,w], at = c(-1,0,1,2), labels = c('<0.1',1,10,100), CI = T,
                          low = alternative_I_low[,w],
                          high = alternative_I_high[,w],
                          plot_main = paste0('Effects for ', w),
                          cex = 1.5, lwd = 1, lwd.means = 2, log = T)
}
dev.off()

