# Full NB model in greta

library(greta)
library(bayesplot)
source('plotting_functions.R')
library(reshape2)
library(ggplot2)
library(openxlsx)


# Get the annotations
data <- openxlsx::read.xlsx("~/Supplementary Table 1. Sample description for C.elegans experiments.xlsx", sheet = 2)
data$Sample <- as.character(data$Sample)
data$Genotype <- as.character(data$Genotype)
CD2Mutant <- sapply(1:nrow(data), function(j) {
  if (is.na(data$Mutagen[j])) 
    return(paste(data$Genotype[j],data$Generation[j],sep=":"))
  else {
    m <- substr(data$Mutagen[j],1,3)
    if (m == 'Ari') m <- 'AA'
    return(paste(data$Genotype[j],m,data$Drug.concentration[j],sep=":"))
  }
})
names(CD2Mutant) <- data$Sample


# Mutation counts
Y <- openxlsx::read.xlsx("~/Supplementary Table 1. Sample description for C.elegans experiments.xlsx", sheet = 3, rowNames = T)


# get rid of reference samples 
CD2Mutant<- CD2Mutant[rownames(Y)] 
data <- data[match(names(CD2Mutant),data$Sample),]


# dimensions
m=119;p=54;r=12;s=196;n=2721


# Adjust generations for 25%/50%/25% pbty of a heterozygous mutation to be fixed/remain/lost
generation.function <- function(N) {
  if (is.na(N)) return(N)
  if (N==0) return(0)
  if (N==1) return(1)
  alpha = sum(sapply(1:N, function(i) 1/(2^(i-1)))) +
    0.25*sum(sapply(1:(N-1), function(i) sum(sapply(1:i, function(j) 1/(2^(j-1))))))
  return(alpha)
}


# Create the design matrix

interactions <- sapply(which(!is.na(data$Mutagen)), function(j) {
  m <- substr(data$Mutagen[j],1,3)
  if (m == 'Ari') m <- 'AA'
  return(paste(data$Genotype[j],m,sep=":"))
})
names(interactions) <- data$Sample[which(!is.na(data$Mutagen))]

# raw dose
doses <- data$Drug.concentration
doses[is.na(doses)] <- 0 # for MA experiments, replace NA dose with 0

mutagens <- data$Mutagen
mutagens[is.na(mutagens)] <- 'no'
mutagens[mutagens=='MMS/EMS/DMS'] <- 'no' # N2 zero-exposure samples

# only indicators
X <- data.frame(sapply(unique(data$Genotype), function(x) as.numeric(data$Genotype == x)),
                sapply(unique(mutagens)[-1], function(x) as.numeric(mutagens == x)),
                sapply(unique(interactions), function(x) sapply(rownames(Y), function(z) 
                  ifelse(is.na(data$Mutagen[match(z,data$Sample)]), 0, as.numeric(interactions[z] == x))))) # 2721 x 274
rownames(X) <- data$Sample

X.tmp <- data.frame(sapply(unique(data$Genotype), function(x) as.numeric(data$Genotype == x)) * sapply(data$Generation,generation.function),
                sapply(unique(mutagens)[-1], function(x) as.numeric(mutagens == x)) * doses,
                sapply(unique(interactions), function(x) sapply(rownames(Y), function(z) 
                  ifelse(is.na(data$Mutagen[match(z,data$Sample)]), 0, as.numeric(interactions[z] == x))))*doses) # 2721 x 274
#X.tmp$name <- rownames(X.tmp)
#X.tmp <- X.tmp[,c(275,1:274)]
write.csv(X.tmp, file = 'yoda1/X.full.csv', row.names = F)
Y.tmp <- data.frame(rowSums(Y[,1:16]),rowSums(Y[,17:32]),
               rowSums(Y[,33:48]),rowSums(Y[,49:64]),
               rowSums(Y[,65:80]),rowSums(Y[,81:96]),
               rowSums(Y[,97:98]),rowSums(Y[,99:104]),
               rowSums(Y[,105:106]),rowSums(Y[,107:112]),
               rowSums(Y[,113:119]))
colnames(Y.tmp) <- c('C>A','C>G','C>T','T>A','T>C','T>G','MNV','D','DI','I','SV')
#Y.tmp$name <- rownames(Y)
#Y.tmp <- Y.tmp[,c(12,1:11)]
write.csv(Y.tmp, file = 'yoda1/spectrum.csv', row.names = F)


# Generations
g = t(t(sapply(data$Generation,generation.function)))


# Genotypes - indicator matrices
G <-  X[,1:p]
G1 <- X[,1:p]


# Interactions and mutagens in N2 - indicator matrix for interactions with N2:mutagen
W2 <- X[,(p+r+1):ncol(X)]
# Genotype per interaction matrix
l = 208
zero_exposure_samples <- data$Sample[!is.na(data$Mutagen) & data$Drug.concentration==0]
for (j in 1:length(zero_exposure_samples)) {
  if (as.character(data$Mutagen[match(zero_exposure_samples[j],data$Sample)])=="MMS/EMS/DMS")
    next
  colname <- interactions[zero_exposure_samples[j]]
  colname <- gsub(x = colname, replacement = '.', pattern = '[ ]')
  colname <- gsub(x = colname, replacement = '.', pattern = '[,]')
  colname <- gsub(x = colname, replacement = '.', pattern = '[:]')
  colname <- gsub(x = colname, replacement = '.', pattern = '[-]')
  colname <- gsub(x = colname, replacement = '.', pattern = '[)]')
  colname <- gsub(x = colname, replacement = '.', pattern = '[(]')
  W2[match(zero_exposure_samples[j],data$Sample),colname] <- 1
}


# Interactions only
X <- X[,-grep('N2.',colnames(X),fixed=T)]
X <- X[,colSums(X)>0]
W <- X[,(p+r+1):ncol(X)]


# Binary mutagen presence matrix
M1 <- X[,(p+1):(p+r)]
M1[M1>0] <- 1


# Mutagen exposure: use median doses as 1
Mall <- X[,(p+1):(p+r)] * doses
avdose <- NULL
for (j in 1:ncol(Mall)) {
  avdose <- c(avdose, mean(Mall[Mall[,j]>0,j]))
  Mall[,j] = Mall[,j] / mean(Mall[Mall[,j]>0,j])
}

# new doses
doses = t(t(rowSums(Mall)))


# there are 3 batches of EMS samples but the first two are mixed and there is no data to resolve it => add a random dose adjustment for them

EMS_unknown_batch <- as.numeric( mutagens == 'EMS' )
EMS_unknown_batch[grep('CD0796',rownames(Y))[1]:length(EMS_unknown_batch)] <- 0
EMS.samples <- names(CD2Mutant)[EMS_unknown_batch>0]

EMS_per_experiment <- rep(0,sum(EMS_unknown_batch))
j = 1
for (experiment in unique(interactions[EMS.samples])) {
  inds <- which(interactions[EMS.samples] == experiment)
  EMS_per_experiment[inds] <- j
  if (experiment == 'parp-1:EMS') {
    j = j+1
    EMS_per_experiment[inds[c(9:length(inds))]] <- j
  }
  else if (experiment == 'him-6:EMS') {
    j = j+1
    EMS_per_experiment[inds[c(1,11:18)]] <- j
  }
  else if (length(inds)>9) {
    j = j+1
    EMS_per_experiment[inds[10:length(inds)]] <- j
  }
  j <- j+1
}
EMS_unknown_batch[EMS_unknown_batch>0] <- EMS_per_experiment

#############################################################################################################################################
################################################################## MODEL ####################################################################

### Step 1: estimate genetic contributions from MA samples

ma.samples <- data$Sample[is.na(data$Mutagen) & data$Generation>0]
ma.X <- X[ma.samples,] * sapply(data$Generation[match(ma.samples, data$Sample)],generation.function)
ma.X <- ma.X[,colSums(ma.X)>0]
ma.Y <- Y[ma.samples,]

# Specify the model

sigma_GH <- variable(lower = 0)

beta_GH = lognormal(meanlog = 0, sdlog = sigma_GH, dim = c(m, ncol(ma.X))) # lognormal with fitted variance

mu_GH = (ma.X %*% t(beta_GH))

size = 100

prob = size / (size + mu_GH)

distribution(ma.Y) = negative_binomial(prob = prob, size = size * matrix(1,nrow = nrow(ma.Y), ncol = ncol(ma.Y)))

ma.model <- model(beta_GH,sigma_GH)

ma.draws <- mcmc(ma.model, warmup = 2000, n_samples = 5000, thin = 2) # takes a lot of memory and time - better do on HPC

# Diagnostics
mcmc_trace(ma.draws[,grep('beta_GH',colnames(ma.draws[[1]]))[sample(1:(m*ncol(ma.X)),9)]])
mcmc_trace(ma.draws[,grep('sigma_GH',colnames(ma.draws[[1]])),drop=F])

# values

beta_GH_greta <- matrix(colMeans((do.call('rbind',ma.draws)[,grep('beta_GH',colnames(ma.draws[[1]]))])), nrow = m, ncol = ncol(ma.X))

beta_GH_greta_low <- matrix(apply((do.call('rbind',ma.draws)[,grep('beta_GH',colnames(ma.draws[[1]]))]),2,quantile,0.025), nrow = m, ncol = ncol(ma.X))

beta_GH_greta_high <- matrix(apply((do.call('rbind',ma.draws)[,grep('beta_GH',colnames(ma.draws[[1]]))]),2,quantile,0.975), nrow = m, ncol = ncol(ma.X))

beta_GH_greta_var <- matrix(apply((do.call('rbind',ma.draws)[,grep('beta_GH',colnames(ma.draws[[1]]))]),2,var), nrow = m, ncol = ncol(ma.X))

colnames(beta_GH_greta) = colnames(beta_GH_greta_low) = colnames(beta_GH_greta_high) =colnames(beta_GH_greta_var) <- colnames(ma.X)


# Plot some signatures

clrs = c("#2EBAED","#000000","#DE1C14","#D4D2D2","#ADCC54","#F0D0CE","brown","#8DD3C7","goldenrod","#BEBADA","darkmagenta")

plot_dnvhugesig_wb((beta_GH_greta[,1:10]), col = clrs, norm = F, ymax = 1)

# check quality 

pdf('beta_GH_quality.pdf',10,10)

plot(rowSums(as.matrix(ma.Y)), rowSums(as.matrix(ma.X) %*% t(beta_GH_greta)), pch = 16, cex = 0.7)

abline(a = 0, b = 1, col = 'red', lty = 2)

dev.off()

# check convergence

pdf('beta_GH_convergence.pdf',12,8)

mcmc_trace(ma.draws[,grep('beta_GH',colnames(ma.draws[[1]]))[sample(1:(m*ncol(ma.X)),9)]])

mcmc_trace(ma.draws[,grep('beta_GH',colnames(ma.draws[[1]]))[sample(1:(m*ncol(ma.X)),9)]])

mcmc_trace(ma.draws[,grep('beta_GH',colnames(ma.draws[[1]]))[sample(1:(m*ncol(ma.X)),9)]])

mcmc_trace(ma.draws[,grep('sigma_GH',colnames(ma.draws[[1]])),drop=F])

dev.off()


#############################################################################################################################################

### Step 2: mutagen contributions and interactions

# 1. The values of MA experiments will go in as gamma priors; for the cases without MA,  we will assume same prior as for N2.

library(MASS)

beta.draws <- do.call('rbind',ma.draws)[,grep('beta_GH',colnames(ma.draws[[1]]))]
shapes <- matrix(sapply(1:ncol(beta.draws), function(i) { print(i); return(fitdistr(beta.draws[,i],'gamma')[[1]][1]) }), nrow = 119, ncol = 51)
rates <- matrix(sapply(1:ncol(beta.draws), function(i) { print(i); return(fitdistr(beta.draws[,i],'gamma')[[1]][2]) }), nrow = 119, ncol = 51)

new_shapes <- matrix(1,nrow = 119, ncol = ncol(G1))
new_rates <- matrix(1,nrow = 119, ncol = ncol(G1))
for (j in 1:ncol(G1)) {
  if (colnames(G1)[j] %in% colnames(ma.X)) {
    new_shapes[,j] <- shapes[,match(colnames(G1)[j],colnames(ma.X))]
    new_rates[,j] <- shapes[,match(colnames(G1)[j],colnames(ma.X))]
  }
  else {
    new_shapes[,j] <- shapes[,match('N2',colnames(ma.X))]
    new_rates[,j] <- shapes[,match('N2',colnames(ma.X))]
  }
}

# Gamma priors fitted to posterior samples from Step 1
beta_GH <- gamma(shape = new_shapes, rate = new_rates, dim = c(m, ncol(G1)))

# Mutagen signatures
sigma_M <- gamma(shape = 1, rate = 1, dim = c(1,r))
beta_M = lognormal(meanlog = matrix(0,nrow = m,ncol=r), sdlog = matrix(1,nrow=m,ncol=1) %*% sigma_M, dim = c(m, ncol(Mall)))

# Adjustment of generations (account for the fact that mutagen exposure samples are not necessarily 
# 1st generation after knockout / after splitting them and their mates with different expousre dose)
sigma_G = 0.5
alpha_G = normal(mean = 0, sd=sigma_G, truncation = c(0,Inf), dim = c(l,1))
r_G = W2 %*% alpha_G

# Dose adjusment for EMS mixed batches
W_EMS <- matrix(0,nrow = n, ncol = length(unique(EMS_unknown_batch))-1)
for (j in 1:ncol(W_EMS))
  W_EMS[which(EMS_unknown_batch == j),j] <- 1
EMS_dose_adjustment = normal(mean = 0, sd = 0.5, dim = c(length(unique(EMS_unknown_batch))-1,1))
r_doses = exp(W_EMS %*% EMS_dose_adjustment)


# genotype aplification by mutagens
sigma_G_M = 0.5
alpha_G_M = normal(mean = 0, sd=sigma_G_M, truncation = c(0,Inf), dim=c(s,1))
r_GM = W %*% alpha_G_M

# Interaction term - regularized
sigma2_I_2 = gamma(shape = 1, rate = 1, dim = c(1,s))
beta_I_2 = laplace(mu = matrix(0,nrow = m,ncol=s), sigma = matrix(1,nrow=m,ncol=1) %*% sigma2_I_2)


# Bring everything together
mu_GH = as.matrix(G1) %*% t(beta_GH) # genetic background
mu_M = M1 %*% t(beta_M) # mutagen contribution
mu1 = exp(W %*% t(beta_I_2)) # multiplicative genotype - mutation interaction
mu =  mu_GH * ((g + r_G + (r_doses*doses) * r_GM) %*% matrix(1,nrow=1,ncol=m)) + mu1 * mu_M * ((doses*r_doses) %*% matrix(1,nrow=1,ncol=m))
size = 100 # account for a bit of overdispersion
prob = size / (size + mu)
distribution(Y) = negative_binomial(size = size, prob = prob) # fit negative binomial distribution to the counts
model_full <- model(beta_GH, beta_M, beta_I_2, alpha_G, alpha_G_M, sigma2_I_2, sigma_M, EMS_dose_adjustment) # get back all the coefficients

# Do the HMC magic
draws.step2 <- mcmc(model_full, warmup = 2000, n_samples = 5000, thin = 2)

# Convergence diagnostics
pdf('MUconvergence.pdf',12,8)
mcmc_trace(draws.step2[,grep('beta_M',colnames(draws.step2[[1]]))[sample(1:(r*m),16)]])
mcmc_trace(draws.step2[,grep('beta_I',colnames(draws.step2[[1]]))[sample(1:(s*m),16)]])
mcmc_trace(draws.step2[,grep('alpha_G',colnames(draws.step2[[1]]))[sample(1:l,16)]])
mcmc_trace(draws.step2[,grep('alpha_G_M',colnames(draws.step2[[1]]))[sample(1:s,16)]])
mcmc_trace(draws.step2[,grep('sigma2_I_2',colnames(draws.step2[[1]]))[sample(1:s,16)],drop=F])
mcmc_trace(draws.step2[,grep('sigma_M',colnames(draws.step2[[1]])),drop=F])
mcmc_trace(draws.step2[,grep('EMS',colnames(draws.step2[[1]]))[sample(1:sum(EMS_unknown_batch),16)]])
dev.off()

# Extract point estimates for parameter
draws_all <- do.call('rbind',draws.step2)
alpha_G_greta <- colMeans(draws_all[,grep('alpha_G',colnames(draws_all))[1:l]])
alpha_G_greta_low <- apply(draws_all[,grep('alpha_G',colnames(draws_all))[1:l]],2,quantile,0.025)
alpha_G_greta_high <- apply(draws_all[,grep('alpha_G',colnames(draws_all))[1:l]],2,quantile,0.975)
alpha_GM_greta <- colMeans(draws_all[,grep('alpha_G_M',colnames(draws_all))])
alpha_GM_greta_low <- apply(draws_all[,grep('alpha_G_M',colnames(draws_all))],2,quantile,0.025)
alpha_GM_greta_high <- apply(draws_all[,grep('alpha_G_M',colnames(draws_all))],2,quantile,0.975)
alpha_GM_greta_var <- apply(draws_all[,grep('alpha_G_M',colnames(draws_all))],2,var)

pdf('Dose_dependent_genotype_amplification_coefficient_greta_090519.pdf',18,8)
par(mar = c(8,4,2,2))
foo <- barplot(t(alpha_GM_greta)[order(alpha_GM_greta)], ylim = c(0,max(alpha_GM_greta_high)),col = 'white')
axis(side = 1, at = foo, labels = colnames(W)[order(alpha_GM_greta)], cex.axis = 0.4, las = 2)
arrows(x0 = foo, y0 = t(alpha_GM_greta_low)[order(alpha_GM_greta)], y1 = t(alpha_GM_greta_high)[order(alpha_GM_greta)], col ='darkgrey', length=0.01, lwd=2, angle=90)
dev.off()

pdf('Additional_generation_per_sample_greta_090519.pdf',18,8)
par(mar = c(8,4,2,2))
foo <- barplot(t(alpha_G_greta)[order(alpha_G_greta)], ylim = c(0,max(alpha_G_greta_high)),col = 'white')
axis(side = 1, at = foo, labels = colnames(W2)[order(alpha_G_greta)], cex.axis = 0.4, las = 2)
arrows(x0 = foo, y0 = t(alpha_G_greta_low)[order(alpha_G_greta)], y1 = t(alpha_G_greta_high)[order(alpha_G_greta)], col ='darkgrey', length=0.01, lwd=2, angle=90)
dev.off()

beta_GH_greta_full <- matrix(colMeans(exp(draws_all[,grep('beta_GH',colnames(draws_all))])), nrow = m, ncol = p)
beta_GH_greta_full_var <- matrix(apply(exp(draws_all[,grep('beta_GH',colnames(draws_all))]),2,var), nrow = m, ncol = p)
beta_GH_greta_full_low <- matrix(apply(exp(draws_all[,grep('beta_GH',colnames(draws_all))]),2,quantile,0.025), nrow = m, ncol = p)
beta_GH_greta_full_high <- matrix(apply(exp(draws_all[,grep('beta_GH',colnames(draws_all))]),2,quantile,0.975), nrow = m, ncol = p)

colnames(beta_GH_greta_full) <- colnames(G1)
colnames(beta_GH_greta_full_low) <- colnames(G1)
colnames(beta_GH_greta_full_high) <- colnames(G1)

clrs = c("#2EBAED","#000000","#DE1C14","#D4D2D2","#ADCC54","#F0D0CE","brown","#8DD3C7","goldenrod","#BEBADA","darkmagenta")
plot_dnvhugesig_wb((beta_GH_greta_full[,1:10]), CI=T, col = clrs, 
                   low = (beta_GH_greta_full_low[,1:10]),
                   high = (beta_GH_greta_full_high[,1:10]), norm = F)

beta_M_greta_full <- matrix(colMeans(draws_all[,grep('beta_M',colnames(draws_all))]), nrow = m, ncol = r)
beta_M_greta_full_var <- matrix(apply(draws_all[,grep('beta_M',colnames(draws_all))],2,var), nrow = m, ncol = r)
beta_M_greta_full_low <- matrix(apply(draws_all[,grep('beta_M',colnames(draws_all))],2,quantile,0.025), nrow = m, ncol = r)
beta_M_greta_full_high <- matrix(apply(draws_all[,grep('beta_M',colnames(draws_all))],2,quantile,0.975), nrow = m, ncol = r)

colnames(beta_M_greta_full) = colnames(beta_M_greta_full_var) = colnames(beta_M_greta_full_high) = colnames(beta_M_greta_full_low) = colnames(Mall)

plot_dnvhugesig_wb((beta_M_greta_full),CI=T, col = clrs, 
                   low = (beta_M_greta_full_low), high = (beta_M_greta_full_high),
                   ymax = 0.3)

beta_I_greta_full <- matrix(colMeans(draws_all[,grep('beta_I',colnames(draws_all))]), nrow = m, ncol = s)
beta_I_greta_full_var <- matrix(apply(draws_all[,grep('beta_I',colnames(draws_all))],2,var), nrow = m, ncol = s)
beta_I_greta_full_low <- matrix(apply(draws_all[,grep('beta_I',colnames(draws_all))],2,quantile,0.025), nrow = m, ncol = s)
beta_I_greta_full_high <- matrix(apply(draws_all[,grep('beta_I',colnames(draws_all))],2,quantile,0.975), nrow = m, ncol = s)

colnames(beta_I_greta_full) = colnames(beta_I_greta_full_var)=colnames(beta_I_greta_full_low) = colnames(beta_I_greta_full_high) = colnames(W)

EMS_dose_adjustment_greta = colMeans(draws_all[,grep('EMS',colnames(draws.step2[[1]]))])
EMS_dose_adjustment_greta_low = apply(draws_all[,grep('EMS',colnames(draws.step2[[1]]))],2,quantile,0.025)
EMS_dose_adjustment_greta_high = apply(draws_all[,grep('EMS',colnames(draws.step2[[1]]))],2,quantile,0.975)

# Check the quality of fit
mu <- (as.matrix(G1) %*% t((beta_GH_greta_full))) * 
  ((g + as.matrix(W2) %*% t(t(alpha_G_greta)) + (doses * exp(W_EMS * EMS_dose_adjustment_greta)) * (as.matrix(W) %*% t(t(alpha_GM_greta)))) 
   %*% matrix(1,nrow = 1,ncol=119)) +
  (as.matrix(M1) %*% t((beta_M_greta_full))) * exp(as.matrix(W) %*% t(beta_I_greta_full)) *
  ((doses * exp(W_EMS * EMS_dose_adjustment_greta)) %*% matrix(1,nrow = 1,ncol=119))

# -logLik
divergence <- function (a,b) {
  return (a * log ( (a+.Machine$double.eps)/(b + .Machine$double.eps)) - a + b)
} # KL distance
sum(divergence(Y,mu))

# plot for per-sample counts
plot(rowSums(Y), rowSums(mu),pch = 16)
abline(a=0,b=1,col='red')



#########################################################################################################

# Save the matrices with coefficients

rownames(beta_GH_greta_full) <- colnames(Y)
rownames(beta_GH_greta_full_low) <- colnames(Y)
rownames(beta_GH_greta_full_high) <- colnames(Y)
rownames(beta_M_greta_full) <- colnames(Y)
rownames(beta_M_greta_full_low) <- colnames(Y)
rownames(beta_M_greta_full_high) <- colnames(Y)
names(alpha_G_greta) <- colnames(W2)
names(alpha_G_greta_low) <- colnames(W2)
names(alpha_G_greta_high) <- colnames(W2)
names(EMS_dose_adjustment_greta) = names(EMS_dose_adjustment_greta_low) = names(EMS_dose_adjustment_greta_high) <- rown
l <- list(beta_GH_greta_full, beta_GH_greta_full_low, beta_GH_greta_full_high,
          beta_M_greta_full, beta_M_greta_full_low, beta_M_greta_full_high,
          data.frame(mean=alpha_G_greta,low=alpha_G_greta_low,high=alpha_G_greta_high),
          data.frame(mean=alpha_G_greta,low=alpha_G_greta_low,high=alpha_G_greta_high))
openxlsx::write.xlsx(l, file="~/Supplementary Table 2. Signatures of DNA repair deficiencies and mutagen exposures.xlsx",
                     colNames = T, rowNames = T,
                     sheetName = c('Genotypes-mean', "Genotypes-low","Genotypes-high",
                                   'Mutagens-mean', "Mutagens-low","Mutagens-high",
                                   'Generation_adjustment', 'EMS_dose_adjustment'))

rownames(beta_I_greta_full) <- colnames(Y)
rownames(beta_I_greta_full_low) <- colnames(Y)
rownames(beta_I_greta_full_high) <- colnames(Y)
names(alpha_GM_greta) <- colnames(W)
names(alpha_GM_greta_low) <- colnames(W)
names(alpha_GM_greta_high) <- colnames(W)
l <- list(beta_I_greta_full, beta_I_greta_full_low, beta_I_greta_full_high,
          data.frame(mean=alpha_GM_greta,low=alpha_GM_greta_low,high=alpha_GM_greta_high))
openxlsx::write.xlsx(l, file="~/Supplementary Table 3. Interaction effects beta_I on the mutagen signatures and alpha_GM on genotype signatures.xlsx",
                     colNames = T, rowNames = T,
                     sheetName = c('Log-Beta-means', "Log-Beta-low","Log-Beta-high",'Alpha'))
