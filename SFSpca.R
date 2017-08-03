# This is meant to be run on Amphiprion. Reduction of data dimensions: Performs a PCA on a giant matrix of 100000 simulations x 52845 SS (195 x 271 SFS)
setwd('~/22_Amarel')

# Load the data
load(file = 'sub.sfs.all.RData')

# Subsetting giant SFS matrix to first 10000 simulations
sub <- sub.sfs.all[1:10000,]

# Need to remove constant columns
constants <- names(sub[, sapply(sub, function(v) var(v, na.rm=TRUE)==0)]) # 33005 columns are constant
sub_noconstant <- sub[,apply(sub, 2, var, na.rm=TRUE) != 0]
dim(sub_noconstant) # 10000 x 19840

# Perform PCA
sfs_pca <- prcomp(sub_noconstant, center = TRUE, scale. = TRUE)
sfs_pca$rotation # loadings
sfs_pca$x # scores, 10000 x 10000

print(sfs_pca)
summary(sfs_pca)

# Plot variance vs PC
pdf('sfs_pca_summary.pdf')
plot(sfs_pca, type = "l")

# To get values proportional to eigenvalues
ev <- sfs_pca$sdev^2
ev.varprop <- ev/sum(ev)
ev.varprop[1:25]

# Plotting PCA
plot(sfs_pca$x)
plot(sfs_pca$x[,2:3])

# Projecting other simulations or observed data into PC space
obs <- read.table("PADEadults_zeros_jointMAFpop1_0.obs", skip = 1, header = TRUE) # read in observed data
obs.sumstats <- as.vector(as.matrix(obs)) # 'flatten' observed data
obs.sumstats.df <- t(data.frame(obs.sumstats)) # needs to be a matrix or dataframe for the predict step
obs.sumstats.df2 <- t(data.frame(obs.sumstats.df[!colnames(obs.sumstats.df) %in% constants])) # need to remove columns that I did for the big SFS matrix; should be 19840 
colnames(obs.sumstats.df2) <- colnames(sub_noconstant) # add column names
obs_proj <- predict(sfs_pca, newdata = obs.sumstats.df2)

sub_rest <- sub.sfs.all[10001:100000,]
sub_rest_proj <- predict(sfs_pca, newdata = sub_rest)

# Plot rest of simulations and observed data onto PC space
pdf("proj_pca.pdf")
plot(sfs_pca$x[,c(1,3)])
points(sub_rest_proj[,c(1,3)], col= 'tomato')
points(obs_proj[,c(1,3)], col = 'blue')
dev.off()

# Save the projected PC's for 90000 simulations and the PC's for the first 10000 simulations in a single file
full_sfs_pc <- rbind(sfs_pca$x, sub_rest_proj) # should be 100000 x 10000
save(full_sfs_pc, file = "~/22_Amarel/full_sfs_pc.RData") # reduced dimensionality of 100000 simulated SFSs

# Save the observed PC's as an R object for use in abc
save(obs_proj, file = "~/22_Amarel/obs_proj.RData") # reduced dimensionality of observed data from 52845 variables to 10000

#### ABC with PCs ####
library(abc)
setwd("~/22_Amarel")
load(file = "~/22_Amarel/full_sfs_pc.RData") # full_sfs_pc
load(file = "~/22_Amarel/obs_proj.RData") # obs_proj

# Read in params file
setwd("~/22_Amarel/params")
params <- read.table("full_params.txt", header = TRUE)

dim(params) # should be 100000 x 3
dim(full_sfs_pc) # should be 100000 x 10000
length(obs_proj) # should be 1 x 10000

# Rejection ABC
calcs.pca <- abc(target = obs_proj, param = params, sumstat = full_sfs_pc, tol = 0.005, method = 'rejection')

# Loclinear ABC
calcs.pca.loclinear <- abc(target = obs_proj[,1:25], param = params, sumstat = full_sfs_pc[,1:25], tol = 0.005, method = 'loclinear',
                           transf = c('log', 'log', 'logit'), logit.bounds = matrix(c(NA, NA, NA, NA, 0, 0.5), ncol = 2, byrow = TRUE))

#### Visualizing rejection ABC posteriors ####
save(calcs.pca, file = "~/22_Amarel/calcs.pca0.005.RData")
save(calcs.pca.loclinear, file = "~/22_Amarel/calcs.pca.loclinear0.005.RData")

# Rejection ABC
# Looking at priors and posteriors on my computer for ABC done on reduced dimension SFS using rejection method
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA")
load(file = "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA/calcs.pca0.005.RData") # calcs.pca

plot(density(log10(calcs.pca$unadj.values[,1])), col='red', xlab = 'log10(POPONE)', main = '')
plot(density(log10(calcs.pca$unadj.values[,2])), col='red', xlab = 'log10(POPTWO)', main = '')
plot(density(calcs.pca$unadj.values[,3]), col='red', xlab = 'DISP', main = '')

# Plot priors and posteriors
library(KScorrect)

pdf('abc_calcs_pca_tol0.005.pdf')
popone <- rlunif(10000, 100, 100000) # exp(10) doesn't seem to make a difference
popone_range <- range(popone)
plot(density(log10(popone), from = log10(popone_range[1]), to = log10(popone_range[2])), ylim = c(0,0.7), xlab = 'log10(POPONE)', main = '')
lines(density(log10(calcs.pca$unadj.values[,1]), from = log10(popone_range[1]), to = log10(popone_range[2])), col='blue')

poptwo <- rlunif(10000, 100, 100000)
poptwo_range <- range(poptwo)
plot(density(log10(poptwo), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), ylim = c(0,0.7), xlab = 'log10(POPTWO)', main = '')
lines(density(log10(calcs.pca$unadj.values[,2]), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), col='blue')

disp <- runif(10000, 0, 0.5)
disp_range <- range(disp)
plot(density(disp, from = disp_range[1], to = disp_range[2]), ylim = c(0, 4), xlab = 'DISP', main = '')
lines(density(calcs.pca$unadj.values[,3], from = disp_range[1], to = disp_range[2]), col='red')
dev.off()

# Correlations between parameters
pdf('POPvDISP_correlatlions_pca.pdf')
plot(log10(calcs.pca$unadj.values[,1]) ~ log10(calcs.pca$unadj.values[,2]), xlab = 'log10(POPTWO)', ylab = 'log10(POPONE)', col = rgb(0,0,0,0.7))
plot(calcs.pca$unadj.values[,3] ~ log10(calcs.pca$unadj.values[,1]), xlab = 'log10(POPONE)', ylab = 'DISP')
lm1 <- lm(calcs.pca$unadj.values[,3] ~ log10(calcs.pca$unadj.values[,1]))
abline(lm1, col = 'tomato')

plot(calcs.pca$unadj.values[,3] ~ log10(calcs.pca$unadj.values[,2]), xlab = 'log10(POPTWO)', ylab = 'DISP')
lm2 <- lm(calcs.pca$unadj.values[,3] ~ log10(calcs.pca$unadj.values[,2]))
abline(lm2, col = 'tomato')
dev.off()

# Loclinear ABC
# Looking at priors and posteriors on my computer for ABC done on reduced dimension SFS using loclinear method
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA")
load(file = "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA/calcs.pca.loclinear0.005.RData") # calcs.pca.loclinear

plot(density(log10(calcs.pca.loclinear$unadj.values[,1])), col='red', xlab = 'log10(POPONE)', main = '', ylim = c(0, 1))
lines(density(log10(calcs.pca.loclinear$adj.values[,1])))
plot(density(log10(calcs.pca.loclinear$unadj.values[,2])), col='red', xlab = 'log10(POPTWO)', main = '', ylim = c(0, 1))
lines(density(log10(calcs.pca.loclinear$adj.values[,2])))
plot(density(calcs.pca.loclinear$unadj.values[,3]), col='red', xlab = 'DISP', main = '')
lines(density(calcs.pca.loclinear$adj.values[,3]))

# Plot priors and posteriors
library(KScorrect)

pdf('abc_calcs_pcaloclinear_tol0.005.pdf')
popone <- rlunif(1000, 100, 100000)
popone_range <- range(popone)
plot(density(log10(popone), from = log10(popone_range[1]), to = log10(popone_range[2])), ylim = c(0,0.7), xlab = 'log10(POPONE)', main = '')
lines(density(log10(calcs.pca.loclinear$unadj.values[,1]), from = log10(popone_range[1]), to = log10(popone_range[2])), col='blue')
lines(density(log10(calcs.pca.loclinear$adj.values[,1]), from = log10(popone_range[1]), to = log10(popone_range[2])), col='red')

poptwo <- rlunif(1000, 100, 100000)
poptwo_range <- range(poptwo)
plot(density(log10(poptwo), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), ylim = c(0,1), xlab = 'log10(POPTWO)', main = '')
lines(density(log10(calcs.pca.loclinear$unadj.values[,2]), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), col='blue')
lines(density(log10(calcs.pca.loclinear$adj.values[,2]), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), col='red')

disp <- runif(1000, 0, 0.5)
disp_range <- range(disp)
plot(density(disp, from = disp_range[1], to = disp_range[2]), ylim = c(0, 4), xlab = 'DISP', main = '')
lines(density(calcs.pca.loclinear$unadj.values[,3], from = disp_range[1], to = disp_range[2]), col='blue')
lines(density(calcs.pca.loclinear$adj.values[,3], from = disp_range[1], to = disp_range[2]), col='red')
dev.off()

# Correlations between parameters
pdf('POPvDISP_correlatlions_pca.pdf')
plot(log10(calcs.pca.loclinear$unadj.values[,1]) ~ log10(calcs.pca.loclinear$unadj.values[,2]), xlab = 'log10(POPTWO)', ylab = 'log10(POPONE)', col = rgb(0,0,0,0.7))
plot(calcs.pca.loclinear$unadj.values[,3] ~ log10(calcs.pca.loclinear$unadj.values[,1]), xlab = 'log10(POPONE)', ylab = 'DISP')
lm1 <- lm(calcs.pca.loclinear$unadj.values[,3] ~ log10(calcs.pca.loclinear$unadj.values[,1]))
abline(lm1, col = 'tomato')

plot(calcs.pca.loclinear$unadj.values[,3] ~ log10(calcs.pca.loclinear$unadj.values[,2]), xlab = 'log10(POPTWO)', ylab = 'DISP')
lm2 <- lm(calcs.pca.loclinear$unadj.values[,3] ~ log10(calcs.pca.loclinear$unadj.values[,2]))
abline(lm2, col = 'tomato')
dev.off()



