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

# To get values proportional to eigenvalues
ev <- sfs_pca$sdev^2
ev.varprop <- ev/sum(ev)
ev.varprop[1:25]

# Cumulative variance
cumulative <- cumsum(ev.varprop)

# Plot variance/cumulative variance vs PC
pdf("sfs_pca_summary.pdf")
par(mfrow = c(1,3))
# plot(sfs_pca, type = 'l')
plot(ev[1:2000])
abline(h = 0, col = 'red')
plot(ev.varprop[1:2000])
abline(h = 0, col = 'red')
plot(cumulative[1:2000])
abline(h = 0, col = 'red')
dev.off()

# Plotting PCA
pdf('pca1to3')
par(mfrow = c(2,2))
plot(sfs_pca$x[,1:2])
plot(sfs_pca$x[,2:3])
plot(sfs_pca$x[,c(1,3)])
dev.off()

# Projecting other simulations or observed data into PC space
obs <- read.table("PADEadults_zeros_jointMAFpop1_0.obs", skip = 1, header = TRUE) # read in observed data
# obs.vector <- as.vector(t(obs)) # 'flatten' the observed data by rows, t converts to a matrix
# names(obs.vector) <- colnames(sub) # name the columns so that I can remove columns that were constant in the simulated data
# obs.vector.df <- t(data.matrix(obs.vector[!names(obs.vector) %in% constants])) # need to remove columns that I did for the big SFS matrix; should be 19840, and convert to a matrix or dataframe for predict.prcomp
# names(obs.vector.df) <- colnames(sub_noconstant) # add column names because they somehow disappear when I convert to a dataframe...
# obs_proj <- predict(sfs_pca, newdata = obs.vector.df) # project observed data into PC space

obs.vector <- as.vector(as.matrix(obs)) # 'flatten' observed data by column
names(obs.vector) <- colnames(sub) # name the columns so that I can remove columns that were constant in the simulated data
obs.vector.df <- t(data.matrix(obs.vector[!names(obs.vector) %in% constants])) # need to remove columns that I did for the big SFS matrix; should be 19840, and convert to a matrix or dataframe for predict.prcomp
names(obs.vector.df) <- colnames(sub_noconstant) # add column names because they somehow disappear when I convert to a dataframe...
obs_proj <- predict(sfs_pca, newdata = obs.vector.df) # project observed data into PC space
obs_proj_df <- as.data.frame(obs_proj)

# obs.sumstats.df <- t(data.frame(obs.sumstats)) # needs to be a matrix or dataframe for the predict step
# obs.sumstats.df2 <- t(data.frame(obs.sumstats.df[!colnames(obs.sumstats.df) %in% constants])) # need to remove columns that I did for the big SFS matrix; should be 19840 
# colnames(obs.sumstats.df2) <- colnames(sub_noconstant) # add column names
# obs_proj <- predict(sfs_pca, newdata = obs.sumstats.df2)

sub_rest <- sub.sfs.all[10001:100000,]
sub_rest_proj <- predict(sfs_pca, newdata = sub_rest)

# Plot rest of simulations and observed data onto PC space
pdf("proj_pca.pdf")
par(mfrow = c(2,2))
plot(sfs_pca$x[,c(1,2)])
points(sub_rest_proj[,c(1,2)], col= 'tomato')
points(obs_proj_df[,c(1,2)], col = 'blue')
plot(sfs_pca$x[,c(1,3)])
points(sub_rest_proj[,c(1,3)], col= 'tomato')
points(obs_proj_df[,c(1,3)], col = 'blue')
plot(sfs_pca$x[,c(2,3)])
points(sub_rest_proj[,c(2,3)], col= 'tomato')
points(obs_proj_df[,c(2,3)], col = 'blue')
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

# Loclinear ABC with 25 retained PC's
calcs.pca.loclinear <- abc(target = obs_proj[,1:25], param = params, sumstat = full_sfs_pc[,1:25], tol = 0.005, method = 'loclinear',
                           transf = c('log', 'log', 'logit'), logit.bounds = matrix(c(NA, NA, NA, NA, 0, 0.5), ncol = 2, byrow = TRUE))

# Loclinear ABC with 350 retained PC's
calcs.pca350.loclinear <- abc(target = obs_proj[,1:350], param = params, sumstat = full_sfs_pc[,1:350], tol = 0.005, method = 'loclinear',
                           transf = c('log', 'log', 'logit'), logit.bounds = matrix(c(NA, NA, NA, NA, 0, 0.5), ncol = 2, byrow = TRUE))

# Loclinear ABC with 250 retained PC's
calcs.pca250.loclinear <- abc(target = obs_proj[,1:250], param = params, sumstat = full_sfs_pc[,1:250], tol = 0.005, method = 'loclinear',
                              transf = c('log', 'log', 'logit'), logit.bounds = matrix(c(NA, NA, NA, NA, 0, 0.5), ncol = 2, byrow = TRUE))

# Loclinear ABC with 100 retained PC's
calcs.pca100.loclinear <- abc(target = obs_proj[,1:100], param = params, sumstat = full_sfs_pc[,1:100], tol = 0.005, method = 'loclinear',
                              transf = c('log', 'log', 'logit'), logit.bounds = matrix(c(NA, NA, NA, NA, 0, 0.5), ncol = 2, byrow = TRUE))

# Neuralnet ABC with 100 retained PC's
calcs.pca100.neuralnet <- abc(target = obs_proj[,1:100], param = params, sumstat = full_sfs_pc[,1:100], tol = 0.005, method = 'neuralnet',
                              transf = c('log', 'log', 'logit'), logit.bounds = matrix(c(NA, NA, NA, NA, 0, 0.5), ncol = 2, byrow = TRUE))

# Rejection ABC with 8000 retained PC's
calcs.pca8000.rejection <- abc(target = obs_proj[,1:8000], param = params, sumstat = full_sfs_pc[,1:8000], tol = 0.005, method = 'rejection')

# Rejection ABC with 1000 retained PC's
calcs.pca1000.rejection <- abc(target = obs_proj[,1:1000], param = params, sumstat = full_sfs_pc[,1:1000], tol = 0.005, method = 'rejection')

# Loclinear ABC with 450 retained PC's
calcs.pca450.loclinear <- abc(target = obs_proj[,1:450], param = params, sumstat = full_sfs_pc[,1:450], tol = 0.005, method = 'loclinear',
                           transf = c('log', 'log', 'logit'), logit.bounds = matrix(c(NA, NA, NA, NA, 0, 0.5), ncol = 2, byrow = TRUE))

# Rejection ABC with 450 retained PC's
calcs.pca450.rejection <- abc(target = obs_proj[,1:450], param = params, sumstat = full_sfs_pc[,1:450], tol = 0.005, method = 'rejection')

#### Visualizing rejection ABC posteriors ####
save(calcs.pca, file = "~/22_Amarel/calcs.pca0.005.RData")
save(calcs.pca.loclinear, file = "~/22_Amarel/calcs.pca.loclinear0.005.RData")
save(calcs.pca350.loclinear, file = "~/22_Amarel/calcs.pca350.loclinear0.005.RData")
save(calcs.pca250.loclinear, file = "~/22_Amarel/calcs.pca250.loclinear0.005.RData")
save(calcs.pca100.loclinear, file = "~/22_Amarel/calcs.pca100.loclinear0.005.RData")
save(calcs.pca100.neuralnet, file = "~/22_Amarel/calcs.pca100.neuralnet0.005.RData")
save(calcs.pca8000.rejection, file = "~/22_Amarel/calcs.pca8000.rejection0.005.RData")
save(calcs.pca1000.rejection, file = "~/22_Amarel/calcs.pca1000.rejection0.005.RData")
save(calcs.pca450.loclinear, file = "~/22_Amarel/calcs.pca450.loclinear0.005.RData")
save(calcs.pca450.rejection, file = "~/22_Amarel/calcs.pca450.rejection0.005.RData")

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

# Plot 2D kernel densities
library(MASS)

p1 <- kde2d(log10(calcs.pca$unadj.values[,1]), calcs.pca$unadj.values[,3])
image(p1, xlab = 'log10(POPONE) - unadjusted', ylab = 'DISP - unadjusted')
p2 <- kde2d(log10(calcs.pca$unadj.values[,2]), calcs.pca$unadj.values[,3])
image(p2, xlab = 'log10(POPTWO) - unadjusted', ylab = 'DISP - unadjusted')

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

# Plot 2D kernel densities
library(MASS)

p1 <- kde2d(log10(calcs.pca.loclinear$unadj.values[,1]), calcs.pca.loclinear$unadj.values[,3])
image(p1, xlab = 'log10(POPONE) - unadjusted', ylab = 'DISP - unadjusted')
p2 <- kde2d(log10(calcs.pca.loclinear$unadj.values[,2]), calcs.pca.loclinear$unadj.values[,3])
image(p2, xlab = 'log10(POPTWO) - unadjusted', ylab = 'DISP - unadjusted')
p3 <- kde2d(log10(calcs.pca.loclinear$adj.values[,1]), calcs.pca.loclinear$adj.values[,3])
image(p3, xlab = 'log10(POPONE) - adjusted', ylab = 'DISP - adjusted')
p4 <- kde2d(log10(calcs.pca.loclinear$adj.values[,2]), calcs.pca.loclinear$adj.values[,3])
image(p4, xlab = 'log10(POPTWO) - adjusted', ylab = 'DISP - adjusted')

#### Plotting PC1 vs parameters to examine regression ####
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA")
load(file = "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA/calcs.pca.loclinear0.005.RData") # calcs.pca.loclinear

par(mfrow = c(2,2))
plot(calcs.pca.loclinear$unadj.values[, 'POPONE'] ~ calcs.pca.loclinear$ss[,1])
points(calcs.pca.loclinear$adj.values[, 'POPONE'] ~ calcs.pca.loclinear$ss[,1], col = 'tomato')
plot(calcs.pca.loclinear$unadj.values[, 'POPTWO'] ~ calcs.pca.loclinear$ss[,1])
points(calcs.pca.loclinear$adj.values[, 'POPTWO'] ~ calcs.pca.loclinear$ss[,1], col = 'tomato')
plot(calcs.pca.loclinear$unadj.values[, 'DISP'] ~ calcs.pca.loclinear$ss[,1])
points(calcs.pca.loclinear$adj.values[, 'DISP'] ~ calcs.pca.loclinear$ss[,1], col = 'tomato')

plot(calcs.pca.loclinear$unadj.values[,'POPONE'] ~ calcs.pca.loclinear$adj.values[,'POPONE'])
plot(calcs.pca.loclinear$unadj.values[,'POPTWO'] ~ calcs.pca.loclinear$adj.values[,'POPTWO'])
plot(calcs.pca.loclinear$unadj.values[,'DISP'] ~ calcs.pca.loclinear$adj.values[,'DISP'])

# PCA of the 500 retained simulations
dim(calcs.pca.loclinear$ss) # 500 x 25
retained_pca <- prcomp(calcs.pca.loclinear$ss, center = TRUE, scale. = TRUE)
retained_pca$rotation # loadings
retained_pca$x # scores
retained_pca_df <- as.data.frame(retained_pca$x)

# Project observed first 25 PCs into new PC space
load(file = "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA/obs_proj.RData") # obs_proj
obs_proj <- obs_proj[,1:25]
newdata <- t(matrix(obs_proj))
colnames(newdata) <- names(obs_proj)
obs_proj2 <- predict(retained_pca, newdata = newdata) # project observed SS into PC space
obs_proj2_df <- as.data.frame(obs_proj2)

# PCA of 500 retained simulations using 25 PCs, with oberved projected onto it 
par(mfrow = c(2,2))
plot(retained_pca_df[,1:2])
points(obs_proj2_df[,c(1,2)], col = 'red', pch = 19) #matrices plot really weird
plot(retained_pca$x[,2:3])
points(obs_proj2_df[,c(2,3)], col = 'red', pch = 19)
plot(retained_pca$x[,c(1,3)])
points(obs_proj2_df[,c(1,3)], col = 'red', pch = 19)

#### Loclinear of 350 retained PCs ####
# Looking at priors and posteriors on my computer for ABC done on reduced dimension SFS using loclinear method & 350 PCs/SS
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA")
load(file = "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA/calcs.pca350.loclinear0.005.RData") # calcs.pca350.loclinear

plot(density(log10(calcs.pca350.loclinear$unadj.values[,1])), col='red', xlab = 'log10(POPONE)', main = '', ylim = c(0, 1))
lines(density(log10(calcs.pca350.loclinear$adj.values[,1])))
plot(density(log10(calcs.pca350.loclinear$unadj.values[,2])), col='red', xlab = 'log10(POPTWO)', main = '', ylim = c(0, 1))
lines(density(log10(calcs.pca350.loclinear$adj.values[,2])))
plot(density(calcs.pca350.loclinear$unadj.values[,3]), col='red', xlab = 'DISP', main = '')
lines(density(calcs.pca350.loclinear$adj.values[,3]))

# Plot priors and posteriors
library(KScorrect)

popone <- rlunif(10000, 100, 100000) # exp(10) doesn't seem to make a difference
popone_range <- range(popone)
plot(density(log10(popone), from = log10(popone_range[1]), to = log10(popone_range[2])), ylim = c(0,1), xlab = 'log10(POPONE)', main = '')
lines(density(log10(calcs.pca350.loclinear$unadj.values[,1]), from = log10(popone_range[1]), to = log10(popone_range[2])), col='blue')
lines(density(log10(calcs.pca350.loclinear$adj.values[,1]), from = log10(popone_range[1]), to = log10(popone_range[2])), col='red')

poptwo <- rlunif(10000, 100, 100000)
poptwo_range <- range(poptwo)
plot(density(log10(poptwo), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), ylim = c(0,1), xlab = 'log10(POPTWO)', main = '')
lines(density(log10(calcs.pca350.loclinear$unadj.values[,2]), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), col='blue')
lines(density(log10(calcs.pca350.loclinear$adj.values[,2]), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), col='red')

disp <- runif(10000, 0, 0.5)
disp_range <- range(disp)
plot(density(disp, from = disp_range[1], to = disp_range[2]), ylim = c(0, 5), xlab = 'DISP', main = '')
lines(density(calcs.pca350.loclinear$unadj.values[,3], from = disp_range[1], to = disp_range[2]), col='blue')
lines(density(calcs.pca350.loclinear$adj.values[,3], from = disp_range[1], to = disp_range[2]), col='red')

# Plotting unadjusted values for the retained 500 simulations
par(mfrow = c(2,3))
plot(log10(calcs.pca350.loclinear$unadj.values[,1]) ~ log10(calcs.pca350.loclinear$unadj.values[,2]), xlab = 'log10(POPTWO) - unadjusted', ylab = 'log10(POPONE) - unadjusted')
plot(calcs.pca350.loclinear$unadj.values[,3] ~ log10(calcs.pca350.loclinear$unadj.values[,1]), xlab = 'log10(POPONE) - unadjusted', ylab = 'DISP - unadjusted')
plot(calcs.pca350.loclinear$unadj.values[,3] ~ log10(calcs.pca350.loclinear$unadj.values[,2]), xlab = 'log10(POPTWO) - unadjusted', ylab = 'DISP - unadjusted')
plot(log10(calcs.pca350.loclinear$adj.values[,1]) ~ log10(calcs.pca350.loclinear$adj.values[,2]), xlab = 'log10(POPTWO) - adjusted', ylab = 'log10(POPONE) - adjusted')
plot(calcs.pca350.loclinear$adj.values[,3] ~ log10(calcs.pca350.loclinear$adj.values[,1]), xlab = 'log10(POPONE) - adjusted', ylab = 'DISP - adjusted')
plot(calcs.pca350.loclinear$adj.values[,3] ~ log10(calcs.pca350.loclinear$adj.values[,2]), xlab = 'log10(POPTWO) - adjusted', ylab = 'DISP - adjusted')

#### Loclinear of 250 retained PCs ####
# Looking at priors and posteriors on my computer for ABC done on reduced dimension SFS using loclinear method & 250 PCs/SS
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA")
load(file = "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA/calcs.pca250.loclinear0.005.RData") # calcs.pca250.loclinear

# Plot priors and posteriors
library(KScorrect)

popone <- rlunif(10000, 100, 100000) # exp(10) doesn't seem to make a difference
popone_range <- range(popone)
plot(density(log10(popone), from = log10(popone_range[1]), to = log10(popone_range[2])), ylim = c(0,1), xlab = 'log10(POPONE)', main = '')
lines(density(log10(calcs.pca250.loclinear$unadj.values[,1]), from = log10(popone_range[1]), to = log10(popone_range[2])), col='blue')
lines(density(log10(calcs.pca250.loclinear$adj.values[,1]), from = log10(popone_range[1]), to = log10(popone_range[2])), col='red')

poptwo <- rlunif(10000, 100, 100000)
poptwo_range <- range(poptwo)
plot(density(log10(poptwo), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), ylim = c(0,1), xlab = 'log10(POPTWO)', main = '')
lines(density(log10(calcs.pca250.loclinear$unadj.values[,2]), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), col='blue')
lines(density(log10(calcs.pca250.loclinear$adj.values[,2]), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), col='red')

disp <- runif(10000, 0, 0.5)
disp_range <- range(disp)
plot(density(disp, from = disp_range[1], to = disp_range[2]), ylim = c(0, 5), xlab = 'DISP', main = '')
lines(density(calcs.pca250.loclinear$unadj.values[,3], from = disp_range[1], to = disp_range[2]), col='blue')
lines(density(calcs.pca250.loclinear$adj.values[,3], from = disp_range[1], to = disp_range[2]), col='red')

#### Loclinear of 100 retained PCs ####
# Looking at priors and posteriors on my computer for ABC done on reduced dimension SFS using loclinear method & 100 PCs/SS
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA")
load(file = "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA/calcs.pca100.loclinear0.005.RData") # calcs.pca100.loclinear

# Plot priors and posteriors
library(KScorrect)

popone <- rlunif(10000, 100, 100000) # exp(10) doesn't seem to make a difference
popone_range <- range(popone)
plot(density(log10(popone), from = log10(popone_range[1]), to = log10(popone_range[2])), ylim = c(0,1), xlab = 'log10(POPONE)', main = '')
lines(density(log10(calcs.pca100.loclinear$unadj.values[,1]), from = log10(popone_range[1]), to = log10(popone_range[2])), col='blue')
lines(density(log10(calcs.pca100.loclinear$adj.values[,1]), from = log10(popone_range[1]), to = log10(popone_range[2])), col='red')

poptwo <- rlunif(10000, 100, 100000)
poptwo_range <- range(poptwo)
plot(density(log10(poptwo), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), ylim = c(0,1), xlab = 'log10(POPTWO)', main = '')
lines(density(log10(calcs.pca100.loclinear$unadj.values[,2]), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), col='blue')
lines(density(log10(calcs.pca100.loclinear$adj.values[,2]), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), col='red')

disp <- runif(10000, 0, 0.5)
disp_range <- range(disp)
plot(density(disp, from = disp_range[1], to = disp_range[2]), ylim = c(0, 5), xlab = 'DISP', main = '')
lines(density(calcs.pca100.loclinear$unadj.values[,3], from = disp_range[1], to = disp_range[2]), col='blue')
lines(density(calcs.pca100.loclinear$adj.values[,3], from = disp_range[1], to = disp_range[2]), col='red')

#### Neuralnet of 100 retained PCs ####
# Looking at priors and posteriors on my computer for ABC done on reduced dimension SFS using neuralnet method & 100 PCs/SS
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA")
load(file = "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA/calcs.pca100.neuralnet0.005.RData") # calcs.pca100.neuralnet

# Plot priors and posteriors
library(KScorrect)

popone <- rlunif(10000, 100, 100000) # exp(10) doesn't seem to make a difference
popone_range <- range(popone)
plot(density(log10(popone), from = log10(popone_range[1]), to = log10(popone_range[2])), ylim = c(0,1), xlab = 'log10(POPONE)', main = '')
lines(density(log10(calcs.pca100.neuralnet$unadj.values[,1]), from = log10(popone_range[1]), to = log10(popone_range[2])), col='blue')
lines(density(log10(calcs.pca100.neuralnet$adj.values[,1]), from = log10(popone_range[1]), to = log10(popone_range[2])), col='red')

poptwo <- rlunif(10000, 100, 100000)
poptwo_range <- range(poptwo)
plot(density(log10(poptwo), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), ylim = c(0,1), xlab = 'log10(POPTWO)', main = '')
lines(density(log10(calcs.pca100.neuralnet$unadj.values[,2]), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), col='blue')
lines(density(log10(calcs.pca100.neuralnet$adj.values[,2]), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), col='red')

disp <- runif(10000, 0, 0.5)
disp_range <- range(disp)
plot(density(disp, from = disp_range[1], to = disp_range[2]), ylim = c(0, 5), xlab = 'DISP', main = '')
lines(density(calcs.pca100.neuralnet$unadj.values[,3], from = disp_range[1], to = disp_range[2]), col='blue')
lines(density(calcs.pca100.neuralnet$adj.values[,3], from = disp_range[1], to = disp_range[2]), col='red')

#### Rejection of 8000 retained PCs ####
# Looking at priors and posteriors on my computer for ABC done on reduced dimension SFS using simple rejection method & 8000 PCs/SS
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA")
load(file = "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA/calcs.pca8000.rejection0.005.RData") # calcs.pca8000.rejection

# Plot priors and posteriors
library(KScorrect)

popone <- rlunif(10000, 100, 100000) # exp(10) doesn't seem to make a difference
popone_range <- range(popone)
plot(density(log10(popone), from = log10(popone_range[1]), to = log10(popone_range[2])), ylim = c(0,1), xlab = 'log10(POPONE)', main = '')
lines(density(log10(calcs.pca8000.rejection$unadj.values[,1]), from = log10(popone_range[1]), to = log10(popone_range[2])), col='blue')

poptwo <- rlunif(10000, 100, 100000)
poptwo_range <- range(poptwo)
plot(density(log10(poptwo), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), ylim = c(0,1), xlab = 'log10(POPTWO)', main = '')
lines(density(log10(calcs.pca8000.rejection$unadj.values[,2]), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), col='blue')

disp <- runif(10000, 0, 0.5)
disp_range <- range(disp)
plot(density(disp, from = disp_range[1], to = disp_range[2]), ylim = c(0, 5), xlab = 'DISP', main = '')
lines(density(calcs.pca8000.rejection$unadj.values[,3], from = disp_range[1], to = disp_range[2]), col='blue')

#### Rejection of 1000 retained PCs ####
# Looking at priors and posteriors on my computer for ABC done on reduced dimension SFS using simple rejection method & 1000 PCs/SS
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA")
load(file = "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA/calcs.pca1000.rejection0.005.RData") # calcs.pca1000.rejection

# Plot priors and posteriors
library(KScorrect)

popone <- rlunif(10000, 100, 100000) # exp(10) doesn't seem to make a difference
popone_range <- range(popone)
plot(density(log10(popone), from = log10(popone_range[1]), to = log10(popone_range[2])), ylim = c(0,1), xlab = 'log10(POPONE)', main = '')
lines(density(log10(calcs.pca1000.rejection$unadj.values[,1]), from = log10(popone_range[1]), to = log10(popone_range[2])), col='blue')

poptwo <- rlunif(10000, 100, 100000)
poptwo_range <- range(poptwo)
plot(density(log10(poptwo), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), ylim = c(0,1), xlab = 'log10(POPTWO)', main = '')
lines(density(log10(calcs.pca1000.rejection$unadj.values[,2]), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), col='blue')

disp <- runif(10000, 0, 0.5)
disp_range <- range(disp)
plot(density(disp, from = disp_range[1], to = disp_range[2]), ylim = c(0, 5), xlab = 'DISP', main = '')
lines(density(calcs.pca1000.rejection$unadj.values[,3], from = disp_range[1], to = disp_range[2]), col='blue')

#### Loclinear of 450 retained PCs ####
# Looking at priors and posteriors on my computer for ABC done on reduced dimension SFS using loclinear method & 450 PCs/SS
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA")
load(file = "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA/calcs.pca450.loclinear0.005.RData") # calcs.pca450.loclinear

# Plot priors and posteriors
library(KScorrect)

popone <- rlunif(10000, 100, 100000) # exp(10) doesn't seem to make a difference
popone_range <- range(popone)
plot(density(log10(popone), from = log10(popone_range[1]), to = log10(popone_range[2])), ylim = c(0,1), xlab = 'log10(POPONE)', main = '')
lines(density(log10(calcs.pca450.loclinear$unadj.values[,1]), from = log10(popone_range[1]), to = log10(popone_range[2])), col='blue')
lines(density(log10(calcs.pca450.loclinear$adj.values[,1]), from = log10(popone_range[1]), to = log10(popone_range[2])), col='red')

poptwo <- rlunif(10000, 100, 100000)
poptwo_range <- range(poptwo)
plot(density(log10(poptwo), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), ylim = c(0,1), xlab = 'log10(POPTWO)', main = '')
lines(density(log10(calcs.pca450.loclinear$unadj.values[,2]), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), col='blue')
lines(density(log10(calcs.pca450.loclinear$adj.values[,2]), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), col='red')

disp <- runif(10000, 0, 0.5)
disp_range <- range(disp)
plot(density(disp, from = disp_range[1], to = disp_range[2]), ylim = c(0, 5), xlab = 'DISP', main = '')
lines(density(calcs.pca450.loclinear$unadj.values[,3], from = disp_range[1], to = disp_range[2]), col='blue')
lines(density(calcs.pca450.loclinear$adj.values[,3], from = disp_range[1], to = disp_range[2]), col='red')

#### Rejection of 450 retained PCs ####
# Looking at priors and posteriors on my computer for ABC done on reduced dimension SFS using rejection method & 450 PCs/SS
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA")
load(file = "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA/calcs.pca450.rejection0.005.RData") # calcs.pca450.rejection

# Plot priors and posteriors
library(KScorrect)

popone <- rlunif(10000, 100, 100000) # exp(10) doesn't seem to make a difference
popone_range <- range(popone)
plot(density(log10(popone), from = log10(popone_range[1]), to = log10(popone_range[2])), ylim = c(0,1), xlab = 'log10(POPONE)', main = '')
lines(density(log10(calcs.pca450.rejection$unadj.values[,1]), from = log10(popone_range[1]), to = log10(popone_range[2])), col='blue')

poptwo <- rlunif(10000, 100, 100000)
poptwo_range <- range(poptwo)
plot(density(log10(poptwo), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), ylim = c(0,1), xlab = 'log10(POPTWO)', main = '')
lines(density(log10(calcs.pca450.rejection$unadj.values[,2]), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), col='blue')

disp <- runif(10000, 0, 0.5)
disp_range <- range(disp)
plot(density(disp, from = disp_range[1], to = disp_range[2]), ylim = c(0, 5), xlab = 'DISP', main = '')
lines(density(calcs.pca450.rejection$unadj.values[,3], from = disp_range[1], to = disp_range[2]), col='blue')

#####################################
#### Code to make a nice plot of posterior distributions ####
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA/results/posteriors.png", width=8, height=5, res=300, units="in")

par(
  mar=c(5, 2, 2, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  mfrow =c(1,3),
  oma=c(2,2,0,0),
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA
)

# Rejection ABC
# Looking at priors and posteriors on my computer for ABC done on reduced dimension SFS using rejection method
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA")
load(file = "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/PCA/calcs.pca0.005.RData") # calcs.pca

# Plot priors and posteriors
library(KScorrect)

popone <- rlunif(10000, 100, 100000) # exp(10) doesn't seem to make a difference
popone_range <- range(popone)
plot(density(log10(popone), from = log10(popone_range[1]), to = log10(popone_range[2])), ylim = c(0,0.7), xlab = 'log10(POPONE)', main = '', lty = 2, col = 'gray60')
lines(density(log10(calcs.pca$unadj.values[,1]), from = log10(popone_range[1]), to = log10(popone_range[2])))
mtext('Density', 2, 2.5)
text(x = 2.1, y = 0.65, labels = 'A', font = 2, ps = 18)

poptwo <- rlunif(10000, 100, 100000)
poptwo_range <- range(poptwo)
plot(density(log10(poptwo), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), ylim = c(0,0.7), xlab = 'log10(POPTWO)', ylab = '', main = '', lty = 2, col = 'gray60')
lines(density(log10(calcs.pca$unadj.values[,2]), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])))
text(x = 2.1, y = 0.65, labels = 'B', font = 2, ps = 18)

disp <- runif(10000, 0, 0.5)
disp_range <- range(disp)
plot(density(disp, from = disp_range[1], to = disp_range[2]), ylim = c(0, 4), xlab = 'DISP', ylab = '', main = '', lty = 2, col = 'gray60')
lines(density(calcs.pca$unadj.values[,3], from = disp_range[1], to = disp_range[2]))
text(x = 0.02, y = 3.75, labels = 'C', font = 2, ps = 18)

dev.off()
