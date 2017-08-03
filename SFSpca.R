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
calcs.pca.loclinear <- abc(target = obs_proj, param = params, sumstat = full_sfs_pc, tol = 0.005, method = 'loclinear',
                           transf = c('log', 'log', 'logit'), logit.bounds = matrix(c(NA, NA, NA, NA, 0, 0.5), ncol = 2, byrow = TRUE))
