# This is meant to run on Amphiprion after having removed low frequency MAFs and downsampling to 663
#### ABC test ####
library(abc)
# setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Migrate-adults/Arlequin") # this was a test
# setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS")
# params <- read.table("migrationABC_DNA1millionby100.params", header = TRUE)
# dat <- read.table("migrationABC_DNA1millionby100_sub200_jointMAFpop1_0.obs", header = TRUE)

# Set working directory on Amphiprion
setwd("~/22_Amarel")

# Read in the chunks of manipulated SFSs that were saved as RData files
load(file = "~/22_Amarel/sub.sfs.all2.1to14267.RData") # sub.sfs.all2
load(file = "~/22_Amarel/sub.sfs.all2.14260to25000.RData") # sub.sfs.all

load(file = "~/22_Amarel/sub.sfs.all1.1to12500RData") # sub.sfs.all1.1
load(file = "~/22_Amarel/sub.sfs.all1.12501to25000RData") # sub.sfs.all1.2

load(file = "~/22_Amarel/sub.sfs.all3.1to12500.RData") # sub.sfs.all3.1
load(file = "~/22_Amarel/sub.sfs.all3.12501to25000.RData") # sub.sfs.all3.2

load(file = "~/22_Amarel/sub.sfs.all4.1to12500.RData") # sub.sfs.all4.1
load(file = "~/22_Amarel/sub.sfs.all4.12501to25000.RData") # sub.sfs.all4.2

# Manipulate the [,,25001:50000] simulations because some of them were downsampled twice
sub.sfs.all2.1 <- sub.sfs.all2[1:14259,]  # 14259 x 52845
sub.sfs.all2.2 <- sub.sfs.all #10741 x 52845
  
# Join all the maniputed SFS chunks together
sub.sfs.all <- rbind (sub.sfs.all1.1, sub.sfs.all1.2, sub.sfs.all2.1, sub.sfs.all2.2, sub.sfs.all3.1, sub.sfs.all3.2, sub.sfs.all4.1, sub.sfs.all4.2)

# Read in PADE SFS
obs <- read.table("PADEadults_zeros_jointMAFpop1_0.obs", skip = 1, header = TRUE)
obs.sumstats <- matrix(obs, nrow = 1, ncol = 52845) # Converts SFS from data.frame to matrix 1 x 52845 so that it can be compared to the simulated SFSs
# obs.sumstats <- as.vector(as.matrix(obs) # Another way to do the same thing

# Read in params file
setwd("~/22_Amarel/params")
params <- read.table("full_params.txt", header = TRUE)

dim(full_params) # should be 100000 x 3
dim(sub.sfs.all) # should be 100000 x 52845

calcs <- abc(target = obs.sumstats, param = full_params, sumstat = sub.sfs.all, tol = 0.01, method = 'rejection', transf = c('log', 'log', 'logit'), logit.bounds = matrix(c(NA, NA, NA, NA, 0, 1), ncol = 2, byrow = TRUE)) # Get error b/c summary statistics in target doesn't match real data. How to input SFS?

#### Looking at some Amarel simulated SFSs that I happened to downsample twice. How different do they look ####
#### Plot matrix ####
sfs <- matrix(unlist(t(sub.sfs.all[1,])), nrow = 195, ncol = 271)
sfs <- matrix(unlist(t(sub.sfs.all2[14260,])), nrow = 195, ncol = 271)

sfs <- matrix(unlist(t(sub.sfs.all[2,])), nrow = 195, ncol = 271)
sfs <- matrix(unlist(t(sub.sfs.all2[14261,])), nrow = 195, ncol = 271)

sfs <- data.matrix(sfs)
min <- min(sfs, na.rm = TRUE)
max <- max(sfs,na.rm = TRUE)

rbPal <- colorRampPalette(c('blue', 'red'))
ColorLevels <- seq(min, max, length=length(rbPal(max +1)))

# Set layout.  We are going to include a colorbar next to plot.
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1),
       heights=c(1,1))

#plotting margins.  These seem to work well for me.
par(mar = c(5,7,2.5,1), font = 2)

# Plot it up!
image(1:ncol(sfs), 1:nrow(sfs), t(sfs),
      col=rbPal(max+1), xlab="Pop_1", ylab="Pop_2",
      axes=FALSE, zlim = c(min, max),
      main= NA, xlim = c(0.5, 158), ylim = c(0.5, 120))
abline(0,1)

# Now annotate the plot
box()
axis(side = 1, at=seq(1,158,1), labels=0:157,
     cex.axis=1.0)
axis(side = 2, at=seq(1,120,1), labels=0:119, las= 1,
     cex.axis=1)

# Add colorbar to second plot region
par(mar = c(3,2.5,2.5,2))
image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=rbPal(max+1),xlab="",ylab="",xaxt="n", las = 1)

