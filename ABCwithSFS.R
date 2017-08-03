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

load(file = "~/22_Amarel/sub.sfs.all1.1to12500.RData") # sub.sfs.all1.1
load(file = "~/22_Amarel/sub.sfs.all1.12501to25000.RData") # sub.sfs.all1.2

load(file = "~/22_Amarel/sub.sfs.all3.1to12500.RData") # sub.sfs.all3.1
load(file = "~/22_Amarel/sub.sfs.all3.12501to25000.RData") # sub.sfs.all3.2

load(file = "~/22_Amarel/sub.sfs.all4.1to12500.RData") # sub.sfs.all4.1
load(file = "~/22_Amarel/sub.sfs.all4.12501to25000.RData") # sub.sfs.all4.2

# Manipulate the [,,25001:50000] simulations because some of them were downsampled twice & rename the [,,25001:50000] SFS so naming convention matches others
sub.sfs.all2.1 <- sub.sfs.all2[1:14259,]  # 14259 x 52845
sub.sfs.all2.2 <- sub.sfs.all #10741 x 52845
  
# Join all the maniputed SFS chunks together
sub.sfs.all <- rbind (sub.sfs.all1.1, sub.sfs.all1.2, sub.sfs.all2.1, sub.sfs.all2.2, sub.sfs.all3.1, sub.sfs.all3.2, sub.sfs.all4.1, sub.sfs.all4.2)
save(sub.sfs.all, file = "~/22_Amarel/sub.sfs.all.RData") # This file is pretty big

library(abc)
setwd("~/22_Amarel")
load(file = "~/22_Amarel/sub.sfs.all.RData")

# Read in PADE SFS
obs <- read.table("PADEadults_zeros_jointMAFpop1_0.obs", skip = 1, header = TRUE)
# obs.sumstats <- matrix(obs, nrow = 1, ncol = 52845) # Converts SFS from data.frame to matrix 1 x 52845 so that it can be compared to the simulated SFSs
obs.sumstats <- as.vector(as.matrix(obs)) # Another way to do the same thing

# Read in params file
setwd("~/22_Amarel/params")
params <- read.table("full_params.txt", header = TRUE)

dim(params) # should be 100000 x 3
dim(sub.sfs.all) # should be 100000 x 52845
length(obs.sumstats) # should be 52845

calcs <- abc(target = obs.sumstats, param = params, sumstat = sub.sfs.all, tol = 0.005, method = 'rejection')

# Let's save this and move to my computer so that I can look at it
# save(calcs, file = "~/22_Amarel/calcs0.01.RData")
save(calcs, file = "~/22_Amarel/calcs0.005.RData")

#### Looking at priors and posteriors on my computer ####
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/22_Amarel_output")
# load(file = "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/22_Amarel_output/calcs0.01.RData") # calcs - tol = 0.01
load(file = "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/22_Amarel_output/calcs0.005.RData") # calcs - tol = 0.005

plot(density(log10(calcs$unadj.values[,1])), col='red', xlab = 'log10(POPONE)', main = '')
plot(density(log10(calcs$unadj.values[,2])), col='red', xlab = 'log10(POPTWO)', main = '')
plot(density(calcs$unadj.values[,3]), col='red', xlab = 'DISP', main = '')

# Plot priors and posteriors
library(KScorrect)

# pdf('abc_calcs_tol0.01.pdf')
pdf('abc_calcs_tol0.005.pdf')
popone <- rlunif(1000, 100, 100000) # exp(10) doesn't seem to make a difference
popone_range <- range(popone)
plot(density(log10(popone), from = log10(popone_range[1]), to = log10(popone_range[2])), ylim = c(0,0.5), xlab = 'log10(POPONE)', main = '')
lines(density(log10(calcs$unadj.values[,1]), from = log10(popone_range[1]), to = log10(popone_range[2])), col='blue')
# plot(density(popone, bw = 2000, adjust = 1, cut = 3))
# lines(density(calcs$unadj.values[,1], bw = 2000, adjust = 1, cut = 3), col='blue')

poptwo <- rlunif(1000, 100, 100000)
poptwo_range <- range(poptwo)
plot(density(log10(poptwo), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), ylim = c(0,0.5), xlab = 'log10(POPTWO)', main = '')
lines(density(log10(calcs$unadj.values[,2]), from = log10(poptwo_range[1]), to = log10(poptwo_range[2])), col='blue')
# plot(density(poptwo, bw = 2000, adjust = 1, cut = 3))
# lines(density(calcs$unadj.values[,2], bw = 2000, adjust = 1, cut = 3), col='blue')

disp <- runif(1000, 0, 0.5)
disp_range <- range(disp)
plot(density(disp, from = disp_range[1], to = disp_range[2]), ylim = c(0, 4), xlab = 'DISP', main = '')
lines(density(calcs$unadj.values[,3], from = disp_range[1], to = disp_range[2]), col='red')
# plot(density(disp, bw = 0.02, adjust = 1, cut =3), ylim = c(0, 4))
# lines(density(calcs$unadj.values[,3], bw = 0.02, adjust = 1, cut = 3), col='red')
dev.off()

# Correlations between parameters
pdf('POPvDISP_correlatlions.pdf')
plot(log10(calcs$unadj.values[,1]) ~ log10(calcs$unadj.values[,2]), xlab = 'log10(POPTWO)', ylab = 'log10(POPONE)', col = rgb(0,0,0,0.7))
plot(calcs$unadj.values[,3] ~ log10(calcs$unadj.values[,1]), xlab = 'log10(POPONE)', ylab = 'DISP')
lm1 <- lm(calcs$unadj.values[,3] ~ log10(calcs$unadj.values[,1]))
abline(lm1)

plot(calcs$unadj.values[,3] ~ log10(calcs$unadj.values[,2]), xlab = 'log10(POPTWO)', ylab = 'DISP')
lm2 <- lm(calcs$unadj.values[,3] ~ log10(calcs$unadj.values[,2]))
abline(lm2)
dev.off()

#### Looking at some Amarel simulated SFSs that I happened to downsample twice. How different do they look ####
# Calculting some differences between the observed and simulated SFSsfs
load('~/Downloads/sub.sfs.all2.1to14267.RData') # sub.sfs.all2
test <- sub.sfs.all2[1:2,]

diffs <- t(t(test) - obs.sumstats) # Calculates differences between simulated and observed SFS

setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/22_Amarel_output/")
full_params <- read.table('full_params.txt', header = TRUE)
full_params <- full_params[25001:39267,]

calcs <- abc(target = obs.sumstats, param = full_params, sumstat = sub.sfs.all2, tol = 0.001, method = 'rejection') 

plot(density(log10(calcs$unadj.values[,1])), col='red')
plot(density(log10(calcs$unadj.values[,2])), col='red')
plot(density(calcs$unadj.values[,3]), col='red')

#### Plot matrix ####
sfs <- matrix(unlist(t(sub.sfs.all[1,])), nrow = 195, ncol = 271)
sfs <- matrix(unlist(t(sub.sfs.all2[14260,])), nrow = 195, ncol = 271)

sfs <- matrix(unlist(t(sub.sfs.all[2,])), nrow = 195, ncol = 271)
sfs <- matrix(unlist(t(sub.sfs.all2[14261,])), nrow = 195, ncol = 271)

sfs <- matrix(unlist(t(obs.sumstats)), nrow = 195, ncol = 271)

sfs <- matrix(unlist(t(sfs)), nrow = 195, ncol = 271)

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

