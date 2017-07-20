# Before beginning, I need to append 200 files with 500 SFS together in a certain order
# To be done on Amphiprion. This code reads in multiple .obs files that are in the same directory
# setwd("~/21_renameparams_test")
# setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/22_Amarel_output/")
setwd("~/22_Amarel") # This needs to be the directory where the Amarel SFS files are

file_list <- list.files(".", recursive = TRUE, pattern = ".obs$") # should contain 200 file names, but picks up the observed SFS too
file_list <- file_list[-201]

# file_list <- file_list[1:20]
file_list <- split(file_list, ceiling(seq_along(file_list)/2)) # creates a list of lists with 10 SFS in each

#### Read in this file so I can easily get column and row names ####
# setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/DNA/DNA_1millionby100")
setwd("~/")

# sfs <- read.table("migrationABC_DNA1millionby100_1_1_jointMAFpop1_0.obs", skip = 1)
sfs <- read.table("jointSFS_fornames", skip = 1) 

####################################################################################
# Make sure I'm in the right directory
setwd("~/22_Amarel") # This needs to be the directory where the Amarel SFS files are
# Set up dimensions within each clump of simulations
# =============
# = Read File =
# =============

# ---- Data Dimensions Needed for Read ----
# See "Setup" Section Above
header_lines <- 0
bonus_top_lines <- 2 # these would be gibberish comments to skip in each matrix
good_lines <- 195 # these are the data lines
n_matrices <- 500

# ---- Read in each data file [200] 1 Matrix at a Time [500] ----
full_read <- list() # Create empty list to hold each simulation chunk as they're read in

# Loop and Read
for(i in 1:length(file_list)){
  data_read <- list() # Create empty list for each simulation chunk
  for(j in 1:n_matrices){
    
    # Determine Line of File to Start Reading for this Iteration
    lines_read_previously <- (j-1)*(bonus_top_lines+good_lines)
    line_start <- header_lines + lines_read_previously + (bonus_top_lines)
    
    # Read in Desired Portion of File
    # The 'sep' argument should refer to what separates values in the numeric part 
    # See the `write.table(make_nm(), ...` piece of code under 'Add to File' section above
    data_read[[j]] <- scan(file_list[i], skip=line_start, nlines=good_lines, what = as.list(c("character", rep("integer",271))))[-1] # [-1] gets rid of row names
    data_read[[j]] <- lapply(data_read[[j]], as.integer) # if problems try as.numeric
    
    # Format to Matrix
    # The value returned from `scan()` is a vector
    # The number of columns is inferred; can be set explicitly for clarity, if known
    data_read[[j]] <- matrix(unlist(data_read[[j]]), nrow=good_lines)
  }
  full_read <- append(full_read, data_read) # should be 100000 x 52845
}

# save(full_read, file = "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/22_Amarel_output/full_read_subset1.RData")
# load(file = "/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/22_Amarel_output/full_read_subset1.RData")

save(full_read, file = "~/22_Amarel/full_read.RData") # saving and loading will make it faster to read in the second time
load(file = "~/22_Amarel/full_read.RData")

# Format to array
data_array <- simplify2array(full_read) #dim should be 195 x 271 x 100000

# Name the columns and rows
colnames(data_array) <- colnames(sfs)
rownames(data_array) <- rownames(sfs)

# Save the data array so that I can load it quickly
save(data_array, file = "~/22_Amarel/data_array.RData")
load(file = "~/22_Amarel/data_array.RData")

#### Removing rare SNPs with MAF < 0.05 ####
# First create a complete theoretical SFS with the MAF counts across populations
x <- 0:270
y <- 0:194

count.matrix <- matrix(nrow = 195, ncol = 271)
colnames(count.matrix) <- 0:270
rownames(count.matrix) <- 0:194
for (i in 1:ncol(count.matrix)){
  count.matrix[,i] <- x[i]+y
}

# Duplicate the count matirx into an array matching the dimensions of the simulatlions
count.array <- array(rep(count.matrix, n_matrices*length(file_list)), c(195,271,n_matrices*length(file_list)))
rownames(count.array) <- rownames(sfs)
colnames(count.array) <- colnames(sfs)

# Save the count array as an R object
save(count.array, file = "~/22_Amarel/count.array.RData")
load(file = "~/22_Amarel/count.array.RData")

#### Now I have the a data array of observed SFSs and  theoretical SNP counts(195 x 271 x 100000) ####
load(file = "~/22_Amarel/data_array.RData")
load(file = "~/22_Amarel/count.array.RData")

n_matrices <- 25000

# Divide the simulations so that R can deal with their size
data_array <- data_array[,,1:25000]
data_array <- data_array[,,25001:50000]
data_array <- data_array[,,50001:75000]
data_array <- data_array[,,75001:100000]

count.array <- count.array[,,1:25000]

# Replace across population counts with zero based on where there are zeros in the simulated SFS matrix
count.array[which(data_array == 0)] <- 0

# Divide the array by the number of alleles sampled to get frequency
freq.array <- count.array/464

# Then replace any values less than 0.05 with zero
freq.array[freq.array < .043] <- 0

# Need to convert the frequncies back to numbers
count.array2 <- freq.array * 464

# And now I want to convert the matrix back to the counts of how many times a certain combination of minor alleles in each population occurs
data_array[which(count.array2 == 0)] <- 0
data_array[1:20,1:20,1]

#### How many minor alleles in each population? ####
pop1 <- rowSums(data_array[,,1]) #~40-45% loci retained most of the time, 32% was lowest
hist(pop1)
plot(pop1)
lines(pop1)

pop0 <- colSums(data_array[,,1])
hist(pop0)
plot(pop0)
lines(pop0)

sum(pop0) # number of loci remaining after filtering out low frequency ones
sum(pop1)

#### Need to make sure each simulation has 663 loci ####

# Creates a list of how many loci each simulatlion has
loci <- list()

for(i in 1:n_matrices){
  loci[i] <- sum(data_array[,,i])
}

# Sampling from each simulation seems problematic because I 1) do I do it based off the whole matrix? Then I might not get combinations of mafs by chance
# 2) Each count of a maf combination on its own?
# How to sample randomly while wanting the sample to add up to 663 exactly?

#### Code to subsample and reformat the SFS for a single simulation ####
dat1 <- data_array[,,6] # Need to indicate which array
nonzeros <- which(dat1 != 0)
nonzeros.dat1 <- dat1[nonzeros]

sum(dat1)
sum(nonzeros.dat1[which(nonzeros.dat1 == 1)]) # these two lines add up to sum(dat1), 248 +53 = 301 for simulation 1
sum(nonzeros.dat1[which(nonzeros.dat1 > 1)])

adds <- vector()
index <- which(nonzeros.dat1 > 1)

for(j in 1:length(index)){
  adds <- append(adds, (rep(nonzeros[index][j], (nonzeros.dat1[index][j])-1)))
}

full.nonzeros <- append(nonzeros, adds)

loci_663 <- sample(full.nonzeros,200)
index.counts <- as.data.frame(table(loci_663))
colnames(index.counts)[1] <- "array.index"

#### Code to format back into the SFS matrix ####
# Creating a dataframe with the number of indices equal to the number of entries in a 195 x 271 matrix (52845 entries), all starting with a count of zero
array.index <- data.frame(array.index = 1:length(as.list(data_array[,,1])), Freq = rep(0, length(as.list(data_array[,,1]))))

# Merge array.index with the randomly sampled index.counts so that I have a dataframe of how many counts at each index
array.counts <- merge(array.index,index.counts, by = "array.index", all.x = TRUE)
array.counts <- array.counts[,-2] # Get rid of all the zeros
array.counts[is.na(array.counts)] <- 0

# Format the list of counts into a matrix
sub.dat1 <- matrix(array.counts$Freq.y, nrow = 195, ncol = 271)
dim(sub.dat1)
sub.dat1[1:25,1:25]
data_array[1:25,1:25,6] #compare to the subsetted data

# Name the columns and rows
colnames(sub.dat1) <- colnames(sfs)
rownames(sub.dat1) <- rownames(sfs)




#### Trying to put it all together so that I can loop through all the simulations ####
num_loci <- 663 # set number of loci to keep for each simulation
sub.sfs.all <- data.frame()

for(i in 1:n_matrices){
  index.one <- which(data_array[,,i] == 1) # indices of elements exactly equal to 1
  index.greaterthanone <- which(data_array[,,i] > 1) # indices of elements greater than 1
  adds <- vector()
  for(j in 1:length(index.greaterthanone)){
    # adds <- vector()
    # adds <- append(adds, (rep(index.greaterthanone[j], (data_array[,,i][index.greaterthanone][j])-1)))
    adds <- append(adds, (rep(index.greaterthanone[j], (data_array[,,i][index.greaterthanone][j])-1)))
  }
  full.nonzeros <- c(index.one, index.greaterthanone, adds)
  sub.sfs <- data.frame()
  sub.sfs <- as.data.frame(table(sample(full.nonzeros,num_loci)))
  
  array.index <- data.frame(Var1 = 1:length(as.list(data_array[,,i]))) # Creating a dataframe with the number of indices equal to the number of entries in a 195 x 271 matrix (52845 entries)
  
  array.counts <- merge(array.index,sub.sfs, by = "Var1", all.x = TRUE) # Merge array.index with the randomly sampled index.counts so that I have a dataframe of how many counts at each index
  array.counts[is.na(array.counts)] <- 0 # Replace all the NA's with zeros
  
  sub.dat <- matrix(array.counts$Freq, nrow = 1, ncol = 52845) # Format the list of counts into a matrix (nrow = 195, ncol = 271) or a single row (nrow = 1, ncol = 52845)
  
  sub.sfs.all <- rbind(sub.sfs.all, sub.dat)
}

# Save each simulation chunk as a R datatype
sub.sfs.all1 <- sub.sfs.all
save(sub.sfs.all1, file = "~/22_Amarel/sub.sfs.all1.RData")
sub.sfs.all2 <- sub.sfs.all
save(sub.sfs.all2, file = "~/22_Amarel/sub.sfs.all2.RData")
sub.sfs.all3 <- sub.sfs.all
save(sub.sfs.all3, file = "~/22_Amarel/sub.sfs.all3.RData")
sub.sfs.all4 <- sub.sfs.all
save(sub.sfs.all4, file = "~/22_Amarel/sub.sfs.all4.RData")

sub.sfs.all <- rbind(sub.sfs.all1, sub.sfs.all2, sub.sfs.all3, sub.sfs.all4)

# colnames(sub.sfs.all) <- colnames(sfs) # Name the columns and rows if in matrix form
# rownames(sub.sfs.all) <- rep(rownames(sfs), 10) #PROBLEM because rownames are the same

# write.table(sub.sfs.all, file = "migrationABC_DNA1millionby100_sub200_jointMAFpop1_0.obs", sep = "\t", row.names = TRUE, col.names = NA)

#### Plot ####
# Convert dataframe to matrix
sfs <- data.matrix(sub.dat)

#### Plot matrix ####
min <- min(sfs, na.rm = TRUE)
max <- max(sfs,na.rm = TRUE)

rbPal <- colorRampPalette(c('blue', 'red'))
ColorLevels <- seq(min, max, length=length(rbPal(5)))

# Set layout.  We are going to include a colorbar next to plot.
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1),
       heights=c(1,1))

#plotting margins.  These seem to work well for me.
par(mar = c(5,7,2.5,1), font = 2)

# Plot it up!
image(1:ncol(sfs), 1:nrow(sfs), t(sfs),
      col=rbPal(5), xlab="Pop_1", ylab="Pop_2",
      axes=FALSE, zlim = c(min, max),
      main= NA, xlim = c(0.5, 158), ylim = c(0.5, 110))
abline(0,1)

# Now annotate the plot
box()
axis(side = 1, at=seq(1,158,1), labels=0:157,
     cex.axis=1.0)
axis(side = 2, at=seq(1,110,1), labels=0:109, las= 1,
     cex.axis=1)

# Add colorbar to second plot region
par(mar = c(3,2.5,2.5,2))
image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=rbPal(5),xlab="",ylab="",xaxt="n", las = 1)

#### Number of loci sampled ####
pop0 <- colSums(sfs) #need sfs as a dataframe
hist(pop0, breaks = 135)
plot(pop0)
lines(pop0)

pop1 <- rowSums(sfs)
hist(pop0, breaks = 97)
plot(pop1)
lines(pop1)

###########################################################
#### Often times the SNP counts are really hard to see because they're all concentrated in the bottom left corner. Here is the same plotting code to zoom in ####
#### Plot matrix ####
min <- min(sfs, na.rm = TRUE)
max <- max(sfs,na.rm = TRUE)

rbPal <- colorRampPalette(c('blue', 'red'))
ColorLevels <- seq(min, max, length=length(rbPal(50)))

# Set layout.  We are going to include a colorbar next to plot.
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1),
       heights=c(1,1))

#plotting margins.  These seem to work well for me.
par(mar = c(5,7,2.5,1), font = 2)

# Plot it up!
image(1:ncol(sfs), 1:nrow(sfs), t(sfs),
      col=rbPal(50), xlab="Pop_1", ylab="Pop_2",
      axes=FALSE, zlim = c(min, max),
      main= NA, xlim = c(0.5, 25), ylim = c(0.5, 25))

# Now annotate the plot
box()
axis(side = 1, at=seq(1,25,1), labels=0:24,
     cex.axis=1.0)
axis(side = 2, at=seq(1,25,1), labels=0:24, las= 1,
     cex.axis=1)

# Add colorbar to second plot region
par(mar = c(3,2.5,2.5,2))
image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=rbPal(50),xlab="",ylab="",xaxt="n", las = 1)

#### ABC test ####
library(abc)
# setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Migrate-adults/Arlequin")
setwd("~/22_Amarel")
obs <- read.table("PADEadults_zeros_jointMAFpop1_0.obs", skip = 1, header = TRUE)
obs.sumstats <- as.vector(t(obs)) # Converts SFS from data.frame to a vector by rows so that it can be compared to the simulated SFS

# setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS")
# params <- read.table("migrationABC_DNA1millionby100.params", header = TRUE)
# dat <- read.table("migrationABC_DNA1millionby100_sub200_jointMAFpop1_0.obs", header = TRUE)

dim(full_params) # should be 100000 x 3
dim(sub.sfs.all) # should be 100000 x 52845

calcs <- abc(target = obs.sumstats, param = full_params, sumstat = sub.sfs.all, tol = 0.01, method = 'neuralnet', transf = c('log', 'log', 'logit'), logit.bounds = matrix(c(NA, NA, NA, NA, 0, 1), ncol = 2, byrow = TRUE)) # Get error b/c summary statistics in target doesn't match real data. How to input SFS?
