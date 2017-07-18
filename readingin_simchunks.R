# Read in 200 joint SFS files into R and append them together in a particular order
# To be done on Amphiprion, Need to think about how to deal with the lines in between the SFS's
setwd("~/22_Amarel") # This needs to be the directory where the Amarel SFS files are

file_list <- list.files(".", recursive = TRUE, pattern = ".obs$") # should contain 200 file names

# SFSFiles <- lapply(Sys.glob("*.obs"), read.table) # I think this would work if they were in the same directory with different names

#### Read in this file so I can easily get column and row names ####
# setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/DNA/DNA_1millionby100")
# 
# sfs <- read.table("migrationABC_DNA1millionby100_1_1_jointMAFpop1_0.obs", skip = 1)

####################################################################################
# Define File Name For Testing
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/")
file_name <- "migrationABC_DNA1millionby100_jointMAFpop1_0.obs2" #10 simulations
file_name2 <- "migrationABC_DNA1millionby100_jointMAFpop1_0.obs" #10 simulations
# file_name3 <- "~/Downloads/test1.jointMAFpop1_0.obs" # 2 simulations
file_list <- c(file_name,file_name2)
####################################################################################
# Set up dimensions within each clump of simulations
# =============
# = Read File =
# =============

# ---- Data Dimensions Needed for Read ----
# See "Setup" Section Above
header_lines <- 0
bonus_top_lines <- 2 # these would be gibberish comments to skip in each matrix
good_lines <- 195 # these are the data lines
n_matrices <- 500 # This should be 500 for amarel

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
  full_read <- append(full_read, data_read)
}

# Format to array
data_array <- simplify2array(full_read) #dim should be 195 x 271 x 100000

# # Name the columns and rows
# colnames(data_array) <- colnames(sfs)
# rownames(data_array) <- rownames(sfs)