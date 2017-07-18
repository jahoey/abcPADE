# This meant to be run on Amphiprion. Should have already moved all the .params files from slurm-out to Amphiprion

# The extra columns should have already also been removed. See cut_cols.sh

# Two ways to do this
setwd("~/22_Amarel/params")
# setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/ABC-adults/SFS/22_Amarel_output/params")

# Way 1
paramsfile_list <- list.files(".", pattern = "new") # should contain 200 file names
tables <- lapply(paramsfile_list, read.table, header = TRUE)
full_params <- do.call(rbind,tables) # data.frame

# Way 2
paramsFiles <- lapply(Sys.glob("new.*.params"), read.table, header = TRUE) # I think this would work if they were in the same directory with different names; but reads in as a list
paramsFiles_full <- do.call(rbind,paramsFiles) # data.frame

write.table(full_params, file = "full_params.txt", sep = "\t", row.names = FALSE)