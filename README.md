# abcPADE
Scripts for manipulating fastsimcoal2 output for ABC

migrationABC_5729_2_0.5dispersal.est & migrationABC_DNA1millionby100.tpl are the input files for generating simulations with fastsimcoal
readparams.R script reads in many parameter files from a single directory and appends them to each other
readingin_simchunks.R script reads in many joint SFS files from a single directory and appends them to each other
fsc_SFSmanipulation.R script also contains code to read in simulated SFS chunks, removes rare alleles (based on MAF), counts number of SNPs in each simulation & downsamples each SFS to match the number of loci in the observed data. Also contains code to plot a joint SFS.
ABCwithSFS.R script were attempts to perform ABC directly on SFSs
SFSpca.R script does PCA on simulated SFS array & some trial ABC and posterior plotting
