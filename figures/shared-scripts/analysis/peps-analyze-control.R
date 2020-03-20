suppressMessages( library(dplyr) )
source("../shared-scripts/analysis/concordance-functions2.R")

# INPUT ARGUMENTS : 
argv <- commandArgs( trailingOnly = TRUE )
file <- argv[1] # beginning of the filename
nsim <- as.numeric( argv[2] )
rvalue <- as.numeric( argv[3] )
outfilename <- argv[4]

# LOOP OVER SIMULATIONS : 
sim.data <- data.frame()
for( s in 1:nsim ){
  filename <- paste0(file,s,".txt")
  
  d <- read.table( filename )
  
  # compute concordance:
  conc <- concordances(d, rvalue, "self", verbose=TRUE)$concordance
  
  # Add to dataframe
  sim.data <- rbind( sim.data, data.frame( concordance=conc))
  
}

write.table( sim.data, file=outfilename, quote=FALSE, row.names=FALSE )
