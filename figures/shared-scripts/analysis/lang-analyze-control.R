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

  rcolumn <- 4
  d2 <- data.frame(  cat = ifelse( d$V1 == "self", 0, 1 ),  
                     pf = ifelse( d[,3] == 0, 0, 10^6*d[,rcolumn]/d[,3] ) )
  
  # compute % nonself in the top 10% pfs
  ntest <- nrow( d2 )
  ntop <- round( 0.1*ntest )

	# order epitopes on precursor frequency rank. Same pf will receive the same rank.
	ep.order <- rank( d2$pf )
	sorted.cats <- d2$cat[ order(ep.order) ]
	sorted.ranks <- ep.order[ order(ep.order) ]
	top_cats <- tail( sorted.cats , ntop )

	# check if there are ties to be resolved:
	# tie breaking using expected value (see below)
	top_ranks_plus1 <- tail( sorted.ranks, ntop + 1 )

	if( top_ranks_plus1[1] == top_ranks_plus1[2] ){

		# find the ambiguous rank
		amb_rank <- top_ranks_plus1[2]

		# which epitopes have this ambiguous rank?
		stretch.indices <- which( ep.order == amb_rank )

		# how many epitopes have a higher rank and should be included in any case?
		high_ranks <- sum( ep.order > amb_rank )

		# count the number of those that are foreign
		high_foreign_indices <- which( ep.order > amb_rank)
		high_foreign <- sum( d2$cat[high_foreign_indices] != 0 )

		# how many of the ambiguous rank should still be included?
		nkeep <- ntop - high_ranks

		# what fraction of epitopes with this rank is foreign? 
		# use to compute expected number of foreign in the top
		frac_foreign <- sum( d2$cat[stretch.indices] != 0 )/length(stretch.indices)
		exp_foreign <- high_foreign + nkeep * frac_foreign

	} else {
		# just count the number of foreign in the top ranks
		exp_foreign <- sum( top_cats != 0 )
	}

	# Final percentage foreign among the best-recognized sequences
	top_foreign <- 100*exp_foreign/ntop
  
  # compute concordance:
  conc <- concordances(d, rvalue, "self")$concordance
  
  # Add to dataframe
  sim.data <- rbind( sim.data, data.frame( top=top_foreign, concordance=conc))
  
}

write.table( sim.data, file=outfilename, quote=FALSE, row.names=FALSE )

