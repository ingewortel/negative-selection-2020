# This script will take as input the files with precursor counts for a specific
# language/t-value/ntrain combination, and generate two outputs.
# the *-calc.txt file contains data from all NSIM simulations, 



# Load packages
suppressMessages( library(dplyr) )
suppressMessages( library(ggplot2) )

# Read command line options
argv <- commandArgs( trailingOnly = TRUE )

r <- as.numeric( argv[1] )			# the threshold t
n <- argv[2]						# the number of training strings/peptides used
l <- argv[3]						# the language/proteome
nsim <- as.numeric( argv[4] )		# the total number of simulations
fileid <- argv[5]					# either "lang" or "peps" depending on type of data.
dir <- argv[6]						# directory where data is stored.

# model type (contiguous or other)
type <- unlist(strsplit( dir, "/" ) ) 
type <- type[length(type)]

# start with empty dataframes for the data
plotdata <- data.frame()
sim1data <- data.frame()

# Find the correct precursor count file for the current "foreign" and read it.
filename <- paste0(dir,"/output-r",r,"-",fileid,"-",l,"-n",n,"-r",r,".txt")
outname <- paste0(dir,"/analysis/analysis-r",r,"-",fileid,"-",l,"-n",n)
d <- read.table( filename )

# add a column with simulation number
# the first simulation is in the first nrow(d)/nsim lines
lines <- nrow(d)/nsim
d$sim <- ceiling(seq(1,nrow(d))/lines)

# Repeat for the file with self test epitopes, then combine the two into a single
# dataframe with data for 50 self and 50 foreign strings.
fileself <- paste0(dir,"/output-r",r,"-",fileid,"-self-n",n,"-r",r,".txt")
dself <- read.table( fileself )
dself$sim <- ceiling(seq(1,nrow(dself))/(nrow(dself)/nsim))
d <- rbind(d,dself)

# Extract relevant columns with TCR counts in total repertoire and for this sequence.
# The total counts are always in column 3, the counts for this specific string always
# in column 4. 
# The resulting dataframe has 4 columns: "cat" (where 0 = "self", 1 = "foreign"),
#	"pc" (the precursor counts for this string), "pf" (precursor frequency for this string,
#	so 10^6 x count divided by repertoire size), "sim" the simulation number.
rcolumn <- 4
d2 <- data.frame( 	cat = ifelse( d$V1 == "self", 0, 1 ),  
	pc = d[,rcolumn],
	pf = 10^6*d[,rcolumn]/d[,3],
	sim = d$sim )

# if d[,3] = 0 (no TCRs in total repertoire), set pf to 0 (to avoid issues with ranking below)
d2$pf[d[,3]==0] <- 0

# Compute median and iqr of both counts and frequencies within each simulation
dtmp <- d2 %>%
	group_by(cat, sim) %>%
	summarise( median_pc = median(pc),
		   median_pf = median(pf),
		   lo_pc = quantile( pc, 0.25, na.rm=FALSE),
		   hi_pc = quantile( pc, 0.75, na.rm=FALSE),
		   lo_pf = quantile( pf, 0.25, na.rm=FALSE),
		   hi_pf = quantile( pf, 0.75, na.rm=FALSE) ) %>%
	as.data.frame()

# The first output is a dataframe with only data for the first simulation; this is used
# to generate plots with recognition of individual strings.
# (distribution of TCR recognition for different sequences within one simulation)
# This file is saved at the bottom of this script.
sim1 <- dtmp[ dtmp$sim == 1, ]
sim1$lang <- l
sim1$r <- r
sim1$ntrain <- as.numeric( n )
sim1$cat <- ifelse( sim1$cat == 0, "self", l ) 
sim1data <- rbind( sim1data, sim1 )
sim1data$type <- type


# The second output is a dataframe with the "discrimination" measure: how many of the
# best-recognized strings are foreign. This requires further computations:

# Loop over the sims and compute self-foreign discrimination across multiple simulations
for( s in unique(d2$sim) ) {
		
	# extract data for this sim
	sim.data <- d2[ d2$sim == s, ]

	# 1 ) Self and foreign precursor counts/frequencies for this sim
	self_pc <- dtmp$median_pc[dtmp$sim==s & dtmp$cat==0]
	foreign_pc <- dtmp$median_pc[ dtmp$sim==s & dtmp$cat == 1 ]
	self_pf <- dtmp$median_pf[dtmp$sim==s & dtmp$cat==0]
	foreign_pf <- dtmp$median_pf[dtmp$sim==s & dtmp$cat==1]

	# 2 ) Percentage foreign in the top 10 %
	# how many epitopes should be in the top 10%?
	ntest <- nrow( sim.data )
	ntop <- round( 0.1*ntest )

	# order epitopes on precursor frequency rank. Same pf will receive the same rank.
	ep.order <- rank( sim.data$pf )
	sorted.cats <- sim.data$cat[ order(ep.order) ]
	sorted.ranks <- ep.order[ order(ep.order) ]
	top_cats <- tail( sorted.cats , ntop )
	
	# 3 ) In principle we now just count the number of foreign strings among the top
	# ranked strings. But in practice, we sometimes get ties; this is why the following
	# code is a bit longer.

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
		high_foreign <- sum( sim.data$cat[high_foreign_indices] != 0 )

		# how many of the ambiguous rank should still be included?
		nkeep <- ntop - high_ranks

		# what fraction of epitopes with this rank is foreign? 
		# use to compute expected number of foreign in the top
		frac_foreign <- sum( sim.data$cat[stretch.indices] != 0 )/length(stretch.indices)
		exp_foreign <- high_foreign + nkeep * frac_foreign

	} else {
		# just count the number of foreign in the top ranks
		exp_foreign <- sum( top_cats != 0 )
	}


	# 4 ) After counting the number of foreign strings among the top-ranked ones,
	# while dealing appropriately with ties, we just convert that to a percentage
	# and add it to the dataframe with all the relevant information.

	# Final percentage foreign among the best-recognized sequences
	top_foreign <- 100*exp_foreign/ntop
				

	# add to dataframe
	dtmp2 <- data.frame( 	lang = l,
				r = r,
				ntrain = as.numeric( n ),
				sim = s,
				self_pc = self_pc,
				foreign_pc = foreign_pc,
				self_pf = self_pf,
				foreign_pf = foreign_pf,
				top = top_foreign )
	plotdata <- rbind( plotdata, dtmp2 )
}	


# Function to compute the standard error of the mean (doesn't exist in R by default)
sem <- function( vec ){

	sd( vec, na.rm=FALSE ) / sqrt( length( vec ) )

}


# We now have computed discrimination numbers for each individual simulation. 
# Here, process them to get summary stats over the replicates of each simulation:
plotdata <- plotdata %>%
		group_by(r, ntrain, lang ) %>%
		summarise( 	m_self_pc = median( self_pc, na.rm=FALSE ),
				lo_self_pc = quantile( self_pc, 0.25, na.rm=FALSE ),
				hi_self_pc = quantile( self_pc, 0.75, na.rm=FALSE ),
				sem_self_pc = sem( self_pc ),
				m_foreign_pc = median( foreign_pc, na.rm=FALSE ),
				lo_foreign_pc = quantile(foreign_pc, 0.25, na.rm=FALSE ),
				hi_foreign_pc = quantile( foreign_pc, 0.75, na.rm=FALSE ),
				sem_foreign_pc = sem( foreign_pc ),
				m_self_pf = median( self_pf, na.rm=FALSE ),
				lo_self_pf = quantile( self_pf, 0.25, na.rm=FALSE ),
				hi_self_pf = quantile( self_pf, 0.75, na.rm=FALSE ),
				sem_self_pf = sem( self_pf ),
				m_foreign_pf = median( foreign_pf, na.rm=FALSE ),
				lo_foreign_pf = quantile(foreign_pf, 0.25, na.rm=FALSE ),
				hi_foreign_pf = quantile( foreign_pf, 0.75, na.rm=FALSE ),
				sem_foreign_pf = sem( foreign_pf ),
				m_top = mean( top ),
				lo_top = m_top - sd(top), #quantile( top, 0.25 ), m_top - sd(top)
				hi_top = m_top + sd(top), #quantile( top, 0.75 ), m_top + sd(top)
				sem_top = sem( top )
			) %>%
		as.data.frame()
plotdata$type <- type


write.table( plotdata, paste0(outname,"-calc.txt"), quote=FALSE, row.names=FALSE, sep = "\t" )
write.table( sim1data, paste0(outname,"-sim1.txt"), quote=FALSE, row.names=FALSE, sep = "\t" )
