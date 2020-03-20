suppressMessages( library( dplyr ) )
library(ggplot2)
source("../shared-scripts/plotting/mytheme.R")

argv <- commandArgs( trailingOnly = TRUE )
datafile <- argv[1]
nsim <- as.numeric( argv[2] )
ntest <- as.numeric( argv[3] )
id <- argv[4]
outfile <- argv[5]

d <- read.table( datafile )
colnames(d) <- c("type","cat","epitope","rep","prec","sim","ntrain","r")


if( id == "lang" ){
	xlab <- "Training set size (# strings)"
	pw <- 5.5
} else if (id == "peps" ){
	xlab <- "Training set size (x1000 peptides)"
	pw <- 5.3 # plot width
	d$ntrain <- d$ntrain/1000
} else {
	stop( paste0("Unknown sequence type: ",id))
}

d <- d %>% filter( cat == "self" )


# Calculate sim number in the 7th column
ntypes <- length(unique(d$type))

# Build lookup list
d0 <- d[d$ntrain==0,]

l <- vector("list",ntypes)
names(l) <- unique(d$type)
for( t in unique(d$type) ){
	
	dc <- d0[d0$type == t,]
	l2 <- vector("list", nsim)
	for( s in seq(1,nsim) ){
		dsim <- dc[dc$sim==s,]
		l2[[s]] <- setNames( dsim$prec, dsim$epitope )
	}

	l[[t]] <- l2

}

# Function to look up precursors at ntrain = 0:
prec0 <- function( l, type, epitope, sim ){
	l[[type]][[sim]][[epitope]]
}

# Apply prec0 to the rows of d to find precursors at n=0 for all epitopes
d$p0 <- apply( d, 1, function(x) prec0( l, x[1], x[3], as.numeric(x[6]) ) )

# Find percentage of surviving precursors for each epitope
d$perc <- 100*d$prec/d$p0

# Mean percentage survival for epitopes within one sim
d2 <- d %>%
	group_by( cat, type, ntrain, sim ) %>%
	summarise( mperc = mean( perc ) ) %>%
	as.data.frame()


sem <- function(x){
	return( sd(x)/sqrt(length(x)) )
}


# Mean over sims & sem
d3 <- d2 %>%
	group_by( cat, ntrain ) %>%
	summarise( m = 100-mean( mperc ), sem = sem(mperc), lo = m - sem( mperc ), hi = m + sem( mperc ) ) %>% #lo = quantile( mperc, 0.25 ) , hi = quantile( mperc, 0.75 ) ) %>%
	as.data.frame()
#d3$m <- 100 - d3$m


if( d$r[1] == 3 ){
	xlims <- c(0,25)
} else {
	xlims <- c(0,260)
}

p <- ggplot( d3, aes( x = ntrain, y = m, colour = cat, fill=cat, linetype=type )) + 
	geom_line() + 
	geom_ribbon( aes( ymin=lo, ymax=hi ), alpha=0.2, colour=NA ) +
	labs( x = xlab, y = "TCR deletion (%)" ) +
	scale_y_continuous( limits=c(0,100), expand=c(0,0)) +
	scale_x_continuous( expand=c(0,0), limits=xlims ) +
	scale_colour_manual( values=pathogencolors ) +
	scale_fill_manual( values=pathogencolors ) +
	scale_linetype_manual( values=mylinetypes ) +
	guides( colour=FALSE, fill=FALSE, linetype=FALSE ) +
	mytheme

ggsave( outfile , width = pw, height = 4, units="cm") 
