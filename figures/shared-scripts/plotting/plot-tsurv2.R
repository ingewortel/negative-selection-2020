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
colnames(d) <- c("cat","epitope","rep","prec","sim","ntrain","r")


if( id == "lang" ){
	xlab <- "# training strings"
	pw <- 5.5
} else if (id == "peps" ){
	xlab <- "# training peptides (x1000)"
	pw <- 5.3 # plot width
	d$ntrain <- d$ntrain/1000
} else {
	stop( paste0("Unknown sequence type: ",id))
}



#nsim <- 30
#ntest <- 50




# Calculate sim number in the 7th column
ncats <- length(unique(d$cat))
#d$sim <- unlist( lapply( seq(1,nsim), function(x) rep(x, ncats*ntest ) ) )

# Build lookup list
d0 <- d[d$ntrain==0,]

l <- vector("list",ncats)
names(l) <- unique(d$cat)
for( c in unique(d$cat) ){
	
	dc <- d0[d0$cat == c,]
	l2 <- vector("list", nsim)
	for( s in seq(1,nsim) ){
		dsim <- dc[dc$sim==s,]
		l2[[s]] <- setNames( dsim$prec, dsim$epitope )
	}

	l[[c]] <- l2

}

# Function to look up precursors at ntrain = 0:
prec0 <- function( l, cat, epitope, sim ){
	l[[cat]][[sim]][[epitope]]
}

# Apply prec0 to the rows of d to find precursors at n=0 for all epitopes
d$p0 <- apply( d, 1, function(x) prec0( l, x[1], x[2], as.numeric(x[5]) ) )

# Find percentage of surviving precursors for each epitope
d$perc <- 100*d$prec/d$p0

# Mean percentage survival for epitopes within one sim
d2 <- d %>% 
	group_by( cat, ntrain, sim ) %>%
	summarise( mperc = mean( perc ) ) %>%
	as.data.frame()


sem <- function(x){
	return( sd(x)/sqrt(length(x)) )
}


# Mean over sims & sem
d3 <- d2 %>%
	group_by( cat, ntrain ) %>%
	summarise( m = mean( mperc ), sem = sem(mperc), lo = m - sem( mperc ), hi = m + sem( mperc ) ) %>% #lo = quantile( mperc, 0.25 ) , hi = quantile( mperc, 0.75 ) ) %>%
	as.data.frame()


p <- ggplot( d3, aes( x = ntrain, y = m, colour = cat, fill=cat )) + 
	geom_line() + 
	geom_ribbon( aes( ymin=lo, ymax=hi ), alpha=0.2, colour=NA ) +
	labs( x = xlab, y = "TCR survival (%)" ) +
	scale_y_continuous( limits=c(-3,100), expand=c(0,0)) +
	scale_x_continuous( expand=c(0.03,0) ) +
	scale_colour_manual( values=pathogencolors ) +
	scale_fill_manual( values=pathogencolors ) +
	guides( colour=FALSE, fill=FALSE) +
	mytheme

ggsave( outfile , width = pw, height = 4, units="cm") 
