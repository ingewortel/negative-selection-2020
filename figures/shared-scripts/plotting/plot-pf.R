library(ggplot2)
source("../shared-scripts/plotting/mytheme.R")

argv <- commandArgs( trailingOnly = TRUE )

datafile <- argv[1]
rvalue <- as.numeric( argv[2] )
lang <- argv[3]
freq <- argv[4] # 0 for plotting precursor counts, 1 for plotting precursor frequency.
outfile <- argv[5]
maxpf <- as.numeric( argv[6] )
id <- argv[7]

if( is.na(id) ){ id <- "lang" }

# Read data
d <- read.table( datafile, header=TRUE )

# Subset data for current r and language/pathogen
d <- d[ d$r == rvalue & d$lang == lang , ]




# data depending on "what"
if( freq == 0 ){
	# use precursor counts (pc)
	d$y <- d$median_pc
	d$lo <- d$lo_pc
	d$hi <- d$hi_pc
	ylab <- "Precursors (#)"
} else if ( freq == 1 ){
	# use precursor frequency (pf)
	d$y <- d$median_pf
	d$lo <- d$lo_pf
	d$hi <- d$hi_pf
	ylab <- ifelse( id == "lang", "reacting motifs/million","reacting TCRs/million")#expression(paste('Specific TCRs / ',10^{'6'}))
} else {
	print( paste0("WARNING : freq should be either 0 or 1. Current value : ", freq, ". Plotting precursor counts.") )
	# use precursor counts (pc)
	d$y <- d$median_pc
	d$lo <- d$lo_pc
	d$hi <- d$hi_pc
	ylab <- "Precursors (#)"
}

if( maxpf == 0 ){
	maxpf <- max( d$hi[d$hi!=Inf] )
}

# adjust x axis label for peptides/languages separately
if( id == "lang" ){
	xlab <- "# training strings"
	pw <- 4.5
} else if (id == "peps" ){
	xlab <- "# training peptides (x 1000)"
	d$ntrain <- d$ntrain/1000
	pw <- 6
} else {
	stop( paste0("Unknown sequence type: ",id))
}


# Make plot
p <- ggplot( d, aes( 	x = ntrain,
			y = y,
			colour = cat,
			fill = cat )) +
	geom_line() +
	#geom_line( aes( y = lo ), alpha=0.2, lty=2 ) +
	#geom_line( aes( y = hi ), alpha=0.2, lty=2 ) +
	geom_ribbon( aes( ymin=lo, ymax=hi ), colour=NA, alpha = 0.2 ) +
	labs( 	x = xlab,
		y = ylab) + 
	scale_y_continuous( limits=c(0, maxpf), expand=c(0,0) ) +
	scale_x_continuous( expand=c(0,0) ) +
	guides(colour=FALSE,lty=FALSE, fill=FALSE)+
	scale_colour_manual( values=pathogencolors ) +
	scale_fill_manual( values=pathogencolors ) +
	mytheme + 
	theme(
    		legend.position=c(0,1),
   		legend.justification=c(0,1),
    		legend.background=element_rect(fill=NA)
  	)
print(p)
# Save plot
ggsave(outfile,width=pw,height=4.2,unit="cm")


