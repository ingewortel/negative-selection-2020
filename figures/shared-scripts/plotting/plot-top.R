library(ggplot2)
source("../shared-scripts/plotting/mytheme.R")

argv <- commandArgs( trailingOnly = TRUE )

datafile <- argv[1]
rvalue <- as.numeric( argv[2] )
lang <- argv[3]
interval <- argv[4] # "iqr", "sem", or "none"
id <- argv[5]
outfile <- argv[6]


# Read data
d <- read.table( datafile, header=TRUE )

# Subset data for current r and language/pathogen
d <- d[ d$r == rvalue & d$lang == lang , ]


# Adjust lo and hi bounds depending on "interval" setting
if( interval == "sd" ){
	# do nothing. lo_auc and hi_auc already contain mean +/- sd
	d <- d
} else if ( interval == "sem" ){
	# interval is mean +- SEM
	d$lo_top <- d$m_top - d$sem_top
	d$hi_top <- d$m_top + d$sem_top
} else if ( interval == "none" ){
	# remove interval by setting low and up bounds to the mean.
	d$lo_top <- d$m_top
	d$hi_top <- d$m_top
} else {
	print( paste( "WARNING - Unknown interval type:", interval, ", setting interval to 'none' ") )
}


# Adjust bounds to the plotting window
d$lo_top[ d$lo_top < 0 ] <- 0
d$hi_top[ d$hi_top > 100 ] <- 100

ylab <- "% foreign in top 10%"


if( id == "lang" ){
	xlab <- "# training strings"
	pw <- 4.8
	if( lang == "xh" ){ylab <- "% Xhosa in top 10%" }
} else if (id == "peps" ){
	xlab <- "# training peptides (x 1000)"
	d$ntrain <- d$ntrain/1000
	if( lang == "hiv" ){ylab <- "% HIV in top 10%" }
	pw <- 5.3
} else {
	stop( paste0("Unknown sequence type: ",id))
}



# Make plot
p <- ggplot( d, aes( 	x = ntrain,
			y = m_top,
			colour = lang,
			fill = lang,
			group = type,
			linetype= type )) +
	annotate("segment", x = 0, xend = max( d$ntrain ), size=0.3,y = 50, yend = 50 ) +
	geom_line() +
	geom_ribbon( aes( ymin=lo_top, ymax=hi_top ), colour = NA, alpha = 0.15 ) +
	labs( 	x = xlab,
		y = ylab) + 
	scale_y_continuous( limits = c(20,100), expand=c(0,0),breaks=c(50,75,100) ) +
 	coord_cartesian(ylim=c(35,100)) +
	scale_x_continuous( expand=c(0,0) ) +
	guides(colour=FALSE,lty=FALSE, fill=FALSE)+
	scale_colour_manual( values=pathogencolors ) +
	scale_fill_manual( values=pathogencolors ) +
	scale_linetype_manual( values=mylinetypes ) +
	mytheme + 
	theme(
    		legend.position=c(0,1),
   		legend.justification=c(0,1),
    		legend.background=element_rect(fill=NA),
		axis.line.x=element_blank()
  	)
gg_table <- ggplot_gtable(ggplot_build(p))
gg_table$layout$clip[gg_table$layout$name=="panel"] <- "off"
grid.draw(gg_table)

# Save plot
ggsave(outfile,gg_table,width=pw,height=4,unit="cm")


