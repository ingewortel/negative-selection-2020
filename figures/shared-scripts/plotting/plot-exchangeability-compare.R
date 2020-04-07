argv <- commandArgs( trailingOnly = TRUE )

datafile <- argv[1]
set1 <- argv[2]
set2 <- argv[3]
outfile <- argv[4]

library( ggplot2 )
library( dplyr, warn.conflicts=FALSE )
library( ggbeeswarm )

source("../shared-scripts/plotting/mytheme.R" )


# Read exchangeabilities and create lookup table
etable <- read.table( "data/exch-all.txt" )
elookup <- etable$V3
names(elookup) <- etable$V1
setnames <- c("all", "optimal", "random" )

# Read the sequences from set 1 and set 2
s1 <- as.character( read.table( set1 )$V1 )
s2 <- as.character( read.table( set2 )$V1 )
s2 <- head( s2, 100000 )
s3 <- sample( s1, length(s2) )

# Get exchangeability using lookup
e1 <- elookup[s1]
e2 <- elookup[s2]
e3 <- elookup[s3]

# Generate a dataframe
d1 <- data.frame( set=setnames[1], ex=e1 )
d2 <- data.frame( set=setnames[2], ex=e2 )
d3 <- data.frame( set=setnames[3], ex=e3 )
d <- rbind( d1, d2, d3 )
d$lex <- log10( d$ex+1 )

d$set <- factor( d$set, levels = c("all","random","optimal"))

d2 <- sample_n( d, 500 )

ybreaks <- seq(0,2)
ylabels <- sapply( ybreaks, function(x) as.expression( bquote( 10^.(x) ) ) )

# Plot
linecolors <- setNames( c(NA,"black","black"),setnames)
linetypes <- setNames( c(0,2,1),setnames)
fillcolors <- setNames( c("gray50","gray50","gray50"),setnames )
#p <- ggplot( d, aes( x = lex, y=..density.., linetype = set, fill = set, group = set ) ) + 
p <- ggplot( d, aes( x = set, y = lex, linetype = set, color = set )) +
	#geom_quasirandom( data=d2, show.legend=FALSE, size = 0.5, alpha = 0.5, pch=16) +
	geom_violin( show.legend=FALSE, lty=0, aes(fill = set), alpha = 0.5 ) +
	geom_boxplot(show.legend=FALSE,lty=1,outlier.shape=NA,fill=NA, color="black") +
	#geom_histogram(position="identity", alpha=0.5,  bins=13) + 
	#geom_freqpoly( position="identity", alpha = 0.5, bins = 30 ) +
	#geom_density( bw = 0.06, show.legend = FALSE  ) +
	labs( x = "training peptides", y = "exchangeability" ) +
	#scale_y_continuous( limits=c(0,100)) +
	#scale_x_continuous( expand=c(0,0) ) +
	scale_y_continuous( limits=c(min(ybreaks)-0.1,2.6), expand=c(0,0), breaks=ybreaks, labels=ylabels ) +
	#scale_y_continuous( expand=c(0,0) ) +
	#scale_colour_manual( values = linecolors ) +
	scale_linetype_manual( values = linetypes ) +
	scale_color_manual( values = fillcolors ) +
	scale_fill_manual( values = fillcolors ) +
	mytheme + theme(
		legend.position = "right"
	)
#print(p)

ggsave(outfile, width = 5.3, height = 4, units="cm")