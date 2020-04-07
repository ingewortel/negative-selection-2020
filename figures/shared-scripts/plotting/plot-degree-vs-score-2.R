argv <- commandArgs( trailingOnly = TRUE )

datafile <- argv[1]
outfile <- argv[2]
type <- argv[3]
maxn <- 1000

library( ggplot2, warn.conflicts = FALSE )
library( dplyr, warn.conflicts = FALSE )
library( data.table )

source("../shared-scripts/plotting/mytheme.R" )

if( type == "lang" ){
	xlab <- "log( letter frequency score + 1 )"
} else {
	xlab <- "log( AA frequency score + 1 )"
}

# read data
data <- fread( datafile, data.table = FALSE )
colnames(data) <- c("peptide","degree","score")
data$logdegree <- log10( data$degree + 1 )
data$logscore <- log10( data$score + 1 )

if( nrow(data) > maxn ){
	data <- sample_n( data, maxn )
}

pearson_r <- round( cor.test( data$logscore, data$logdegree )$estimate, 3)
sink("plots/cortest.txt")
print( cor.test( data$logscore, data$logdegree ) )
sink()

# plotting


ybreaks <- seq(0,3)
ylabels <- sapply( ybreaks, function(x) as.expression( bquote( 10^.(x) ) ) )

xpos <- min(data$logscore)
ypos <- 0.95*max(ybreaks)

p <- ggplot( data, aes( x = logscore, y = logdegree ) ) +
	geom_point( size = 0.2, alpha = 0.3) +
	labs( 	x = xlab,
		y = "exchangeability" ) +
	stat_smooth( method = "loess", se = FALSE, color = "red", size=0.8 ) +
	scale_y_continuous( limits=c(min(ybreaks)-0.1,max(ybreaks)), expand=c(0,0), breaks=ybreaks, labels=ylabels ) +
	annotate("text",x=xpos,y=ypos,label=paste("Pearson's r =",pearson_r),hjust=0,size=2.7) +
	mytheme

ggsave(outfile, width=5.6,height=4,units="cm")
