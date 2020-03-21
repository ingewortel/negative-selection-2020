library( dplyr, warn.conflicts=FALSE )
library( ggplot2, warn.conflicts=FALSE )
library( scales )
source( "../shared-scripts/plotting/mytheme.R" )

argv <- commandArgs( trailingOnly = TRUE )
file <- argv[1]
ntrain <- as.numeric( argv[2] )
traintype <- argv[3]
outplot <- argv[4]

# Read file and preprocess
d <- read.table( file )
colnames(d) <- c( "r","n","type","sim","surv")

# Preprocess
surv0 <- d$surv[d$n==0][1]	# the number of TCRs before selection
d2 <- d %>%
	filter( n >1 ) %>% #filter( n == ntrain ) %>%
	filter( type == traintype ) %>%
	group_by( r, n ) %>% 
	summarise( s = mean(surv) ) %>%
	as.data.frame()

# Percentage TCRs remaining
d2$premain <- 100*d2$s/surv0
d2$pdeleted <- 100-d2$premain

yax_label <- function( x ){

	xout <- x
	for( i in seq_along(x) ){

		if( x[i] >= 1 ){
			xout[i] <- x[i]
		} else {
			pow <- log10(x[i])
			xout[i] <- parse( text = paste0("10^",pow) )
		}

	}
	return( xout )

}

p <- ggplot( d2, aes( x = r, y = d2$pdeleted, color = log2(n), group = n ) ) +
	#geom_point() +
	geom_line() +
	labs( x = "threshold t", y = "TCRs deleted (%)", color = "# training strings" ) +
	scale_x_continuous( limits=c(0.9,6.1), expand=c(0,0), breaks=seq(1,6) ) +
	scale_y_log10( expand=c(0,0), breaks = c(0.000001, 0.0001, 0.01, 1, 100 ), labels = yax_label ) +
	scale_color_gradient( limits = c( log2(min(d2$n)/1.5), log2(max(d2$n)*1.5)), breaks = log2( unique( d2$n ) ), labels = unique( d2$n ) ) +
	coord_cartesian( ylim = c(0.000001,105)) +
	mytheme + theme(
		axis.text.y = element_text( hjust = 0 ),
		legend.position = "right",
		legend.title = element_text()
	)

ggsave( outplot, width = 10, height = 5, units = "cm" )