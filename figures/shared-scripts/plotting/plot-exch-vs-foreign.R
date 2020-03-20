argv <- commandArgs( trailingOnly = TRUE )

datafile <- argv[1]
outfile <- argv[2]
nbins <- 10
nsamplepoints <- 250

library( ggplot2 )
library( dplyr, warn.conflicts=FALSE )
library( ggbeeswarm )

source("../shared-scripts/plotting/mytheme.R" )

# Read data
d <- read.table( datafile )

# columns: peptides, class (self=0,foreign=1), 
# nself/nforeign/ntot (number of self,foreign,tot neighbors)
colnames(d) <- c("peptide", "cls", "nself", "nforeign", "ntot" )

# compute exchangeability
d$exch <- log10( d$nself + 1 )

# take a sample
dwith <- d[ d$nforeign > 0 , ]
if( nrow(dwith) > nsamplepoints ){ dwith <- sample_n( dwith, nsamplepoints ) }
dwithout <- d[ d$nforeign == 0, ]
#dwithout <- dwithout[ sample(1:nrow(dwithout), nrow(dwith) ) , ]
if( nrow(dwithout) > nsamplepoints ){dwithout <- sample_n( dwithout, nsamplepoints ) }
dsample <- rbind( dwith, dwithout )

dsample$foreign <- dsample$nforeign > 0

# PLOT 1 : peptides wwo foreign neighbors, plot exchangeability
#p1 <- ggplot( dsample, aes( x = foreign, y = exch, group = nforeign ) ) +
#	geom_quasirandom( size = 0.1 ) +
#	labs( x = "foreign neighbors?",y  = "exchangeability" ) +
#	mytheme

#ggsave(paste0(outfile,"-violin.pdf"), width = 5, height = 5, units = "cm" )

# PLOT 2 : binned exchangeability vs % of peptides with a foreign neighbor
qq <- quantile( d$exch, prob = seq(0,1, length.out=nbins+1 ) )
d$bin <- sapply( d$exch, function(x) which( x <= qq )[1] - 1 )
d$bin[d$bin==0] <- 1
bnames <- names(qq)[-1]#sapply( 1:(length(qq)-1), function(x) paste0( round(qq[x],1),"-",round(qq[x+1],1) ) )
d$bname <- bnames[d$bin]

d2 <- d %>% 
	group_by( bin ) %>%
	summarise( nwithforeign = sum( nforeign > 0 ), ntot = n() ) %>%
	as.data.frame()

d2$pwithforeign <- 100*d2$nwithforeign/d2$ntot
#d2$perc <- as.numeric( sub( "%","",d2$bname ) )
p2 <- ggplot( d2, aes( x = bin, y = pwithforeign ) ) +
	#geom_line( show.legend = FALSE ) + 
	geom_bar( stat = "identity", fill="black" ) +
	scale_x_continuous( breaks=seq(1,nbins))+#limits=c(0,100),expand=c(0,0),breaks=c(0,20,40,60,80,100)) +
	labs( x = "exchangeability\n(decile)", #) +
	#labs( x = "exchangeability",
		y = "% resembling HIV" ) +
	mytheme + theme(
		axis.line.x = element_blank()
	)

ggsave(paste0(outfile,".pdf"), width = 4.7, height = 4, units = "cm" )