source("../shared-scripts/plotting/mytheme.R")
library( ggplot2 )
library( ggbeeswarm )
library( dplyr )

argv <- commandArgs( trailingOnly = TRUE )

selffile <- argv[1]
foreignfile <- argv[2]
nsim <- as.numeric( argv[3] )
r <- as.numeric( argv[4] )
outname <- argv[5]

print(r)

dself <- read.table( selffile )
dself <- head( dself, nrow(dself)/nsim )
dforeign <- read.table( foreignfile )
dforeign <- head( dforeign, nrow(dforeign)/nsim )

d <- rbind( dself, dforeign )
colnames(d) <- c( "id", "string", paste0("r", seq(1,6) ) )

# keep only count columns
d2 <- d[,-c(1,2)]

# these columns contain the number of TCRs with affn >= r.
# Compute the number of TCRs/million with affn == r:
d3 <- d2
for( i in seq(1,ncol(d2)-1) ){
	d3[,i] <- ( d2[,i] - d2[,i+1] )
}

# check if peptides are in the top 10% of best-recognized peptides
# under the current threshold.
qq <- quantile( d[,r], seq( 0, 1, length.out=11 ) )
d$top <- ifelse( d[,r] >= qq[10], "top 10%", "bottom 90%" )
d3$top <- d$top

# count sums per group
dsummary <- d3 %>%
	group_by( top ) %>%
	summarise( 	r1 = sum( as.numeric(r1) ), r2 = sum( as.numeric(r2) ), r3 = sum( as.numeric(r3) ), 
			r4 = sum( as.numeric(r4) ), r5 = sum( as.numeric(r5) ), r6 = sum( as.numeric(r6) )  ) %>%
	as.data.frame()

dplot <- data.frame()
for( ri in 1:6 ){
	dtmp <- data.frame( top = dsummary$top, affn = ri )
	dtmp$count <- dsummary[ , paste0("r",ri)  ]
	dplot <- rbind( dplot, dtmp)
}

totals <- setNames( c(0,0), unique( dplot$top ) )
for( group in unique( dplot$top ) ){
	totals[ group ] <- sum( dplot$count[ dplot$top == group & dplot$affn >= r ] ) 
}
print(totals)

# counts per group per affinity divided by total counts per group
dplot <- dplot %>%
	filter( affn >= r ) %>% 
	group_by( affn, top ) %>%
	summarise( count = count, csum = totals[top], freq = (count)/csum ) %>%
	as.data.frame()

print(dplot)
	

# compute mean affinity from this
#d$affn <- apply( d3, 1, function(x) sum( x*seq(1,6) )/sum(x) )
#meantop <- mean( d$affn[d$top=="top 10%"] )
#meanbottom <- mean( d$affn[ d$top=="bottom 90%"] )

ybreaks <- c(0,3,6)
ylabels <- sapply( ybreaks, function(x) as.expression( bquote( 10^.(x) ) ) )

# plot affinities
p <- ggplot( dplot, aes( x = affn, y = 10^6*freq, group = top, fill = top ) ) +# ggplot( d, aes( x = top, y = affn, group = top, color = id )) +
	geom_bar( stat = "identity", position="dodge", color = "black", size = 0.3 ) +
	labs( x = "affinity", y = "TCRs/million" ) +
	scale_y_log10( limits = c( 1, 10^6) ,breaks = 10^ybreaks, labels = ylabels ) + 
	scale_fill_manual( values = c( "top 10%" = "white", "bottom 90%" = "black" ) ) + 
	#geom_quasirandom(show.legend=FALSE,size=0.3) +
	#scale_x_discrete( labels=c("bottom\n90%","top\n10%") )+
	#annotate("segment", x = 2 - 0.3, 
	#	xend = 2 + 0.3, y = meantop, yend = meantop ) +
	#annotate("segment", x = 1 - 0.3, 
	#	xend = 1 + 0.3, y = meanbottom, yend = meanbottom ) +
	#scale_color_manual( values = pathogencolors ) +
	#labs( x = " ", y = "mean TCR affinity" ) +
	mytheme + theme(
		legend.position = "right"	
	)

ggsave(outname,width=8,height=5,unit="cm")