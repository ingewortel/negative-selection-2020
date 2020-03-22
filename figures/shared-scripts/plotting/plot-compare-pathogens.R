library(ggplot2)
suppressMessages(library(dplyr))
source("../shared-scripts/plotting/mytheme.R")

argv <- commandArgs( trailingOnly = TRUE )
ntrain <- as.numeric( argv[1] )
datafile <- argv[2]
outfile <- argv[3]

d <- read.table( "data/precursors-bias-calc.txt", header=TRUE )

d <- d[ d$r == 4 & d$ntrain == ntrain, ]

# show sem, not sd
d$lo_top <- d$m_top - d$sem_top
d$hi_top <- d$m_top + d$sem_top

data <- d %>% select( lang, m_top, lo_top, hi_top, type )
ncats <- length( unique( data$lang ) ) 

# remake x axis
pathogen2num <- setNames( seq( 1, ncats )*5 - 3, levels( data$lang ) )
num2pathogen <- unlist( lapply( levels(data$lang), function(x) rep(x,5) ))
pathogenlookup <- c( ebola="Ebola", hcmv="HCMV", hepb="Hep. B", hepc="Hep. C", hiv="HIV", lis="Listeria", mal="Malaria", vac="Vaccinia", zika="Zika")
num2pathogen <- pathogenlookup[ num2pathogen ]

data$lang2 <- as.numeric( data$lang )*5 - 3
data$lang2[ data$type == "R" ] <- data$lang2[ data$type == "R" ] - 0.7
data$lang2[ data$type == "B5" ] <- data$lang2[ data$type == "B5" ] + 0.7

x_breaks <- seq( 2, max( data$lang2 ) , by = 5)
x_labels <- num2pathogen[x_breaks]


ggplot( data, aes( x = lang2, y = m_top, fill = type ) ) +
	annotate( "segment", x=0,xend=max( data$lang2) + 1, y=50, yend=50, size=0.3) + 
	geom_point( aes( colour=type )) + 
	#geom_bar( stat = "identity" , colour="black", width=0.6, size=0.3) +
	scale_x_continuous( limits=c(0,max(data$lang2)+1), expand=c(0,0), breaks=x_breaks, labels=x_labels ) +
	geom_errorbar( aes( ymin=lo_top, ymax = hi_top ), width=0.6, size=0.3 ) +
	scale_y_continuous( limits=c(25,101), expand=c(0,0), breaks=c(25,50,75,100) ) +
	scale_colour_manual( values=c(R="gray70",B1="black", B5="dodgerblue3")) + 
	guides( fill = FALSE , colour = FALSE ) +
	labs( x= "", y = "% foreign in top 10%") +
	mytheme +
	theme(
		axis.text.x=element_text(angle=0,hjust=0.5, vjust=0.5),
		axis.title.x=element_blank(),
		axis.line.x=element_blank() )

title <- outfile #paste0("plots/pathogens-bias-n", ntrain, ".pdf")
ggsave( title, width = 12.5, height= 4, units="cm")
