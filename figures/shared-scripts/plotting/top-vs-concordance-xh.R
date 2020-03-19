source("../shared-scripts/plotting/mytheme.R")
source("../shared-scripts/plotting/repositionlabel.R")
library( ggplot2 )

argv <- commandArgs( trailingOnly = TRUE )

concordance_file <- argv[1]
calc_file <- argv[2]
n <- as.numeric( argv[3] )
outname <- argv[4]

dconc <- read.table( concordance_file, header=FALSE )
colnames(dconc) <- c("r","conc")
dcalc <- read.table( calc_file, header=TRUE )

# Combine concordance and model data
data <- data.frame()
for( r in unique( dconc$r ) ){
  concordance <- dconc$conc[dconc$r==r]
  top10 <- dcalc$m_top[ dcalc$r==r & dcalc$lang=="xh" & dcalc$ntrain==n ]
  
  
  data <- rbind( data, data.frame( r=r, concordance=concordance, top=top10 ) )
  
}


data$label <- paste( "t =", data$r)

#shift label position
data <- data[order(data$top),]
miny <- 38
data$labely <- data$top
if( data$labely[1] < miny ){ data$labely[1] <- miny }
data$labely <- repositionlabel(data$labely,3.5)

data <- data[order(data$r),]

data$xend <- c(data$concordance[2:nrow(data)],NA)
data$yend <- c(data$top[2:nrow(data)],NA)

p <- ggplot( data, aes( x = concordance,
                   y = top,
                   label = label)) +
    geom_segment( aes( xend=xend,yend=yend), color=pathogencolors["self"]) +
    geom_point( size = 1) +
    geom_text(aes( y = labely ), size=2.7,nudge_x=0.02, hjust=-0 ) +   
    scale_x_continuous( limits=c(0.45,1), breaks=seq(0.5,1,by=0.1), expand=c(0,0)) +
    scale_y_continuous( limits=c(45,100 ) , expand=c(0,0)) +
    #scale_colour_manual( values=lookup_colours) +
    guides( colour=FALSE )+
    labs( x="concordance", y = "% Xhosa in top 10%") +
    mytheme

ggsave(outname, width=6,height=6,units="cm")
