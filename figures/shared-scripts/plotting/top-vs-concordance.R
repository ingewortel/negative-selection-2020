source("../shared-scripts/plotting/mytheme.R")
source("../shared-scripts/plotting/repositionlabel.R")
library( ggplot2 )

argv <- commandArgs( trailingOnly = TRUE )

concordance_file <- argv[1]
calc_file <- argv[2]
control_file <- argv[3]
n <- as.numeric( argv[4] )
rvalue <- as.numeric( argv[5] )
outfile <- argv[6]

dconc <- read.table( concordance_file, header=TRUE )
dcalc <- read.table( calc_file, header=TRUE )
dcontrol <- read.table( control_file, header=TRUE)

# Combine concordance and model data
data <- data.frame()
for( class in unique( dcalc$lang ) ){
  concordance <- dconc$conc[dconc$class==class]
  top10 <- dcalc$m_top[ dcalc$r==rvalue & dcalc$lang==class & dcalc$ntrain==n ]
  
  
  data <- rbind( data, data.frame( class= class, concordance=concordance, top=top10 ) )
  
}

# Remove en_t2
data <- data[data$class!="en_t2",]

# Add control data
control.data <- data.frame( class="self", concordance=mean(dcontrol$concordance), top = mean(dcontrol$top))
data <- rbind(data, control.data)
pearson_r <- round( cor.test( data$concordance, data$top )$estimate, 3)

sink("data/topvsconc_correlation_coefficient_confint.txt")
print(cor.test(data$concordance, data$top)$conf.int)
sink()

lookup_classes <- c(xh="Xhosa",ta="Tagalog",hi="Hiligaynon",
                    la="Latin",pd="Plautdietsch",me="Medieval English",
                    en="English (different book)",
                    en2="English (different book)",
                    self="English")

data$label <- lookup_classes[as.character(data$class)]
data <- data[ !(data$class %in% c("en") ),]

#shift label position
data <- data[order(data$top),]
miny <- 38
data$labely <- data$top
if( data$labely[1] < miny ){ data$labely[1] <- miny }
data$labely <- repositionlabel(data$labely,3.5)



#data$yshift=0
#data$yshift[data$class=="la"] <- 0.7
#data$yshift[data$class=="pd"] <- -0.7
#data$yshift[data$class=="ta"] <- 0.7
#data$yshift[data$class=="hi"] <- -0.7
#data$yshift[data$class=="self"] <- -1.6

lookup_colours <- setNames(rep("black",length(unique(data$class))), unique(data$class))
lookup_colours["xh"] <- pathogencolors["xh"] # pathogencolors defined in ../scripts/plotting/mytheme.R
lookup_colours["self"] <- pathogencolors["self"]
#lookup_colours["la"] <- pathogencolors["la"]

p <- ggplot( data, aes( x = concordance,
                   y = top,
                   label = label, 
                   colour=class)) +
    #stat_ellipse( data=dcontrol, inherit.aes=FALSE, aes(x=concordance,y=top), lty=2, colour=pathogencolors["self"], size=0.3, level=0.67, segments=100, type="t")+
    #stat_smooth(method = "lm", color="black", size=0.3, lty=2, se=FALSE) +
    geom_point( size = 1) +
    #geom_label(aes( y = top + yshift ), size=2.7,nudge_x=0.02, hjust=0, label.padding=unit(0.1,"lines"),colour="white" ) +
    geom_text(aes( y = labely ), size=2.7,nudge_x=0.02, hjust=-0 ) +   
    scale_x_continuous( limits=c(0.45,1), breaks=seq(0.5,1,by=0.1), expand=c(0,0)) +
    scale_y_continuous( limits=c(45,100 ) , expand=c(0,0)) +
    scale_colour_manual( values=lookup_colours) +
    annotate("text",x=0.48,y=96,label=paste("Pearson's r =",pearson_r),hjust=0,size=2.7) +
    guides( colour=FALSE )+
    labs( x="concordance", y = "% foreign in top 10%") +
    mytheme

ggsave(outfile, width=6,height=6,units="cm")
