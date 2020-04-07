library( ggplot2 )
source( "../shared-scripts/plotting/mytheme.R" )

argv <- commandArgs( trailingOnly = TRUE )

concfile <- argv[1]
outfile <- argv[2]

d <- read.table( concfile )
colnames(d) <- c("r","conc")


p <- ggplot( d, aes(x = r,
                    y = conc )) +
  geom_col(  position="dodge", fill="black" ) +
  scale_y_continuous( limits=c(0,1), expand=c(0,0), breaks=seq(0.5,1,by=0.1)) +
  scale_x_continuous( breaks=seq(1,6)) +
  coord_cartesian( ylim = c(0.48,1)) +
  labs( x = "affinity threshold (t)", y = "concordance") +
  mytheme +
  theme(
    #axis.title.x= element_blank(),
    axis.line.x=element_blank(),
    #axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)
  ) 


 

ggsave(outfile, plot=p, width=4.2, height=3.4, units="cm")
