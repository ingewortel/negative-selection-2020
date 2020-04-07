suppressWarnings( library( ggplot2 ) )
suppressWarnings( library( ggbeeswarm ) )
source( "../shared-scripts/plotting/mytheme.R")

# Read command line arguments
argv <- commandArgs( trailingOnly = TRUE )
dir <- argv[1]		# directory to find TCR count file
fileid <- argv[2]	# peps or lang
r <- argv[3]		# affinity threshold t
l <- argv[4]		# name of the language / organism
n <- argv[5]		# number of training sequences
nsim <- as.numeric( argv[6] )	# number of simulations in the file
outfile <- argv[7]		# name of the plot file to generate
maxpf <- as.numeric( argv[8] )	# max value on the y-axis


# Read the file for the current parameter values
filename <- paste0(dir,"/output-r",r,"-",fileid,"-",l,"-n",n,"-r",r,".txt")
d <- read.table( filename )

# nrow(d)/nsim lines for every simulation
lines <- nrow(d)/nsim

# add a column with simulation number
d$sim <- ceiling(seq(1,nrow(d))/lines)

# add the file with self test epitopes
fileself <- paste0(dir,"/output-r",r,"-",fileid,"-self-n",n,"-r",r,".txt")
dself <- read.table( fileself )
dself$sim <- ceiling(seq(1,nrow(dself))/(nrow(dself)/nsim))
d <- rbind(d,dself)

# look at the first sim
d <- d[ d$sim == 1,]

# compute pf (TCRs for this sequence per million total TCRs)
d$pf <- 10^6 * d[,4]/d[,3]

# If the max y-axis value is set to 0, plot just the range of the data
if( maxpf==0 ){ maxpf <- max(d$pf) }

# adjust x and y-axis labels	
ylab <- "reacting TCRs/million"		
if( fileid == "lang" ){
	lookup = c( 	xh="Xhosa",
			self="English" )
	ylab <- "reacting motifs/million"
}
d$name <- factor( lookup[d$V1], levels=lookup[c("self",l)] )
d$cat <- d$V1

# compute medians
median_self <- median( d$pf[d$cat=="self"])
median_foreign <- median( d$pf[d$cat != "self"])

nums <- c(1,2)

# plot
p <- ggplot( d, aes( 	x = name,
			y = pf,
			colour=cat,
			fill = cat ) ) +
	geom_violin( alpha = 0.3, colour=NA ) +
	annotate("segment", x = 2 - 0.3, 
		xend = 2 + 0.3, y = median_foreign, yend = median_foreign ) +
	annotate("segment", x = 1 - 0.3, 
		xend = 1 + 0.3, y = median_self, yend = median_self ) +
	scale_y_continuous( limits=c(0,maxpf), expand=c(0,0) ) +
	geom_quasirandom(size=0.3, groupOnX=TRUE )+
	guides(colour=FALSE,lty=FALSE, fill=FALSE) +
	labs( 	x = " ", 
		y = ylab ) +
	scale_colour_manual( values = pathogencolors ) +
	scale_fill_manual( values= pathogencolors ) +
	mytheme + 
	theme(
		axis.line.x = element_blank() 
	)

ggsave(outfile,width=4.1,height=4.25,unit="cm")
