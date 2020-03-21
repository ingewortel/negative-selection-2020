library( ggplot2, warn.conflicts=FALSE )
source( "../shared-scripts/plotting/mytheme.R")
argv <- commandArgs( trailingOnly = TRUE )

reffile <- argv[1]
testfile <- argv[2]
outfile <- argv[3]


# Read files
ref_peps <- as.character( read.table( reffile )$V1 )
test_peps <- as.character( read.table( testfile )$V1 )

# split up in a single vector of amino acids
ref_aa <- unlist( strsplit( ref_peps, "" ) )
test_aa <- unlist( strsplit( test_peps, "" ) )

# count the amino acids
ref_table <- table( ref_aa )
test_table <- table( test_aa )

# convert to frequencies
ref_freq <- ref_table / sum( ref_table )
test_freq <- test_table / sum( test_table )

# compute enrichment
AA <- sort( unique( c( names(ref_freq), names(test_freq ) ) ) )
enrichment <- sapply( AA, function(x) log( test_freq[x] / ref_freq[x] ) )


# To a dataframe for plotting
d <- data.frame(
		AA = factor( AA, levels = AA[ order( ref_freq ) ] ),
		enrichment = enrichment
	)

d <- d[ d$AA != "X", ]

ylab <- ""

ggplot(d, aes( x = AA, y = enrichment ) ) +
	geom_bar( stat="identity" , fill="black") +
	labs( 	x = "AAs (from rare to common)", y = ylab) +
	scale_y_continuous( limits=c(-0.4,0.4), expand=c(0,0) ) +
	mytheme +
	theme(
		axis.line.x = element_blank() 
	)

# width was 6.3 - 7.6
ggsave( outfile, width = 6.3, height=4, units="cm")
