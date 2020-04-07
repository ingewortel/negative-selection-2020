# DESCRIPTION:
# Convert AA frequency scores of self peptides to a probability score
#
# INPUT
#   o   scorefile     file with the peptides scored on AA frequency (V1=peps, V2=score)
#   o   outfilename   name of the file to write
#  
# OUTPUT
#   o   outfile       peptides + scores (with n=1)

library( data.table )

# 1. COLLECT INPUT
argv <- commandArgs( trailingOnly=TRUE ) 
scorefile <- argv[1]
outfilename <- argv[2]
epsilon <- as.numeric( argv[3] )

# 2. CALCULATIONS
# Read the data into a dataframe, extract peptides and scores.
data <- fread(scorefile)
peptides <- as.character(data$V1)
scores <- data$V2


# Get the amino acid frequencies again:
AA.matrix <- t(sapply(peptides,function(x) unlist(strsplit(x,""))))
aa.counts <- table(as.vector(AA.matrix))
aa.freq <- aa.counts/sum(aa.counts)

# Calculate minimal and maximal possible score
min.freq <- min(aa.freq)
max.freq <- max(aa.freq)
motif.length <- ncol(AA.matrix)
min.score <- motif.length * min.freq
max.score <- motif.length * max.freq * ( 1 + epsilon )

# Function to calculate score for each peptide
score.to.prob <- function(score,min=min.score,max=max.score){

	# Use formula to calculate probability
	prob <- ( ( max.score - score ) / ( max.score - min.score ) )

}

# Now generate the bias scores
probs <- score.to.prob(scores)

# If probability equals the max probability, score sometimes 
# becomes slightly negative instead of 0 because of imprecision in
# storing numbers. If difference with zero is less than the machine
# epsilon, set prob = 0.
# probs[ abs( probs - 0 )  < .Machine$double.eps ]  <- 0

df <- data.frame(peps=peptides,prob=probs)

# 3. WRITE OUTPUT TO FILE
write.table(df,file=outfilename,quote=FALSE,sep="\t", row.names=FALSE,col.names=FALSE)
