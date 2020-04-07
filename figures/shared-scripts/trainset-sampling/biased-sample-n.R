# DESCRIPTION:
# Sample a set of N "seen self" peptides based on a probability score.
#
# INPUT
#   o   scorefile     file with the peptides scored on AA frequency (V1=peps, V2=score)
#   o   N             the number of peptides to sample.
#   o   n 	      the power for the score (tunes bias strength)
#  
# OUTPUT
#   o   stdout        the N chosen peptides.

library(data.table)

# 1. COLLECT INPUT
argv <- commandArgs( trailingOnly=TRUE ) 
scorefile <- argv[1]
N <- as.numeric(argv[2])
n <- as.numeric(argv[3])

# 2. CALCULATIONS
# Read the data into a dataframe, extract peptides and scores.
data <- fread(scorefile)
peptides <- as.character(data$V1)
scores <- data$V2

# Raise score to the correct power
scores <- scores^n

# Now generate the biased sample:
peps <- sample(peptides, N, prob=scores)


# 3. WRITE OUTPUT TO STDOUT
cat(peps,sep="\n")
