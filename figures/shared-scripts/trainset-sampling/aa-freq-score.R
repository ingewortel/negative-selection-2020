# DESCRIPTION:
# Calculate a score for each peptide based on the frequencies of the amino acids it contains.
# This score can also be calculated using weights for each position in the peptide to create a
# positional bias. 
#
# INPUT
#   o   selffile      file with the peptides to score
#   o   outfilename   the name of the file to output
#   o   pos.bias      (numeric) vector with weights for each position in the peptide
#  
# OUTPUT
#   o   outfile       file with peptide in first column and score in the second.


# 1. COLLECT INPUT
argv <- commandArgs( trailingOnly=TRUE ) 
selffile <- argv[1]
outfilename <- argv[2]
pos.bias <- as.numeric(unlist(strsplit(argv[3]," ")))
print(pos.bias)

# 2. CALCULATION
# Load file with peptides, convert to a matrix with row for each peptide and column for
# each AA.
peptides <- as.character(read.table(selffile)$V1)
AA.matrix <- t(sapply(peptides,function(x) unlist(strsplit(x,""))))

# Determine AA frequencies for each AA
aa.counts <- table(as.vector(AA.matrix))
aa.freq <- aa.counts/sum(aa.counts)

# New matrix, now with the AA frequency instead of the AA itself in the elements:
freq.matrix <- apply(AA.matrix,2,function(x) aa.freq[x])
for( c in ncol(freq.matrix) ){
	freq.matrix[,c] <- freq.matrix[,c]*pos.bias[c]
}


# Frequency score can now be obtained via rowSums()
freq.score <- rowSums(freq.matrix)


# 3. WRITE OUTPUT
# Make dataframe and sort peptides on frequency score
df <- data.frame(peps=peptides,score=freq.score)
df <- df[order(df$score),]

# Write to file with name 'outfilename'
write.table(df,file=outfilename,quote=FALSE,sep="\t", row.names=FALSE,col.names=FALSE)
