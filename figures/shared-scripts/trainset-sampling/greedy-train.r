argv <- commandArgs( trailingOnly = TRUE )


greedyfile <- argv[1]
testfile <- argv[2]
outfile <- argv[3]

library( data.table )


greedy <- fread( greedyfile, data.table = FALSE, header = FALSE )
test <- fread( testfile, data.table = FALSE, header = FALSE )

greedy <- as.character( greedy$V1 )
greedy <- greedy[ !is.element(greedy,test$V2) ]


write.table( greedy, file=outfile, row.names=FALSE, col.names=FALSE, quote=FALSE )
