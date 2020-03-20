suppressMessages( library( dplyr )  )
source("../shared-scripts/analysis/concordance-functions2.R")

argv <- commandArgs( trailingOnly = TRUE )

print(argv)

testepdir <- argv[1]
rvalue <- as.numeric( argv[2] )
nsim <- as.numeric( argv[3] )
outfile <- argv[4]

# Function to extract r-mers from strings
make.rmermatrix <- function(peptides,rvalue){
  
  n.rmers <- nchar(peptides[1])-rvalue+1 #Number of rmers per peptide (depending on peptide length and the value of r)
  n.peps <- length(peptides)
  
  
  #Function to get all possible rmers out of a peptide and express them as p[X]r, where p is the position number where the rmer starts and [X]r the rmer sequence.
  get.rmers <- function(peptide,r){
    pep.length <- nchar(peptide)        #peptide length
    max.position <- pep.length - r + 1  #last position where an rmer can start in the peptide
    positions <- seq(1,max.position)    #all possible starting positions of rmers
    
    #Get the substrings of length r for each possible starting position:
    rmers <- sapply(positions,function(x) paste0(x,substring(peptide,x,x+r-1)))
    rmers
  }
  
  #Use get.rmers() function to create one character vector with all rmers from all peptides (loop over peps with sapply)
  all.rmers <- as.character(sapply(peptides,get.rmers,rvalue))
  
  #Generate a matrix where each peptide is defined as a set of rmers (and each rmer is an integer number)
  un.rmers <- unique(all.rmers)                             #Get the vector of unique rmers (one rmer may occur in various peptides)
  rmer.indices <- 1:length(un.rmers)                        #Generate an integer vector of the same length
  names(rmer.indices) <- un.rmers                           #Link this integer vector to the unique rmers as a lookup table.
  int.rmers <- rmer.indices[all.rmers]                      #Translate the 'all.rmers' to an integer vector using this lookup table.
  
  pep.subsets <- matrix(int.rmers,ncol=n.rmers,byrow=TRUE)
  pep.subsets
}


# loop over classes & sims
data <- data.frame()

for( s in 1:nsim){

  if(nsim==1){
    testepfile <- testepdir
    verbose <- TRUE
  } else {
    testepfile <- paste0( testepdir, "/all-sim",s,".txt")
    verbose <- TRUE
  }
  
    # Read the file with strings
    d <- read.table( testepfile )

    sim.data <- concordances(d, rvalue, "self", verbose)
    print('hi2')
    sim.data$sim <- s
    
    #if( anyNA( sim.data$concordance )){
    #  print(sim.data)
    #  concordances(d,verbose=TRUE)
    #}
    
    # Add to dataframe
    data <- rbind( data, sim.data)
}

# Mean and sem over sims
data <- data %>%
  group_by( class ) %>%
  summarise( conc = mean(concordance), sd=sd(concordance),sem=sd(concordance)/sqrt(n()) )

write.table( data, file=outfile, quote=FALSE, row.names=FALSE )
