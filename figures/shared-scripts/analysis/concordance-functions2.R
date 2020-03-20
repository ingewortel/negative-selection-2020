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


# Function to create sparse matrix
construct.problem <- function(rmermatrix,rvalue){
  library(Matrix)
  require(methods)
  
  # For conversion to a sparse matrix: put into a table with 
  # peptide in the first column and rmer in the second
  rmervec <- as.vector(t(rmermatrix))
  pepvec <- numeric(length(rmervec))
  for(i in 1:nrow(rmermatrix)){
    start <- ncol(rmermatrix)*(i-1)+1
    pepvec[seq(start,start+ncol(rmermatrix)-1)]<-rep(i,ncol(rmermatrix))
  }
  rmertable <- matrix(c(pepvec,rmervec),ncol=2)
  
  #Convert to sparse matrix
  M <- sparseMatrix( rmertable[,1], rmertable[,2] )
  solution <- rep(NA,dim(rmermatrix)[1])
  
  list(rmermatrix=rmermatrix,M=M,s=solution)
}

neighbor.table <- function( string.data, rvalue ){
	
	# get labels and peptides/strings
	classes <- as.character( string.data$V1 )
	strings <- as.character( string.data$V2 )

	# Create the rmermatrix
	rm <- make.rmermatrix( strings, rvalue )
	
	# Create the sparse matrix
	M <- construct.problem( rm, rvalue )$M

	# Create a table to list for each string:
	# - id
	# - label
	# - # neighbors of each class
	unique.classes <- unique(classes)
	n.classes <- length( unique.classes )
	n.strings <- length( strings )
	neighbors <- matrix(0, nrow = n.strings, ncol = n.classes )
	colnames(neighbors) <- unique.classes

	# loop over strings, check neighbors and add them to table
	for( s in 1:nrow(M) ){

		# First list the rmers of this string
    		rmers <- summary(M[s,,drop=FALSE])$j
    
  		# For every of these rmer, find all other peptides that contain it:
    		rmer.strings <- lapply( rmers, function(x) summary(M[,x,drop=FALSE])$i )
    
    		# Neighborhood: take the unique union of the strings in this list.
    		# Don't count the current string
    		neighborvec <- unique( unlist( rmer.strings ) )
    		neighborvec <- neighborvec[ neighborvec != s ]

		# Find the classes that the neighbors belong to
		neighbor.classes <- classes[ neighborvec ]

		# Make a table
		t <- table( factor( neighbor.classes, levels = unique.classes ) )

		# Add to global table
		neighbors[ s, ] <- t

	}

	neighbor.table <- as.data.frame( neighbors )
	neighbor.table$class <- classes

	return(neighbor.table)

}

normalized.concordance <- function( p, p_exp ){

  # extreme case exception (prevent dividing by 0)
  if( p == 1){

    p_norm <- 1

  } else {
    
    # Correction using proportion of expected
    p_transformed <- log( p / (1-p) ) - log( p_exp / (1 - p_exp ) )
    
    # logit transformation back
    p_norm <- exp( p_transformed ) / ( exp( p_transformed ) + 1 )
    
  }
  
  return( p_norm )
}

neighbor.concordance <- function( neighbor.table, selfclass, foreignclass ){
	
	# check that selfclass and foreignclass occur in the data
	if( !is.element( selfclass, neighbor.table$class ) ){
		stop( paste( "Invalid selfclass. Choose from:", paste( unique(neighbor.table$class), collapse="/" ) ) )
	} else if( !is.element( foreignclass, neighbor.table$class )) {
		stop( paste( "Invalid foreignclass. Choose from:", paste( unique(neighbor.table$class), collapse="/" ) ) )
	}

	# select only the strings of selfclass and foreignclass
	self.strings <- neighbor.table[ neighbor.table$class == selfclass, c(selfclass, foreignclass ) ]
	foreign.strings <- neighbor.table[ neighbor.table$class == foreignclass , c(selfclass, foreignclass ) ]
	nself <- nrow( self.strings )
	nforeign <- nrow( foreign.strings )

	# compute total self + foreign neighbors
	self.strings$total <- rowSums( self.strings )
	foreign.strings$total <- rowSums( foreign.strings )

	# combine again
	self.strings$class <- selfclass
	foreign.strings$class <- foreignclass
	strings <- rbind( self.strings, foreign.strings )

	# compute expected percentage of self and foreign neighbors (equal to the ratio of self and foreign strings)
	p.exp <- setNames( c(nself,nforeign)/(nself+nforeign), c( selfclass, foreignclass ) )

	# convert to proportion of self/foreign neighbors. 
	#If strings have no neighbors at all (total = 0), they get the expected percentage self/foreign neighbors.
	strings[ ,selfclass ] <- ifelse( strings$total == 0 ,
					p.exp[ selfclass ],
					strings[, selfclass]/strings$total )
	strings[, foreignclass ] <- ifelse( strings$total == 0,
					p.exp[ foreignclass ],
					strings[, foreignclass]/strings$total )

	# concordance = proportion of neighbors from the class the string is from
	strings$conc <- ifelse( strings$class == selfclass , strings[ , selfclass ] , strings[ , foreignclass ] )

	# output
	list( conc = strings$conc, nvec = strings$total, pexp = p.exp )

}


concordances <- function(string.data, rvalue, selfclass, verbose = FALSE ){
  
  	
  	nt <- neighbor.table( string.data, rvalue )
	save( nt, file = paste0( "data/concordances/neighbortable-r",rvalue,".Rdata" ) )
	if( verbose ) print('...neighbor table generated')
	classes <- nt$class
	nonself.classes <- unique(classes)[unique(classes)!=selfclass ]
	out.data <- data.frame()

	for( c in nonself.classes ){
		if( verbose ){
			print( paste0( "... concordances ", c ) )
		}
		# compute neighbor concordance per string (unnormalized)
		concdata <- neighbor.concordance( nt, selfclass, c )
		c.classes <- classes[ classes == selfclass | classes == c ]

		# Get mean concordance of self and foreign strings, respectively (still unnormalized)
		self.conc <- mean( concdata$conc[ c.classes == selfclass ] )
		nonself.conc <- mean( concdata$conc[ c.classes == c ] )

		# normalize means based on expectation
		p.exp <- concdata$pexp
		self.norm.conc <- normalized.concordance( self.conc, p.exp[selfclass] )
		nonself.norm.conc <- normalized.concordance( nonself.conc, p.exp[c] )

		# weighted average normalized concordance over classes		
		nself <- sum( classes == selfclass )
		nnonself <- sum( classes == c )
		conc <- ( nself*self.norm.conc + nnonself*nonself.norm.conc ) / (nself + nnonself)

		# output
		out.data <- rbind( out.data, data.frame( class=c, 
				concordance=conc, 
				selfc = self.norm.conc, 
				nonselfc = nonself.norm.conc ) ) 
	}

	# remove rownames and add a column with the number of foreign strings the concordances are
	# based on.
	rownames(out.data) <- NULL
	out.data$n <- as.character( table( classes )[ as.character( out.data$class ) ] )
	
	#print(out.data)
 	return(out.data)
}


