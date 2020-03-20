#DESCRIPTION:
#Set of functions to choose an optimal set of self-peptides:
#     -   make.rmermatrix           generates a matrix of rmers (integer numbers), with a line for each peptide in self and a column for all the positions.
#                                     INPUT:
#                                         o   selffile        file with self-peptides to optimize for
#                                         o   rvalue          the value of r (or k) to use  
#                                     OUTPUT: 
#                                         o   rmermatrix
#     -   construct.problem         also makes this rmermatrix, but then converts the problem to a sparse matrix M. Also makes an (empty) solution vector s.
#                                     INPUT:
#                                         o   selffile        file with self-peptides to optimize for
#                                         o   rvalue          the value of r (or k) to use  
#                                     OUTPUT: 
#                                         o   list with $rmermatrix, $s and $M
#     -   optimize.setcover         algorithm to simplify the setcover problem. Iteratively (1) includes peptides that cover unique rmers, (2) removes peptides that are a subset of others.
#                                     INPUT:
#                                         o   M               sparse matrix of the problem to optimize for
#                                         o   solution        vector of size npeps to store solution in (can already contain info)
#                                         o   rvalue          the value of r (or k) to use  
#                                     OUTPUT: 
#                                         o   list with $M (the matrix of the simplified problem), $s (the solution vector with the new information on included/excluded peptides) after optimization
#     -   mps.from.sparseM          generates MPS file based on matrix M provided. 
#                                     INPUT:
#                                         o   M               sparse matrix of the problem to optimize for
#                                         o   outname         name of the mps file to write, eg "test.mps"
#                                         o   rvalue          the value of r (or k) to use  
#                                     OUTPUT: 
#                                         o   value NULL, but also generates the .mps file. 
#     -   greedyPeptides2           solves problem with greedy algorithm.
#                                     INPUT:
#                                         o   rmermatrix      see make.rmermatrix
#                                         o   threshold       rmers with a frequency < threshold will not (necessarily) be covered.
#                                     OUTPUT: 
#                                         o   solution        solution vector after the greedy algorithm. 
#     -   greedyTolerance           returns a solution of N peptides that maximizes tolerance
#                                     INPUT:
#                                         o   selffile        file with self-peptides to optimize for
#                                         o   rvalue          the value of r (or k) to use  
#                                         o   peps            the number of peptides the solution may contain
#                                     OUTPUT: 
#                                         o   list with $s, the solution vector, and data, a dataframe with included peptides and tolerance achieved after each step.


#GENERATING MATRIX M
make.rmermatrix <- function(selffile,rvalue){
  #########################
  #READ DATA & EXTRACT INFO
  #########################
  peptides <- as.character(read.table(selffile)$V1) 
  n.rmers <- nchar(peptides[1])-rvalue+1 #Number of rmers per peptide (depending on peptide length and the value of r)
  n.peps <- length(peptides)
  
  #########################
  #CALCULATE RMERS AND CORRESPONDING INTEGERS
  #########################
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
construct.problem <- function(selffile,rvalue){
  library(Matrix)
  
  #For conversion to a sparse matrix: put into a table with peptide in the first column and rmer in the second
  rmermatrix <- make.rmermatrix(selffile,rvalue)
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

#OPTIMIZATION:
optimize.setcover <- function(M,solution,rvalue){
  library(Matrix)
  op.step <- 0
  
  message("\nStarting optimization ... \n ")
  
  rmers.per.pep <- 7-rvalue
  
  begin.M <- sum(M)
  end.M <- 0
  
  while(begin.M!=end.M){
    
    #Progress:
    op.step <- op.step + 1
    begin.problem.rmers <- sum(colSums(M)!=0)
    begin.problem.peptides <- sum(rowSums(M)!=0)
    begin.M <- sum(M)
    
    cat("\n----------")
    message(paste0("\nOptimization step ",op.step,": "))
    cat(paste0("\nOptimization step ",op.step,": "))
    cat(paste0("\nCurrent problem size of ",begin.problem.peptides," peptides and ",begin.problem.rmers," rmers. Beginning optimization."))
    
    
    #STEP 1 ---- First include all peptides that contain unique rmers. Find these uniquely covered rmers:
    start.time <- Sys.time()
    uniquely.covered.rmers <- which(colSums(M)==1)
    cat(paste0("\n  ",length(uniquely.covered.rmers)," rmers occur only in one peptide. Looking for the corresponding peptides...\n"))
    
    included.peps.before <- sum(solution,na.rm=TRUE)
    num.rmers.before <- sum(colSums(M)!=0)
    
    
    #Find the peptides to which these rmers belong:
    for( r in uniquely.covered.rmers ){
      #list the coordinates (in this case only one) in which column r (the rmer we are looking at) is not zero >> get the peptides. 
      p <- summary(M[,r,drop=FALSE])$i              #drop=FALSE so the matrix is kept sparse. i gives the row number.
      #Now find the rmers covered by this peptide and set the corresponding columns to FALSE (remove them from the matrix).
      for( i in summary(M[p,,drop=FALSE])$j ){    
        M[,i] <- FALSE
      }
      #Also set the peptide to TRUE in the solution vector, because the peptide is included.
      solution[p]<-TRUE
    }
    
    included.peps.after <- sum(solution,na.rm=TRUE)
    num.rmers.after <- sum(colSums(M)!=0)
    
    #Print progress
    cat("     >>")
    print(Sys.time()-start.time)
    cat(paste0("     ",included.peps.after-included.peps.before, " new peptides included, now total ",included.peps.after, " peptides included. ",num.rmers.before-num.rmers.after, " rmers have been removed.\n"))
    
    # print( table( colSums(M) ) )
    # print( table( rowSums(M) ) )
    
    #STEP 1B ---- Remove any peptides that are now redundant
    #Look for rows(peptides) that are now empty in the M matrix (=rmers completely covered), but that are still NA in the solution vector.
    #These don't need to be included so set them to FALSE.
    
    cat("  Checking for peptides that are now redundant...\n")
    
    M.covered.peptides <- rowSums(M)==0
    solution.NA.peptides <- is.na(solution)
    redundant.peptides <- M.covered.peptides & solution.NA.peptides
    solution[redundant.peptides] <- FALSE
    
    cat("   ",sum(redundant.peptides),"redundant peptides found & excluded.\n")
    
    #STEP 2 ---- Remove peptides that are now a subset of other peptides:
    start.time <- Sys.time()
    
    cat("  Checking for peptides that are a subset of other peptides...\n")
    
    #track number of peptides removed.
    n.rem <- 0
    i <- 1 #for tracking which peptide is currently being checked.
    iN <- sum(rowSums(M) %in% 1:(rmers.per.pep-1 )) #Number of peptides to check: the peptides that have between 1-ncol(rmermatrix) rmers left.
    
    message(paste0("Comparing ",iN," peptides with their potential supersets, currently checking peptide number: \n"))
    
    #Loop over the the different options for number of rmers left (1-ncol(rmermatrix))
    for( k in 1:(rmers.per.pep-1 )){
      message("\nk=",k, "\n")
      #print(table(rowSums(M)))    
      #Find the rows (peptide) that have this number of rmers left.
      covers.k <- which(rowSums(M)==k)
      
      #Loop over these peptides:
      for( p in covers.k ){
        
        #List the rmers that are covered by this peptide
        covered.rmers <- summary(M[p,,drop=FALSE])$j
        
        #Look for peptides that cover one or more of the same rmers, and count how many of the same rmers they cover.
        covering <- table(summary(M[,covered.rmers,drop=FALSE])$i)
        
        #Now get the number of the peptides for which the number of overlapping rmers is at least k.
        #(so this will give peptides that are either the same, or a superset.)
        other.covers <- as.integer(names(covering[covering==k]))
        
        #You will always find the peptide itself. The peptide is a subset/equal to another peptide only if other.covers contains more than one element.
        if( length( other.covers ) > 1 ){
          #Remove the row from the matrix. Set the peptide to FALSE in the solution, and keep track of the number of removed peptides.
          M[p,] <- FALSE
          solution[p] <- FALSE
          n.rem <- n.rem + 1
        }
        
        #Print progress:
        i <- i + 1
        if( i %% 1000 == 0 ){
          message("\n         ",i, " of ", iN)
        }
      }
    }
    
    end.M <- sum(M)
    
    #Track the time this step took:
    cat("\n     >>")
    print(Sys.time() - start.time)
    
    #Print some concluding info.
    cat(paste0("     ",n.rem," peptides are a subset of others and have been excluded."))
    
    #print progress:
    
    end.problem.rmers <- sum(colSums(M)!=0)
    end.problem.peptides <- sum(rowSums(M)!=0)
    cat(paste0("\nProblem reduced to ",end.problem.peptides," peptides and ",end.problem.rmers, " rmers."))
    cat(paste0("\nTotal: ",sum(solution==TRUE,na.rm=TRUE)," peptides included, ",sum(solution==FALSE,na.rm=TRUE)," excluded and ",end.problem.peptides, " to go."))
  }
  
  cat("\n__________________________________")
  
  if(sum(M)==0){
    cat("\n Solution found! ",sum(solution==TRUE,na.rm=TRUE)," peptides needed.")
  } else {
    cat("\n Problem cannot be simplified any further. ",end.problem.peptides," peptides and ",end.problem.rmers, " rmers remaining.")
  }
  
  
  
  list(M=M,s=solution)
}

#MPS FILE FOR CBC
mps.from.sparseM <- function(M,filename,rvalue){
  library(Matrix)

  mpsname <- paste0("R",rvalue,"N",nrow(M))
  mps.header <- sprintf("%-13s %s","NAME",mpsname)
  
  #To generate the rows & rhs section: loop over the rmers (one constraint for each rmer: >0. so add a ROW for each rmer with "G"(greater than) & rmername, and an RHS of 0)
  #column section & bounds section: loop over peptides, add a line for each rmer they contain (=each constraint they are in) and one extra line for the objective function.
  
  #Initialize document:
  rows <- "\nROWS"
  rows <- c(rows,sprintf("\n %-3s%s","N","obj")) #first line for the objective function
  sink(file=filename)
  cat(mps.header,rows)
  
  #all the rmers for which constraints should be added (only rmers that haven't been covered during any pre-optimization)
  rmer.indices <- which(colSums(M)!=0)
  
  #Rows section: one constraint per rmer (rmer > 0), so each line has a "G" (greater than) and the rmer number.
  for(i in rmer.indices){
    rmername <- paste0("r",i)
    newline.rows <- sprintf("\n %-3s%s","G",rmername)
    cat(newline.rows)
  }
  
  #Column section: loop over peptides and add a line for all the rmers they contain (each constraint they are in), and extra line for the objective function.
  cols <- "\nCOLUMNS"; cat(cols)
  for(i in 1:nrow(M)){
    pepname <- paste0("p",i)
    rmerset <- summary(M[i,,drop=FALSE])$j
    newlines.cols <- sapply(rmerset,function(x) sprintf("\n    %-10s%-10s%d",pepname,paste0("r",x),1))
    newline.cols.obj <- sprintf("\n    %-10s%-10s%d",pepname,"obj",1)
    cat(newlines.cols,newline.cols.obj)
  }
  
  #RHS section: one constraint per rmer, so add a one for each constraint.
  rhs <- "\nRHS"; cat(rhs)
  
  for(i in rmer.indices){
    rmername <- paste0("r",i)
    newline.rhs <- sprintf("\n    %-10s%-10s%d","rhs",rmername,1)
    cat(newline.rhs)
  }
  
  #Bounds section: indicate for each peptide that it has to be binary.
  bounds <- "\nBOUNDS"; cat(bounds)
  for(i in 1:nrow(M)){
    pepname <- paste0("p",i)
    newline.bounds <- sprintf("\n BV %-10s%-10s","BND",pepname)
    cat(newline.bounds)
  }
  
  #the final line & write the file
  mps.end <- "\nENDATA"
  cat(mps.end)
  
  sink()
  
  
}

#GREEDY ALGORITHMS
greedyPeptides2 <- function(rmermatrix,threshold=1){
  cat("\nStarting Greedy Algorithm\n")
  solution <- rep(FALSE,nrow(rmermatrix))
  
  #preprocessing of rmermatrix based on threshold: set all rmers that have a too low frequency to NA beforehand.
  rmerfreq <- table(as.numeric(rmermatrix))
  rmernames <- as.numeric(dimnames(rmerfreq)[[1]])
  
  rmers.below.threshold <- rmernames[which(rmerfreq<threshold)]
  rmermatrix[is.element(rmermatrix,rmers.below.threshold)]<-NA
  
  low.bound <- max(apply(rmermatrix,2,function(x) length(unique(x[!is.na(x)]))))
  
  
  
  iteration <- 0
  while(sum(!is.na(rmermatrix))!=0){
    
    iteration <- iteration + 1
    included.peptides.begin <- sum(solution,na.rm=TRUE)
    
    if(iteration%%10==0){
      cat("Iteration ",iteration,": ",included.peptides.begin," peptides included, ",100*sum(!is.na(rmermatrix))/(ncol(rmermatrix)*nrow(rmermatrix)),"% of positions left to cover. \n")
    }
    
    
    
    #which peptides contain the most rare rmers?
    rmerfreq2 <- table(as.numeric(rmermatrix,na.rm=TRUE))
    rmernames2 <- as.numeric(dimnames(rmerfreq2)[[1]])
    rarest.rmers <- rmernames2[which(rmerfreq2==min(rmerfreq2,na.rm=TRUE))]
    
    rarest.positions <- matrix(is.element(as.numeric(rmermatrix),rarest.rmers),ncol=ncol(rmermatrix))
    rarecoverpeps <- (rowSums(rarest.positions,na.rm=TRUE)>0)
    
    
    
    if(min(rmerfreq2,na.rm=TRUE)==1){
      #include all peptides with a unique rmer
      covered.rmers <- as.numeric(rmermatrix[rarecoverpeps,])
      which.covered <- matrix(is.element(as.numeric(rmermatrix),covered.rmers),ncol=ncol(rmermatrix))
      rmermatrix[which.covered]<-NA
      solution[rarecoverpeps] <- TRUE
    } else{
      
      
      #which values are NA for these peptides? find the peptides with the least NA (so that cover the most uncovered peptides)
      NAmatrix <- is.na(rmermatrix[rarecoverpeps,])
      maxcoverpeps <- rmermatrix[rarecoverpeps,][rowSums(NAmatrix)==min(rowSums(NAmatrix,na.rm=TRUE)),]
      maxcoverpeps <- matrix(maxcoverpeps,ncol=ncol(rmermatrix))
      
      
      #find peptides that do not contain any duplicated rmers (NAs don't count as duplicated)
      dup.matrix <- matrix(duplicated(as.numeric(maxcoverpeps)),ncol=ncol(rmermatrix))
      dup.matrix[is.na(maxcoverpeps)]<-FALSE
      
      #now take peptides that do not have any duplicated rmers. If there are none, take one peptide with as few as possible duplicated rmers. 
      include.peps <- rowSums(dup.matrix)==0
      if(sum(include.peps,na.rm=TRUE)==0){
        cat("Warning: no peptides without duplicated rmers! There should be at least one... \n")
        include.peps <- rep(FALSE,nrow(dup.matrix))
        include.peps[which(rowSums(dup.matrix)==min(rowSums(dup.matrix)))[1]]<-TRUE
      }
      
      #check which rmers are covered by the now chosen set, change these to NA in the matrix and set the peptides to TRUE in the solution. 
      
      covered.rmers <- as.numeric(maxcoverpeps[include.peps,])
      which.covered <- matrix(is.element(as.numeric(rmermatrix),covered.rmers),ncol=ncol(rmermatrix))
      rmermatrix[which.covered]<-NA
      solution[rarecoverpeps][rowSums(NAmatrix)==min(rowSums(NAmatrix,na.rm=TRUE))][include.peps] <- TRUE
    }
    
    included.peptides.after <- sum(solution,na.rm=TRUE)
    if(included.peptides.begin==included.peptides.after){
      stop("No new peptides included in this iteration. Something went wrong.")
    }
  }
  
  cat("\n Done. ",sum(solution,na.rm=TRUE)," peptides included after ",iteration," iterations. (Low bound: ",low.bound,")")
  solution
  
}
greedyTolerance <- function(selffile,rvalue,peps){
  
  #convert the problem:
  message("Constructing problem matrices")
  
  p <- construct.problem(selffile,rvalue)
  M <- p$M
  s <- rep(FALSE,nrow(M))
  sv <- rep(0,nrow(M)) #for tracking the step in which peptides are included. 0 if not included at all.
  rmerm <- p$rmermatrix
  
  #empty dataframe to track data
  df <- data.frame()
  
  #Loop while the number of peptides included is lower than 'peps' (and while there is no 100%tol yet)
  step <- 0
  while((sum(s)<peps) & (sum(M)!=0)){
    step <- step + 1
    
    
    #1    how many peps left to include?
    peps.needed <- peps-sum(s)
    
    
    #3    calculate rmerfreqs & pepscores
    
    rmerfreq <- colSums(M)
    names(rmerfreq)<-seq(1,ncol(M))
    
    #create 'frequency matrix' from rmermatrix rm by 'looking up' the frequencies in the rmerfreq vector
    freqm <- rmerm
    for(i in seq(1,ncol(rmerm))){#build columnwise
      freqm[,i] <- rmerfreq[rmerm[,i]]
    }
    #now calculate the pepscores with rowSums:
    pscores <- rowSums(freqm)
    
    
    #4    Select highest scoring peptides
    if(sum(pscores==max(pscores))>1){
      chunks <- sum(pscores==max(pscores))
    } else {
      chunks <- 2
    }
    highscores <- head(order(pscores,decreasing=TRUE),chunks)
    
    #5    Remove peptides with duplicate rmers from the highscores set and include the rest
    
    #build rmermatrix for the subproblem
    rmerm0 <- rmerm[highscores,]
    
    #find peptides with duplicated rmers, include the ones that don't have any.
    dupM0 <- apply(rmerm0,2,duplicated)
    choosepeps <- (rowSums(dupM0)==0)
    
    #find back indices. Include them all or the top ones (if set is larger than needed.)
    pep.indices <- highscores[choosepeps]
    
    if(length(pep.indices)>peps.needed){
      pick.peps <- head(pep.indices,peps.needed)
    } else {
      pick.peps <- pep.indices
    }
    
    covered.rmers <- unique(as.numeric(rmerm[pick.peps,]))
    s[pick.peps]<- TRUE
    M[,covered.rmers]<-FALSE
    sv[pick.peps] <- step
    
    rmertol <- round(100*sum(colSums(M)==0)/ncol(M),2)
    peptol <- round(100*sum(rowSums(M)==0)/nrow(M),2)
    
    #add info to df
    df <- rbind(df,data.frame(step=step,npeps=sum(s),rtol=rmertol,ptol=peptol))
    
    #Print progress:
    if(step%%10==0){
      message("Step ",step,": ",sum(s)," peptides included, rmertol=",rmertol,"%, peptol=",peptol,"%")
    }
  }
  message("Done in ",step," steps: ",sum(s)," peptides included, rmertol=",rmertol,"%, peptol=",peptol,"%\n")
  list(s=s,data=df,sv=sv)
}

