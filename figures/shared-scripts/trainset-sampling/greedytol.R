#DESCRIPTION
#creates an optimal set of size Npeps
#INPUT (from stdin):
#         1   selffile      file with peptides in self set to optimize for
#         2   rvalue        value of r (or k) to use in model
#         3   Npeps         the number of peptides the solution should (maximally) contain
#DEPENDENCIES
#         -   scripts/problem-functions.R                                     contains functions for constructing/optimizing/writing problem
#OUTPUT
#         peptides that are part of the solution are passed to STDOUT. 


source("../shared-scripts/trainset-sampling/problem-functions.R")

#self set and r value from input arguments
argv <- commandArgs( trailingOnly=TRUE ) 
selffile <- argv[1]             #the total set of self to optimize for
rvalue <- as.numeric(argv[2])   #the value of r (or k)
Npeps <- as.numeric(argv[3])    #the number of peptides the solution should (maximally) contain

if(Npeps > 0){
  #Get greedy solution:
  no.text.output <- capture.output({gt <- greedyTolerance(selffile,rvalue,Npeps)})
  
  sol <- gt$s
  steps <- gt$sv

  #convert back to peptides:
  all.peptides <- as.character(read.table(selffile)$V1) 
  sol.peptides <- all.peptides[sol]
  

  #order based on the step they were included.
  steps <- steps[steps>0]
  step.order <- rank(steps, ties.method="average")
  sol.peptides <- sol.peptides[order(step.order)]

  #output to stdout
  cat(sol.peptides,sep="\n")
} else {
  cat("")
}

