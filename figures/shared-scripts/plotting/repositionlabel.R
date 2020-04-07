#---------------------------DESCRIPTION
#REPOSITIONLABEL replaces values that are too close together.
#Input arguments:
#           - datavec: a vector with the values that need repositioning
#           - threshold: specifying how far values should at least be apart




#---------------------------CODE
repositionlabel <- function(datavec,threshold){
  
  #First sort the values in ascending order:
  data.order <- order(datavec)                        #Determines the order of the values for sorting
  sort.data <- datavec[data.order]                    #Subsetting with the data.order is the same as sorting
  inv.data.order <- order(order(datavec))             #The inverse of data.order, so 'sorting' the sorted data 
                                                      #with this vector will return original data. 
  
  #For each point in the sorted data: calculate distances to left and right neighbor, respectively:
  dist.left.neighbor <- c(Inf,diff(sort.data))
  dist.right.neighbor <- c(diff(sort.data),Inf)
  
  #Identify stretches of positions that are too close (distance below the threshold):
  
  stretchvec <- rep(0,length(datavec))                #Create a vector with zeros of the same length as datavec.

  stretch.index <- 1
  for (i in 1:length(stretchvec)){
    if((dist.left.neighbor[i]<threshold) | (dist.right.neighbor[i]<threshold)){ #Part of a stretch if at least one neighbor is too close.
      stretchvec[i]<- stretch.index                   #Stretchvec gets value != 0
      if(dist.right.neighbor[i]>=threshold){          #so at the final position of each single stretch: 
        stretch.index <- stretch.index + 1            #increase stretch.index with 1.
      }
    }
  }

  
  S <- rle(stretchvec)                                #Finds the stretches and returns a vector with stretch lengths and stretch values
  stretch.num <- length(S$lengths)                    #Determines the number of stretches
  
  #Reassign the data for the positions where at least one neighbor is too close (so code 1 in stretchvec)
  new.values <- sort.data                             #First create a new vector, start with the unchanged sorted data values
  
  for (i in seq(1,stretch.num)){                      #Loop over the stretches
    if (S$values[i]!=0){                              #Only do something for the positions with a nonzero code
      
      #Get length of the current stretch and its position indices within the sorted data:
      stretch.length <- S$lengths[i]
      stretch.start <- sum(S$lengths[1:i]) - S$lengths[i] + 1
      stretch.end <- stretch.start + stretch.length - 1

      #Calculate the new values: each a threshold apart and surrounding the mean value of the stretch
      new.stretch.values <- (seq(1,stretch.length) - mean(seq(1,stretch.length)))*threshold + mean(sort.data[stretch.start:stretch.end])
      
      #Insert the new values into the total new.values vector at the correct position:
      new.values[stretch.start:stretch.end] <- new.stretch.values
    }
  }
  
  #Now use the previously defined inv.data.order to 'unsort' the new values back into the original order of datavec
  unsort.data <- new.values[inv.data.order]
  
  
  #Return this unsorted vector with new values as function output
  unsort.data 
  
}

