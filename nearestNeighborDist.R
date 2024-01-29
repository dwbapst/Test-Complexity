
nearestNeighborDist <- function(distMat){
   #returns a vector of NNDs for each taxon in a distance matrix
   #this function is included in paleotree mainly for pedagogical use
   if(!inherits(distMat,"dist")){
      if(is.matrix(distMat)){
         if(any(sapply(diag(distMat),function(x) !is.na(x) & x!=0))){
            stop("Diagonal is nonzero and not NA, may be a similarity matrix, not a distance matrix")
            }
         if(!isSymmetric(distMat)){
            stop("Not a Symmetric Distance Matrix?")
            }
         distMat <- distMat
      }else{
         stop("Not a matrix, not a 'dist' object, what is it?")
         }
   }else{
      distMat <- as.matrix(distMat)
      }
   NND <- sapply(1:nrow(distMat),function(x) min(distMat[x,-x]))
   names(NND) <- labels(distMat)[[1]]
   return(NND) #return as a per-taxon 
   }
   

data <- read.table(file="")
