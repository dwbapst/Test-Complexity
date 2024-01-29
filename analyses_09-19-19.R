
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
   
#############################################

data <- read.csv(file="workspace/Test-Complexity/MasterList (vstat 1.05).csv")

###############


### This is the full list of characters for the analysis (inc) and the list for the pca (inc2)
inc <- c(#"w", 
   #"mtheta", 
   #"lw", 
   #"lh", 
   #"ic1", 
   #"ic2", 
   #"ic3", 
   "fcirc", 
   "t", 
   "numcham", 
   "expans", 
   #"height", 
   #"length", 
   #"fcangle", 
   "area", 
   "fcarea", 
   "clava", 
   #"chamwl", 
   "keel", 
   #"bidors", 
   #"biven", 
   "biconvex", 
   "lobe"
   #"double"
   #"depth"
)

inc2 <- c(#"w", 
   "mtheta", 
   "lw", 
   "lh", 
   #"ic1", 
   #"ic2", 
   #"ic3", 
   "fcirc", 
   "t", 
   "numcham", 
   "expans", 
   #"height", 
   #"length", 
   "fcangle", 
   "area", 
   #"fcarea", 
   "clava", 
   "chamwl", 
   "keel", 
   #"bidors", 
   #"biven", 
   "biconvex", 
   "lobe", 
   "double"
   #"depth"
)

includeAll <- c(
   "w", 
   "mtheta", 
   "lw", 
   "lh", 
   "ic1", 
   "ic2", 
   "ic3", 
   "fcirc", 
   "t", 
   "numcham", 
   "expans", 
   "height", 
   "length", 
   "fcangle", 
   "area", 
   "fcarea", 
   "clava", 
   "chamwl", 
   "keel", 
   "bidors", 
   "biven", 
   "biconvex", 
   "lobe", 
   "double",
   "depth"
   )

includedTraits <- includeAll

############

FAD <- data$origin
LAD <- data$extin

# 0.25 Mya bins
time <- seq(100,0,by=-0.25)
midDate <- sapply(2:length(time), function(i) mean(time[i],time[i-1]))

withinInterval_list <- lapply(2:length(time),
                              function(t) {
                                 int_end <- time[t]
                                 int_start <- time[t-1]
                                 originateBeforeEnd <- FAD > int_end
                                 extinctAfterEnd <- LAD < int_start
                                 which(originateBeforeEnd & extinctAfterEnd)
                              }
                           )

orig_withinInterval_list <- lapply(2:length(time),
                              function(t) {
                                 int_end <- time[t]
                                 int_start <- time[t-1]
                                 originateBeforeEnd <- FAD > int_end
                                 originateAfterStart <- FAD < int_start                                 
                                 extinctAfterEnd <- LAD < int_start
                                 which(originateBeforeEnd & extinctAfterEnd)
                              }
)

# control box for ss analyses
#
# number of reps
nRep <- 10000
# min sample size to rarefy to
rareSS <- 15

# valid taxa, m pw dist, nnd
m_pw_dist <- numeric()
m_NND <- numeric()
nValidTaxon <- numeric()
#
centroidList <- list()
#
innovation<-numeric()
#
# sub sample values
m_pw_dist_ssList <- list()
m_NND_ssList <- list()

for(i in 1:length(withinInterval_list)){
   traitsInInterval <- data[withinInterval_list[[i]],includedTraits]
   validTaxon <- apply(traitsInInterval, 1, function(x) all(!is.na(x)))
   traitsValid <- scale(traitsInInterval[validTaxon,])
   nValidTaxon[i] <- sum(validTaxon)
   #
   # calculate centroid morph for z-scaled trait data
   centroidList[[i]] <- apply(traitsValid,2,mean)
   ##
   # non-standardized disparity measures
   distmatrix <- dist(traitsValid, diag = FALSE)
   m_pw_dist[i]<- mean(distmatrix)
   m_NND[i] <- mean(nearestNeighborDist(distmatrix)) 
   #
   ###################################
   # calculate innovation
   if(i>1 & length(orig_withinInterval_list[[i]])>0){
      traitsBeforeInterval <- data[withinInterval_list[[i-1]],includedTraits]
      validTaxonBefore <- apply(traitsBeforeInterval, 1, function(x) all(!is.na(x)))
      traitsValidBefore <- scale(traitsBeforeInterval[validTaxonBefore,])
      #
      traitsOrig <- data[orig_withinInterval_list[[i]],includedTraits]
      validOrig <- apply(traitsOrig, 1, function(x) all(!is.na(x)))
      traitsValidOrig <- scale(traitsOrig[validOrig,,drop=FALSE])
      #
      if(nrow(traitsValidOrig) > 0){
         traitsValidBefore_And_Orig <- rbind(traitsValidOrig, traitsValidBefore)
         dist_orig<-as.matrix(dist(traitsValidBefore_And_Orig))
         diag(dist_orig) <- NA
         dist_orig <- dist_orig[1:nrow(traitsValidOrig),,drop=FALSE]
         min_dist_orig <- apply(dist, 1, min, na.rm=FALSE)
         innovation[i] <- mean(min_dist_orig)
      }else{
         innovation[i]<-NA
         } 
   }else{
      innovation[i]<-NA
      }   
   ########################################
   # subsample to rareSS, over nRep
   m_pw_dist_ss <- numeric()
   m_NND_ss <- numeric()
   for(rep in 1:nRep){
      pickedTaxa <- sample(1:nrow(traitsValid), size=rareSS, replace = FALSE)
      traitsPicked <- traitsValid[pickedTaxa,]
      distmatrix <- dist(traitsPicked, diag = FALSE)
      #
      m_pw_dist_ss[rep]<- mean(distmatrix)
      m_NND_ss[rep] <- mean(nearestNeighborDist(distmatrix))
      }
   m_pw_dist_ssList[[i]] <- m_pw_dist_ss
   m_NND_ssList[[i]] <- m_NND_ss
   }

meanSS_m_pw_dist <- sapply(m_pw_dist_ssList, mean)
u95_SS_m_pw_dist <- sapply(m_pw_dist_ssList, quantile, 0.75)
l05_SS_m_pw_dist <- sapply(m_pw_dist_ssList, quantile, 0.25)
ylims <- c(min(l05_SS_m_pw_dist), max(u95_SS_m_pw_dist))

plot(midDate,meanSS_m_pw_dist,type="l",ylim=ylims)
abline(v=66, lty=3, lwd = 2)
lines(midDate, u95_SS_m_pw_dist, lty=3, lwd=1)
lines(midDate, l05_SS_m_pw_dist, lty=3, lwd=1)

#

meanSS_m_NND <- sapply(m_NND_ssList, mean)
u95_SS_m_NND <- sapply(m_NND_ssList, quantile, 0.75)
l05_SS_m_NND <- sapply(m_NND_ssList, quantile, 0.25)
ylims <- c(min(l05_SS_m_NND), max(u95_SS_m_NND))

plot(midDate,meanSS_m_NND,type="l",ylim=ylims)
abline(v=66, lty=3, lwd = 2)
lines(midDate, u95_SS_m_NND, lty=3, lwd=1)
lines(midDate, l05_SS_m_NND, lty=3, lwd=1)

#


###########################################################

meanSS_m_pw_dist <- sapply(m_pw_dist_ssList, mean)
stderrSS_m_pw_dist <- sapply(m_pw_dist_ssList, sd) / sqrt(nRep)
u95_SS_m_pw_dist <- meanSS_m_pw_dist + (1.96 * stderrSS_m_pw_dist)
l05_SS_m_pw_dist <- meanSS_m_pw_dist - (1.96 * stderrSS_m_pw_dist)
max_SS_pw_dist <- sapply(m_pw_dist_ssList,max)
min_SS_pw_dist <- sapply(m_pw_dist_ssList,min)
ylims <- c(min(c(min_SS_pw_dist,l05_SS_m_pw_dist)),
           max(c(max_SS_pw_dist,u95_SS_m_pw_dist)))
#
plot(midDate, m_pw_dist, 
     type="l", ylim=ylims, col="red")
abline(v=66, lty=3, lwd = 2)
lines(midDate,meanSS_m_pw_dist, lwd=3)
lines(midDate, u95_SS_m_pw_dist, lty=3, lwd=1)
lines(midDate, l05_SS_m_pw_dist, lty=3, lwd=1)
lines(midDate, max_SS_pw_dist, col="blue",lty=4)
lines(midDate, min_SS_pw_dist, col="blue",lty=4)


meanSS_m_NND <- sapply(m_NND_ssList, mean)
min_SS_m_NND <- sapply(m_NND_ssList, min)
max_SS_m_NND <- sapply(m_NND_ssList, max)
stderrSS_m_NND <- sapply(m_NND_ssList, sd) / sqrt(nRep)
u95_SS_m_NND <- meanSS_m_NND + (1.96 * stderrSS_m_NND)
l05_SS_m_NND <- meanSS_m_NND - (1.96 * stderrSS_m_NND)
ylims <- c(min(c(min_SS_m_NND,l05_SS_m_NND)), 
           max(c(max_SS_m_NND,u95_SS_m_NND)))
#
plot(midDate, m_NND, 
     type="l", ylim=ylims, col="red")
abline(v=66, lty=3, lwd = 2)
lines(midDate,meanSS_m_NND, lwd=3)
lines(midDate, u95_SS_m_NND, lty=3, lwd=1)
lines(midDate, l05_SS_m_NND, lty=3, lwd=1)
lines(midDate, max_SS_m_NND, 
      col="blue",lty=4)
lines(midDate, min_SS_m_NND, 
      col="blue",lty=4)



