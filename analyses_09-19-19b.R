
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
                                 which(originateBeforeEnd & extinctAfterEnd & originateAfterStart)
                              }
)

# control box for ss analyses
#
# number of reps
nRep <- 20
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

traitScale <- scale(data[,includedTraits])

for(i in 1:length(withinInterval_list)){
   traitsInInterval <- traitScale[withinInterval_list[[i]],]
   validTaxon <- apply(traitsInInterval, 1, function(x) all(!is.na(x)))
   traitsValid <- traitsInInterval[validTaxon,]
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
      traitsBeforeInterval <- traitScale[withinInterval_list[[i-1]],]
      validTaxonBefore <- apply(traitsBeforeInterval, 1, function(x) all(!is.na(x)))
      traitsValidBefore <- traitsBeforeInterval[validTaxonBefore,]
      #
      traitsOrig <- traitScale[orig_withinInterval_list[[i]],,drop=FALSE]
      validOrig <- apply(traitsOrig, 1, function(x) all(!is.na(x)))
      traitsValidOrig <- traitsOrig[validOrig,,drop=FALSE]
      #
      if(nrow(traitsValidOrig) > 0){
         traitsValidBefore_And_Orig <- rbind(traitsValidOrig, traitsValidBefore)
         dist_orig<-as.matrix(dist(traitsValidBefore_And_Orig))
         diag(dist_orig) <- NA
         dist_orig <- dist_orig[1:nrow(traitsValidOrig),,drop=FALSE]
         min_dist_orig <- apply(dist_orig, 1, min, na.rm = TRUE)
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
stderrSS_m_pw_dist <- sapply(m_pw_dist_ssList, sd) / sqrt(nRep)
u95_SS_m_pw_dist <- meanSS_m_pw_dist + (1.96 * stderrSS_m_pw_dist)
l05_SS_m_pw_dist <- meanSS_m_pw_dist - (1.96 * stderrSS_m_pw_dist)
max_SS_pw_dist <- sapply(m_pw_dist_ssList,max)
min_SS_pw_dist <- sapply(m_pw_dist_ssList,min)
ylims <- c(min(c(min_SS_pw_dist,l05_SS_m_pw_dist)),
           max(c(max_SS_pw_dist,u95_SS_m_pw_dist)))
#
plot(midDate, m_pw_dist, 
     type="l", ylim=ylims, lwd =2, col="red")
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
     type="l", ylim=ylims, col="red",lwd=2)
abline(v=66, lty=3, lwd = 2)
lines(midDate,meanSS_m_NND, lwd=3)
lines(midDate, u95_SS_m_NND, lty=3, lwd=1)
lines(midDate, l05_SS_m_NND, lty=3, lwd=1)
lines(midDate, max_SS_m_NND, 
      col="blue",lty=4)
lines(midDate, min_SS_m_NND, 
      col="blue",lty=4)

#################################


plot(midDate,innovation, type="l")
lines(loess(innovation~midDate))
abline(v=66,lty=3,lwd=2)


################################

survData <- read.table("workspace/Test-Complexity/SurviorData.txt", header = TRUE)

# use 'inc' from Fraas's rmarkdown to pick the chars they used
refCentroid <- survData[,inc]
refCentroid <- apply(refCentroid, 2, mean)

# calculate for each species
diff_ref <- abs(t(apply(data[,inc], 1, function(x) x - refCentroid)))
# rescale to 0 - 1
diff_ref <- apply(diff_ref, 2, function(x) x/max(x, na.rm = TRUE) )
diff_ref <- apply(diff_ref, 2, function(x) x - min(x, na.rm = TRUE) )

boxplot(diff_ref)

# get distance for each species by summing across all included traits
dist_ref <- apply(diff_ref, 1, sum)

# per-interval 
mean_dist_ref_5surv <- numeric()
mean_dist_ref_ssList <- list()

nRep <-100

for(i in 1:length(withinInterval_list)){
   dist_ref_inInterval <- dist_ref[withinInterval_list[[i]]]
   dist_ref_inInterval <- dist_ref_inInterval[!is.na(dist_ref_inInterval)]
   mean_dist_ref_5surv[i] <- mean(dist_ref_inInterval)
   #
   mean_dist_ref_ss<-numeric()
   #
   for(rep in 1:nRep){
      pickedTaxa <- sample(1:length(dist_ref_inInterval),
                           rareSS, replace=FALSE)
      dist_ref_picked <- dist_ref_inInterval[pickedTaxa]
      mean_dist_ref_ss[rep] <- mean(dist_ref_picked)   
   }
   mean_dist_ref_ssList[[i]] <- mean_dist_ref_ss
}

#plot(midDate, mean_dist_ref_5surv, type="l",xlim=c(80,50))
#abline(v=66,col="red",lty=3,lwd=2)   

# sample size standardized
meanSS_mean_dist_ref <- sapply(mean_dist_ref_ssList,mean)
plot(midDate, meanSS_mean_dist_ref, type="l",xlim=c(80,50))
lines(midDate,mean_dist_ref_5surv, lty=2, col="blue", lwd=2)
abline(v=66,col="red",lty=3,lwd=2)   

############

simpleTaxa <- sapply(data$species, 
       function(x) x == "holmdelensis" | x =="monmouthensis")
simpleTaxa <- data[simpleTaxa,]
# use 'inc' from Fraas's rmarkdown to pick the chars they used
refCentroid <- simpleTaxa[,inc]
refCentroid <- apply(refCentroid, 2, mean)

# calculate for each species
diff_ref <- abs(t(apply(data[,inc], 1, function(x) x - refCentroid)))
# rescale to 0 - 1
diff_ref <- apply(diff_ref, 2, function(x) x/max(x, na.rm = TRUE) )
diff_ref <- apply(diff_ref, 2, function(x) x - min(x, na.rm = TRUE) )

boxplot(diff_ref)

# get distance for each species by summing across all included traits
dist_ref <- apply(diff_ref, 1, sum)

# per-interval 
mean_dist_ref_2simple <- numeric()
mean_dist_ref_ssList <- list()

nRep <-100

for(i in 1:length(withinInterval_list)){
   dist_ref_inInterval <- dist_ref[withinInterval_list[[i]]]
   dist_ref_inInterval <- dist_ref_inInterval[!is.na(dist_ref_inInterval)]
   mean_dist_ref_2simple[i] <- mean(dist_ref_inInterval)
   #
   mean_dist_ref_ss<-numeric()
   #
   for(rep in 1:nRep){
      pickedTaxa <- sample(1:length(dist_ref_inInterval),
                           rareSS, replace=FALSE)
      dist_ref_picked <- dist_ref_inInterval[pickedTaxa]
      mean_dist_ref_ss[rep] <- mean(dist_ref_picked)   
      }
   mean_dist_ref_ssList[[i]] <- mean_dist_ref_ss
   }

#plot(midDate, mean_dist_ref, type="l",xlim=c(80,50))
#abline(v=66,col="red",lty=3,lwd=2)   

# sample size standardized
meanSS_mean_dist_ref <- sapply(mean_dist_ref_ssList,mean)

ylims <- c(meanSS_mean_dist_ref,
         mean_dist_ref_2simple,
         mean_dist_ref_5surv)
ylims <- c(min(ylims),max(ylims))
plot(midDate, meanSS_mean_dist_ref, type="l",
     xlim=c(80,50), ylim=ylims)
lines(midDate,mean_dist_ref_2simple, lty=2, col="blue", lwd=2)
lines(midDate,mean_dist_ref_5surv, lty=2, col="orange", lwd=3)
abline(v=66,col="red",lty=3,lwd=2)   

