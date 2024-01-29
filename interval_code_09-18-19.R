    
    
    
    FAD <- 
    LAD <- 
    
    # 0.25 Mya bins
    time <- seq(100,0,by=-0.25)
    
    withinInterval_list <- lapply(2:length(time),
           function(t) {
             int_end <- time[t] 
             int_start <- time[t-1]
             originateBeforeEnd <- FAD > int_end
             extinctAfterEnd <- LAD < int_start
             which(originateBeforeEnd & extinctAfterEnd)
           }
           )
    
    
    # and then to get the traits (for example) of taxa in each interval...
    
    for(i in 1:length(withinInterval_list)){
      traitsInInterval <- traits[withinInterval_list[[i]],]
      # something
      # something
    }
    
    
    
    