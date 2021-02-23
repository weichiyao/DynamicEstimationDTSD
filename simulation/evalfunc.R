#### Update July 4th
############################################################################################
evalfunc <- function(h, hhat, method = c("alor", "cindex", "adist")){
  ###########################################################
  # GOAL # 
  # Evaluate the accuracy of the estimation
  
  # INPUT #
  # h: true hazard rates
  # hhat: estimated hazard rates
  # method: "alor" -- the absolute log odds ratio (ALOR)
  #         "cindex -- estimated concordance probabilities
  
  # OUTPUT #
  # Values of the corresponding evaluation measures
  ###########################################################
  
  nh <- length(h)
  if (nh != length(hhat)) stop("length of h and length of hhat should be the same!")
  
  if (method == "alor") {
    hhat[hhat == 0] <- 0.01
    hhat[hhat == 1] <- 0.99
    
    h[h == 0] <- 0.01
    h[h == 1] <- 0.99
    ret <- mean(abs(log((hhat*(1 - h))/((1 - hhat)*h))))
  } else if (method == "cindex"){
    oh <- order(h, decreasing = TRUE)
    ho <- h[oh]
    hhato <- hhat[oh]
    
    # After ordering, number of pairs that are with I{h_i > h_j} = 1 => denominator
    cpd <- nh * (nh - 1) / 2
    
    # hhat is ordered following h, therefore, count the number of the right ordering 
    cpn <- sapply(1:(nh - 1), function(hi) {
      sum(hhato[hi] > hhato[(hi + 1):nh])
    })

    ret <- sum(cpn) / cpd
  } else if (method == "adist"){
    ret = mean(abs(hhat - h))
  } else {
    stop("Wrong evaluation method is given!")
  }
  
  ret
}
