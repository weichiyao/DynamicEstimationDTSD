#### Updated July 11th: censoring rate 0
#####################======== Large number of pseudo-subjects ===============#####################
traindtv_autocorr_gnrt <- function(nsub = 200, 
                                   model = c("linear", "nonlinear", "interaction"), 
                                   distribution = c("Exp", "Gtz", "WI"), 
                                   nperiod = c(4, 8),
                                   matsigma = NULL,
                                   censor.rate = c("10%", "50%"), 
                                   SNR = c("low", "high"),
                                   scenario = c("0TI2TV", "0TI4TV", "1TI4TV", "2TI1TV", "2TI4TV")){
  
  nstime <- 1000
  
  RET <- tvstimegnrt(nsub = nstime,
                     scenario = scenario,
                     model = model,
                     distribution = distribution,
                     nperiod = nperiod,
                     matsigma = matsigma,
                     SNR = SNR)
  Coeff <- RET$coeff
  chngpt <- findsurvint(y = sort(RET$survtime),
                        nper = nperiod, 
                        rate = 0.10)
  rm(RET)
  gc()
  
  Data <- as.data.frame(matrix(NA, nperiod * nsub, 13))
  names(Data) <- c("ID", "X1", "X2", "X3", "X4", "X5", "X6",
                   "StopT", "Time", "Event", "Theta", "S", "h")
  
  Data$Time <- rep(1:nperiod, nsub)
  Data$ID <- rep(1:nsub, each = nperiod)
  Data[, 2:7] <- genvar(nsub = nsub, 
                        nperiod = nperiod, 
                        matsigma = matsigma,
                        scenario = scenario)
  Data$StopT <- rep(chngpt, nsub)
  
  Data$Theta <- create_theta(model = model, 
                             data = as.matrix(Data[, 2:7]), 
                             coeff = Coeff, 
                             SNR = SNR, 
                             scenario = scenario)
  
  Hfunc <- switch(distribution,
                  "Exp" = ExpHfunc,
                  "WI" = WHfunc,
                  "Gtz" = GtzHfunc
  )
  tfunc <- switch(distribution,
                  "Exp" = Exptfunc,
                  "WI" = Wtfunc,
                  "Gtz" = Gtztfunc
  )
  
  TS <- as.vector(rep(c(0, chngpt), nsub))
  
  tlen <- length(TS)
  seqt2 <- (nperiod + 1) * c(1:(tlen / (nperiod + 1))) # each : -nperiod
  seqt1 <- (nperiod + 1) * c(0:((tlen - 1) / (nperiod + 1))) + 1 # each : -1
  
  # Temp cumulative hazards at each change point
  R <- Hfunc(ts1 = TS[-seqt1], 
             ts2 = TS[-seqt2], 
             theta = Data$Theta, 
             coeff = Coeff)
  TS <- TS[-seqt2]
  rm(seqt2)
  rm(seqt1)
  # Temp cumulative hazards at each change point: Each row belongs to a subject
  R <- matrix(R, ncol = nperiod, byrow = TRUE)
  
  # Cumulative hazards at each change point: Each column belongs to a subject
  R <- apply(R, 1, cumsum)
  
  # Survprob: each column belongs to a subject
  survprob <- exp(-R)
  Data$h <- as.vector(sapply(1:nsub, function(ni) (c(1, survprob[, ni][-nperiod]) - survprob[, ni]) / c(1, survprob[, ni][-nperiod])))
  Data$S <- as.vector(survprob)
  rm(survprob)
  
  U <- runif(nsub)
  
  TS <- as.vector(rep(c(0, chngpt[1:(nperiod - 1)]), nsub))
  for (Count in 1:nsub) {
    idxC <- which(Data$ID == Count)
    VEC <- c(0, R[, Count], Inf)
    rID <- findInterval(-log(U[Count]), VEC)
    if (rID == 1){
      Data[idxC, ][1, ]$Event <- 1
    } else if (rID <= nperiod){
      Data[idxC, ][1:(rID - 1), ]$Event <- 0
      Data[idxC, ][rID, ]$Event <- 1
    } else {
      Data[idxC, ]$Event <- 0
    }
  }
  
  Data$X7 <- rep(runif(nsub, 0, 1), each = nperiod)
  Data <- Data[!is.na(Data$Event), ]
  rm(R)
  gc()
  
  if (length(unique(Data$ID)) != nsub){
    stop("ID length NOT equal to nsub")
  }
  
  nobs <- nrow(Data)
  Data$X8 <- sample(c(0, 1), nobs, replace = TRUE)
  Data$X9 <- sample(1:5, nobs, replace = TRUE)
  Data$X10 <- runif(nobs, 0, 2)
  
  
  Censor.time <- create_ctime(nsub = nsub, 
                              SNR = SNR, 
                              model = model, 
                              distribution = distribution,
                              censor.rate = censor.rate,
                              nperiod = nperiod,
                              scenario = scenario)
  ###================== Add Censoring =========================================
  for (j in 1:nsub ){
    idxj <- which(Data$ID == j)
    Vec <- c(0, Data[idxj, ]$StopT, Inf)
    ID <- findInterval(Censor.time[j], Vec)
    
    if( ID <= nrow(Data[idxj, ]) ){
      Data[idxj, ][ID, ]$Event <- 0
      Data[idxj, ][ID, ]$StopT <- Censor.time[j]
      
      TALLjj <- c(0, chngpt)[ID:(ID + 1)]
      TALLjj[2] <- Censor.time[j]
      
      Sjj <- c(1, Data[idxj, ]$S)
      
      Hjj <- Hfunc(ts1 = TALLjj[2], ts2 = TALLjj[1], theta = Data[idxj, ][ID, ]$Theta, coeff = Coeff)
      Data[idxj, ][ID, ]$S <- Sjj[ID] * exp(-Hjj)
      Data[idxj, ][ID, ]$h <- (Sjj[ID] - Data[idxj, ][ID, ]$S) / Sjj[ID]
      if( ID != length(idxj) ){
        Data[idxj, ][(ID + 1):(length(idxj)), ]$Event <- NA
      }
    } 
  }
  Data <- Data[!is.na(Data$Event), ]
  if (length(unique(Data$ID)) != nsub){
    stop("ID length NOT equal to nsub")
  }
  
  
  Data$I <- 1:nrow(Data)
  Data$Theta <- NULL
  RET <- NULL
  RET$fullData <- Data
  RET$Info = list(Model = model, 
                  DRate = sum(Data$Event) / nsub, 
                  Dist = distribution, 
                  Coeff = Coeff)
  
  rm(Data)
  gc()
  return(RET)
}