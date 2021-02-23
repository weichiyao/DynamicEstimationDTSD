#####################################################################################################
create_theta <- function(data, scenario, model, coeff, SNR){
  a <- 1.5
  if (SNR == "high") a <- 3
  if (model == "linear"){
    Fstar <- data %*% coeff$Beta1 + coeff$Beta0[1]
  } else if (model == "nonlinear"){
    Fstar <- a * log(abs(data %*% coeff$Beta2 + coeff$Beta0[2])) + coeff$Beta0[1]
  } else if (model == "interaction"){
    Fstar <- itct_term(data = data, scenario = scenario, coeff = coeff, SNR = SNR)
  } else {
    stop("Wrong model is set.")
  }
  return(exp(Fstar))
}
## Compute the cumulative hazards function
ExpHfunc <- function(ts1, ts2, theta, coeff){
  coeff$Lambda * theta * (ts1 - ts2)
}

WHfunc <- function(ts1, ts2, theta, coeff){
  coeff$Lambda * theta * (ts1 ^ coeff$V - ts2 ^ coeff$V)
}

GtzHfunc <- function(ts1, ts2, theta, coeff){
  coeff$Lambda * theta / coeff$Alpha * 
    (exp(coeff$Alpha * ts1) - exp(coeff$Alpha * ts2))
}


## Compute the continuous times
Exptfunc <- function(tall, theta, coeff, t0, rid){
  t0 / coeff$Lambda / theta[rid] + tall[rid]
}

Wtfunc <- function(tall, theta, coeff, t0, rid){
  (t0 / coeff$Lambda / theta[rid] + tall[rid] ^ coeff$V) ^ (1 / coeff$V)
}

Gtztfunc <- function(tall, theta, coeff, t0, rid){
  t0 <- coeff$Alpha * t0 / coeff$Lambda / theta[rid]
  log(t0 + exp(coeff$Alpha * tall[rid])) / coeff$Alpha
}


findsurvint <- function(y, nper, rate) {
  # INPUT # 
  # y = true survival times (not censored) from the existing DGP used in the previous paper
  # nper = the number of intervals / periods we want
  # rate = proportion of censoring at the end that we want 
  #        (we will also add censoring at any moment later with the DGP itself)
  
  # OUTPUT # 
  # The interval limits to define the periods
  
  int <- quantile(y, probs = seq((1 - rate) / nper, 1 - rate, length.out = nper))
  
  int
  
}

#####################======== MAIN FUNCTION ===============#####################
tvstimegnrt <- function(nsub = 200, 
                        scenario = c("0TI2TV", "0TI4TV", "1TI4TV", "2TI1TV", "2TI4TV"), 
                        model = c("linear", "nonlinear", "interaction"), 
                        distribution = c("Exp", "Gtz", "WI"), 
                        nperiod = c(4, 8),
                        matsigma = NULL, 
                        SNR = c("high", "low")){
  
  Data <- matrix(NA, nperiod * nsub, 7)
  colnames(Data) <- c("ID","X1","X2","X3","X4","X5","X6")
  Data[, 1] <- rep(1:nsub, each = nperiod)
  Data[, 2:7] <- genvar(nsub = nsub, 
                        nperiod = nperiod, 
                        matsigma = matsigma,
                        scenario = scenario)
  
  ## Set the coefficients and compute the Theta = exp(f(X))
  coeffTS <- create_coeff(scenario = scenario, 
                          model = model, 
                          distribution = distribution, 
                          SNR = SNR, 
                          nsub = nsub,
                          nperiod = nperiod)
  Coeff <- coeffTS$Coeff
  TS <- coeffTS$TS
  rm(coeffTS)
  
  Theta <- create_theta(model = model, 
                        data = Data[, 2:7], 
                        coeff = Coeff,
                        scenario = scenario,
                        SNR = SNR)

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
  
  tlen <- length(TS)
  seqt2 <- nperiod * c(1:(tlen / nperiod)) # each : -nperiod
  seqt1 <- nperiod * c(0:((tlen - 1) / nperiod)) + 1 # each : -1
  
  R <- Hfunc(ts1 = TS[-seqt1], 
             ts2 = TS[-seqt2], 
             theta = Theta[-seqt2], 
             coeff = Coeff)
  # each row belongs to a subject
  R <- matrix(R, ncol = nperiod - 1, byrow = TRUE)
  
  U <- runif(nsub)
  survtime <- rep(0, nsub)
  survnrow <- rep(0, nsub)
  for (Count in 1:nsub) {
    idxC <- which(Data[, "ID"] == Count)
    VEC <- c(0, cumsum(R[Count, ]), Inf)
    R.ID <- findInterval(-log(U[Count]), VEC)
    TT <- -log(U[Count]) - VEC[R.ID] 
    survnrow[Count] <- R.ID
    survtime[Count] <- tfunc(tall = TS[idxC], 
                             theta = Theta[idxC], 
                             coeff = Coeff, 
                             t0 = TT, 
                             rid = R.ID)
  }
  rm(R)
  rm(Theta)
  rm(U)
  gc()
  RET = list(survtime = survtime,
             coeff = Coeff)
  return(RET)
}

