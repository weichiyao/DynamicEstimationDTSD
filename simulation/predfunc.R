predfunc <- function(scenario = c("0TI2TV", "0TI4TV", "1TI4TV", "2TI1TV", "2TI4TV"),
                     SNR = c("low", "high"), 
                     num.train = c(200, 1000, 5000),
                     num.test = 1000,
                     censor.rate = c("10%", "50%"),
                     autocorrtype = c("weak", "strong"),
                     distribution = c("Exp", "WI", "Gtz"),
                     relationship = c("linear", "nonlinear", "interaction"),
                     measurement = c("adist", "alor", "cindex"),
                     num.period = c(4, 8)){
  
  ################################################################################
  ## == GOAL:
  ## This function is produce the measurment values 
  ##
  ## == INPUT:
  ## scenario     -- different TI/TV combinations: c("0TI2TV", "0TI4TV", "1TI4TV", "2TI1TV", "2TI4TV")
  ## SNR          -- signal to noise ratio: c("low", "high")
  ## censor.rate  -- censoring rate: c("10%", "50%"),
  ## autocorrtype -- autocorrelation type: c("weak", "strong"),
  ## distribution -- survival distribution: c("Exp", "WI", "Gtz"),
  ## relationship -- survival relationship: c("linear", "nonlinear", "interaction"),
  ## measurement  -- evaluation measure: c("adist", "alor", "cindex")
  ## num.period   -- number of periods: c(4, 8) 
  ## nsim         -- number of simulations
  ## num.train    -- number of training samples
  ## num.test     -- number of test samples
  ##
  ## == OUTPUT:
  ## Measurement values for six methods: benchmark (bench), Separate with Hellinger splitting rule (sephlgr)
  ## SKHT with Hellinger splitting rule (skhthlgr), Superpp with Hellinger splitting rule (supphlgr),
  ## Superpp with Hellinger splitting rule based on baseline values (supp0hlgr),
  ## Superpp with DTPO (suppDTPO).
  ################################################################################
  
  cscen <- match(scenario, c("0TI2TV", "0TI4TV", "1TI4TV", "2TI1TV", "2TI4TV"), nomatch = 0)
  csignal <- match(SNR, c("low", "high"), nomatch = 0)
  ccrate <- match(censor.rate, c("10%", "50%"), nomatch = 0)
  cautocorr <- match(autocorrtype, c("weak", "strong"), nomatch = 0)
  cmodel <- match(relationship, c("linear", "nonlinear", "interaction"), nomatch = 0)
  cmeas <- match(measurement, c("adist", "alor", "cindex"), nomatch = 0)
  cnperiod <- match(num.period, c(4, 8), nomatch = 0)
  cdist <- match(distribution, c("Exp", "WI", "Gtz"), nomatch = 0)
  if (cscen == 0) stop("Wrong scenario is given! ") 
  if (csignal == 0) stop("Wrong SNR is given! ")
  # if (cN == 0) stop("Wrong num.train is given! ")
  if (cautocorr == 0) stop("Wrong autocorrtype is given! ")
  if (ccrate == 0) stop("Wrong censor.rate is given! ")
  if (cmodel == 0) stop("Wrong relationship is given! ")
  if (cdist == 0) stop("Wrong distribution is given! ")
  if (cmeas == 0) stop("Wrong measurement is given! ")
  if (cnperiod == 0) stop("num.period can only set to be either 4 or 8.")
  
  scen <- c("0TI2TV", "0TI4TV", "1TI4TV", "2TI1TV", "2TI4TV")[cscen]
  signal <- c("low", "high")[csignal]
  crate <- c("10%", "50%")[ccrate]
  autocorr <- c("weak", "strong")[cautocorr]
  model <- c("linear", "nonlinear", "interaction")[cmodel]
  dist <- c("Exp", "WI", "Gtz")[cdist]
  meas <- c("adist", "alor", "cindex")[cmeas]
  nperiod <- c(4, 8)[cnperiod]
  ntrain <- num.train
  ntest <- num.test
  
  if (nperiod == 8) {
    if (scen != "2TI4TV") stop("Plots for num.period = 8 only show for scenario '2TI4TV'")
  }
  
  testN <- ntest * 20
  ## Generate the autocorrelation matrix
  matsigma <- create_matsigma(scenario = scen, autocorrtype = autocorr)
  ##### ==== Generate the testdata ==== #####
  testData <- traindtv_autocorr_gnrt(nsub = testN, 
                                     model = model, 
                                     distribution = dist, 
                                     censor.rate = "10%", 
                                     nperiod = nperiod,
                                     matsigma = matsigma,
                                     scenario = scen,
                                     SNR = signal)$fullData
  testData <- testData[testData$Time <= nperiod, ]
  testData <- testdtv_gnrt(data = testData, ntest = ntest, id = "ID", period = "Time", y = "Event")
  
  ## Generate the training data
  RET <- traindtv_autocorr_gnrt(nsub = ntrain, 
                                model = model, 
                                distribution = dist, 
                                nperiod = nperiod,
                                matsigma = matsigma,
                                censor.rate = crate,
                                scenario = scen,
                                SNR = signal)
  trainData <- RET$fullData
  trainData <- trainData[trainData$Time <= nperiod, ]
  rm(RET)
  gc()
  
  ## Set the formula
  CFsep <- Event ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10
  CFsupp <- Event ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + current + Time 
  CFSKHT <- Event ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + Time
  timefixedcov <- c("X1","X2","X7")
  timevaryingcov = c("X3","X4","X5","X6","X8","X9","X10")
  ntree <- 500L
  
  ##### ==== Create the test dataset for each pair of (t, u) ==== #####
  testDatatu <- rep(list(0), nperiod * (nperiod + 1)/2)
  kk <- 1
  h <- rep(list(0), nperiod * (nperiod + 1)/2)
  utlist <- data.frame(matrix(0, nrow = nperiod * (nperiod + 1)/2, ncol = 2))
  names(utlist) = c("t", "u")
  
  for (currentperiod in 0:(nperiod - 1)){
    for (predperiod in (currentperiod + 1):nperiod){
      # For each pair of (t, u)
      testDatatu[[kk]] <- createseparate(data = testData[[predperiod]], id = "ID", 
                                         timefixedcov = timefixedcov, 
                                         timevaryingcov = timevaryingcov, 
                                         y = "Event", period = "Time", hazard = "h",
                                         currentperiod = currentperiod,
                                         predperiod = predperiod)
      h[[kk]] <- testDatatu[[kk]]$h
      utlist$t[kk] <- currentperiod
      utlist$u[kk] <- predperiod
      kk <- kk + 1
    }
  }
  rm(kk)
  rm(testData)
  gc()
  
  ################# == 1. Current benchmark with no covariates == #################
  # The ratio of the number of events to the total number of subjects at risk of experiencing the event.
  
  # npj -- u
  #         (t, u)
  # u = 1   (0, 1)
  # u = 2   (0, 2)  (1, 2)
  # u = 3   (0, 3)  (1, 3)  (2, 3)
  # u = 4   (0, 4)  (1, 4)  (2, 4)  (3, 4)
  
  SKHTtrainData0 <- createpoolugivent(data = trainData, id = "ID", 
                                      timefixedcov = c("X1","X2","X7"), 
                                      timevaryingcov = c("X3","X4","X5","X6","X8","X9","X10"), 
                                      y = "Event", period = "Time", hazard = "h",
                                      currentperiod = 0)
  
  hhat <- sapply(1:nperiod, function(ti){
    dat <- SKHTtrainData0[SKHTtrainData0$Time == ti, ]
    h_hati <- sum(dat$Event) / NROW(dat)
    h_hati
  })
  
  # The i-th cell: t = i - 1, u = i, i+1, ..., nperiod
  hbench <- lapply(1:nrow(utlist), function(kk){
    rep(hhat[utlist$u[kk]], ntest)
  })
  rm(SKHTtrainData0)
  rm(hhat)
  gc()
  
  print("Done -- 1. Current benchmark with no covariates")
  
  
  ################# == 2. Separate with Hellinger (t, u) == #################
  hsephlgr <- rep(list(0), nperiod * (nperiod + 1)/2)
  kk <- 1
  for (currentperiod in 0:(nperiod - 1)){
    for (predperiod in (currentperiod + 1):nperiod){
      # Create the training dataset
      septrainData <- createseparate(data = trainData, id = "ID", 
                                     timefixedcov = timefixedcov, 
                                     timevaryingcov = timevaryingcov, 
                                     y = "Event", period = "Time", hazard = "h",
                                     currentperiod = currentperiod,
                                     predperiod = predperiod)
      
      
      if (sum(septrainData$Event) == 0){
        hsephlgr[[kk]] <- rep(0, ntest)
      } else if (all(septrainData$Event == 1)) {
        hsephlgr[[kk]] <- rep(1, ntest)
      } else {
        septrainData$Event <- factor(septrainData$Event, levels = c("0", "1"))
        
        # Train the model with Hellinger splitting rule
        RFsephlgr <- ranger(formula = CFsep, data = septrainData, replace = FALSE, probability = TRUE,
                            splitrule = "hellinger", num.trees = ntree, num.threads = 1, oob.error = FALSE, verbose = FALSE)
        # Prediction
        hsephlgr[[kk]] <- predict(RFsephlgr, data = testDatatu[[kk]][, c(timefixedcov, timevaryingcov)], num.threads = 1, verbose = FALSE)$predictions[, 2]
      } 
      
      kk <- kk + 1
    }
  }
  rm(septrainData)
  rm(RFsephlgr)
  rm(kk)
  gc()
  print("Done -- 2. Separate with Hellinger (t, u)")
  
  
  ################# == 3. SKHT with Hellinger (t) == #################
  hSKHThlgr <- rep(list(0), nperiod * (nperiod + 1)/2)
  kk <- 1
  for (currentperiod in 0:(nperiod - 1)){
    # Generate training data for poolugivent
    SKHTtrainData <- createpoolugivent(data = trainData, id = "ID", 
                                       timefixedcov = c("X1","X2","X7"), 
                                       timevaryingcov = c("X3","X4","X5","X6","X8","X9","X10"), 
                                       y = "Event", period = "Time", hazard = "h",
                                       currentperiod = currentperiod)
    if (sum(SKHTtrainData$Event) == 0){
      for (predperiod in (currentperiod + 1):nperiod){
        hSKHThlgr[[kk]] <- rep(0, ntest)
        kk <- kk + 1
      }
    } else if (all(SKHTtrainData$Event == 1)){
      for (predperiod in (currentperiod + 1):nperiod){
        hSKHThlgr[[kk]] <- rep(1, ntest)
        kk <- kk + 1
      }
    } else {
      SKHTtrainData$Event <- factor(SKHTtrainData$Event, levels = c("0", "1"))
      # Train -- Hellinger splitting rule
      RFSKHThlgr <- ranger(formula = CFSKHT, data = SKHTtrainData, replace = FALSE, probability = TRUE, 
                           splitrule = "hellinger", num.trees = ntree, num.threads = 1, oob.error = FALSE, verbose = FALSE)
      
      for (predperiod in (currentperiod + 1):nperiod){
        # Prediction
        hSKHThlgr[[kk]] <- predict(RFSKHThlgr, data = testDatatu[[kk]][, c("Time", timefixedcov, timevaryingcov)], num.threads = 1, verbose = FALSE)$predictions[, 2]
        
        kk <- kk + 1
      }
    }
  }
  
  rm(SKHTtrainData)
  rm(RFSKHThlgr)
  rm(kk)
  gc()
  
  
  print("Done -- 3. SKHT with Hellinger (t)")
  
  
  ################# == 4.5.6. Superpp with Baseline Hellinger & Hellinger & DTPOmain == #################
  hsupp0hlgr <- rep(list(0), nperiod * (nperiod + 1)/2)
  hsupphlgr <- rep(list(0), nperiod * (nperiod + 1)/2)
  hsuppDTPO <- rep(list(0), nperiod * (nperiod + 1)/2)
  
  # Generate training data for Superpp 
  supptrainData <- createsuperpp(data = trainData, id = "ID", 
                                 timefixedcov = c("X1","X2","X7"), 
                                 timevaryingcov = c("X3","X4","X5","X6","X8","X9","X10"), 
                                 y = "Event", period = "Time", hazard = "h")
  supptrainData$Event <- factor(supptrainData$Event, levels = c("0", "1"))
  
  
  supp0trainData <- createpoolugivent(data = trainData, id = "ID", 
                                      timefixedcov = c("X1","X2","X7"), 
                                      timevaryingcov = c("X3","X4","X5","X6","X8","X9","X10"), 
                                      y = "Event", period = "Time", hazard = "h",
                                      currentperiod = 0)
  supp0trainData = subset(supp0trainData, select = -current)
  supp0trainData <- createsuperpp(data = supp0trainData, id = "ID", 
                                  timefixedcov = c("X1","X2","X7"), 
                                  timevaryingcov = c("X3","X4","X5","X6","X8","X9","X10"), 
                                  y = "Event", period = "Time", hazard = "h")
  supp0trainData$Event <- factor(supp0trainData$Event, levels = c("0", "1"))
  
  
  # Train the model
  RFsupp0hlgr <- ranger(formula = CFsupp, data = supp0trainData, replace = FALSE, probability = TRUE,
                        splitrule = "hellinger", num.trees = ntree, num.threads = 1, oob.error = FALSE, verbose = FALSE)
  RFsupphlgr <- ranger(formula = CFsupp, data = supptrainData, replace = FALSE, probability = TRUE,
                       splitrule = "hellinger", num.trees = ntree, num.threads = 1, oob.error = FALSE, verbose = FALSE)
  benchDTPO <- glm(formula = CFsupp, data = supptrainData, family = "binomial")
  
  kk <- 1
  for (currentperiod in 0:(nperiod - 1)){
    for (predperiod in (currentperiod + 1):nperiod){
      # Prediction
      hsupphlgr[[kk]] <- predict(RFsupphlgr, data = testDatatu[[kk]][, c("current", "Time", timefixedcov, timevaryingcov)], num.threads = 1, verbose = FALSE)$predictions[, 2]
      hsuppDTPO[[kk]] <- as.numeric(predict(benchDTPO, newdata = testDatatu[[kk]][, c("current", "Time", timefixedcov, timevaryingcov)], 
                                            type = "response"))
      hsupp0hlgr[[kk]] <- predict(RFsupp0hlgr, data = testDatatu[[kk]][, c("current", "Time", timefixedcov, timevaryingcov)], num.threads = 1, verbose = FALSE)$predictions[, 2]
      kk <- kk + 1
    }
  }
  
  rm(supptrainData)
  rm(RFsupp0hlgr)
  rm(RFsupphlgr)
  rm(benchDTPO)
  rm(kk)
  gc()
  
  print("Done -- 4.5.6. Superpp with Baseline Hellinger & Hellinger & DTPOmain")
  rm(trainData)
  rm(testDatatu)
  gc()
  
  ################# == Evaluation == #################
  acc <- sapply(1:nrow(utlist), function(kk){
    v <- rep(0, 6)
    v[1] <- evalfunc(h = h[[kk]], hhat = hbench[[kk]], method = meas)
    v[2] <- evalfunc(h = h[[kk]], hhat = hsephlgr[[kk]], method = meas)
    v[3] <- evalfunc(h = h[[kk]], hhat = hSKHThlgr[[kk]], method = meas)
    v[4] <- evalfunc(h = h[[kk]], hhat = hsupphlgr[[kk]], method = meas)
    v[5] <- evalfunc(h = h[[kk]], hhat = hsuppDTPO[[kk]], method = meas)
    v[6] <- evalfunc(h = h[[kk]], hhat = hsupp0hlgr[[kk]], method = meas)
    names(v) <- c("bench", "sephlgr", "skhthlgr", "supphlgr", "suppDTPO", "supp0hlgr")
    return(v)
  }) %>% 
    t() %>% 
    as_tibble() %>% 
    mutate(config_index = 1:nrow(utlist)) %>%
    left_join(utlist %>% rowid_to_column('config_index')) %>% 
    select(-"config_index")
  print("Done -- Evaluation")
  return(acc)
}