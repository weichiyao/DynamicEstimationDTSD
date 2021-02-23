##########################################################################################################
testdtv_gnrt = function(data, ntest, id = NULL, period = NULL, y = NULL){
  ############################################################
  # GOAL #
  # generates the test data sets, one for each t (given T > t)
  
  # INPUT #
  # data = person-period data (pseudo-subject observations)
  # ntest = size of the test data set (for each t = 0, ..., maxt - 1), a sample of size ntest is produced.
  # period = name of the variable representing the period/(discrete) time.
  # id = name of subject identifier.
  # y = name of the target (the 0-1 event indicator).
  
  # OUTPUT #
  # a list of maxt elements
  # element j (j = 1, ..., maxt) is the data frame of the test points when t = (j - 1)
  ############################################################
  
  maxt <- max(data[, period])    # number of periods
  allid <- unique(data[, id])    # unique id from the data
  ind <- sample(allid, ntest)    # sample of subject (their id)
  
  inddata <- which(data[, id] %in% ind)
  
  out <- list(data[inddata, ])
  
  for (j in 2:maxt) {
    tempdata <- data[data[, period] == j, ]
    indc <- tempdata[, id]
    # print(length(unique(indc)))
    
    # number of subjects in the sample provided (must be large enough to make sure that we will
    # have at least ntest subjects for each j=1, ..., maxt
    if (length(unique(indc)) < ntest) stop(sprintf("Number of subjects in %1.0f-th set from the provided sample: %1.0f < ntest!", 
                                                   j, length(unique(indc))))
    
    ind <- sample(indc, ntest)          # sample of subject (their id)    
    inddata <- which(data[, id] %in% ind)
    out[[j]] <- data[inddata, ]
  }
  out
}
