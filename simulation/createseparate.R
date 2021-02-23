###### GOAL
# Create the data set for the seperate approach -- to only use local information at each unique combination of (t, u)

###### DETAILS:
# For a given pair (t, u), we use only the subjects that are still alive and not censored at time u − 1 
# to build the forest to estimate the hazard at time u;
# we only use all the covariates information available up to time t to build the forests

# maxperiod = T
# currentperiod = t: the subject is currently alive at time t, t = {0, 1, ..., T - 1}
# predperiod = u: we are interested in estimating the hazard function at time u for the subject 
#                 who's currently alive at time t, u = {t + 1, ..., T}
# t < u!

###### INPUT
# data = data frame in person-period (i.e. augmented, one line per subject per period) format.
# timefixedcov = vector of names of the time fixed covariates.
# timevaryingcov = vector of names of the timevarying covariates.
# y = name of the target (the 0-1 event indicator).
# period = name of the variable representing the period.
# id = name of the subject identifier.
# lagtouse = vector indicating which lag to use for the time-varying covariates, same for each subject
#            (0 means only the current values, (0,1) means the current and lag 1 values etc.).
#            Default NULL indicating only the current values; if it is larger the number of available lags, 
#            then all the available lags are used (baseline to the current period).
# missinglag = What to do with lags that are before baseline. (1 = input the baseline value ; 0 = put a NA)
# hazard = name of the variable representing the true hazard (if available).

###### OUTPUT

# A data frame to be used for the separate approach

createseparate  <- function(data, 
                            id,
                            timefixedcov = NULL,
                            timevaryingcov = NULL,
                            y,
                            period,
                            hazard = NULL,
                            currentperiod,
                            predperiod,
                            lagtouse = 0,
                            missinglag = 1)
{
  
  maxperiod <- max(data[, period])
  
  if (predperiod > maxperiod) stop(" predperiod must be <= maximum observed period ")
  if (predperiod <= currentperiod) stop(" predperiod must be > currentperiod ")
  
  # use only the subjects that 
  # 1) are currently alive at time t -- automatically satisfied
  #    This is because in the person-period dataset, we only gives subject observations at time if it is still alive
  # 2) are still alive and not censored at time u − 1 -- we need to substract this part
  
  # number of subjects in the output dataset
  finalid <- unique(data[data[, period] == predperiod, id])
  nfinal <- length(finalid)
  
  if (is.null(timevaryingcov) | lagtouse == 0) {
    lagtouse <- 0
    ctvnam <- NULL 
  } 
  
  if (!is.null(timevaryingcov) & lagtouse > 0) {
    ctvnam <- as.vector(sapply(1:length(timevaryingcov), 
                               function(j) c(timevaryingcov[j], paste(timevaryingcov[j], "LAG", lagtouse[-1], sep = ""))))
  }
  
  if (!is.null(timevaryingcov) & lagtouse == 0){
    ctvnam <- as.vector(sapply(1:length(timevaryingcov), 
                               function(j) timevaryingcov[j]))
    
  }
  
  datafinal <- data.frame(matrix(0, nrow = nfinal, 
                                 ncol = length(timevaryingcov) * (max(lagtouse) + 1) + length(timefixedcov) + 5))
  names(datafinal) <- c("I", id, "current", period, y, ctvnam, timefixedcov)
  
  datafinal[, "I"] <- 1:nfinal
  datafinal[, id] <- finalid
  datafinal[, "current"] <- currentperiod
  datafinal[, period] <- predperiod
  
  datafinal[, c(timefixedcov, y)] <- data[data[, period] == predperiod, c(timefixedcov, y)]  # This gives subjects that are still alive and not censored at (u - 1)v
  # This is sufficient since if it is censored or dead at u - 1, then we will not have the new observations X(u - 1)
  
  if (!is.null(timevaryingcov)) {
    lagfinal <- lapply(1:nfinal, function(i) {
      lagfinali <- currentperiod + 1 - lagtouse
      lagfinali[lagfinali <= 0] <- ifelse(missinglag, 1, NA) # if missinglag = 1, then input the baseline value; if = 0, then NA
      lagfinali
    })
    for (i in 1:nfinal){
      # Lags to be used for each subject
      datafinal[i, ctvnam] <- unlist(data[data[, id] == finalid[i], timevaryingcov, drop = FALSE][lagfinal[[i]], ])
    }
  }
  
  if(!is.null(hazard)) {
    datafinal[, hazard] <- data[data[, period] == predperiod, hazard]
  }
  
  datafinal
  
}
