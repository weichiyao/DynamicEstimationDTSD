###### GOAL
# Create the data set for the pooling over u given t approach -- pooling all u values together for a given t

###### DETAILS
# for a given t, we use all the covariates information available up to time t to build a forest 
# to get estimations of the hazard function at u = t + 1, ..., T
# 1) Bou-Hamad et al. (2009) and Bou-Hamad et al. (2011a), 2) Schmid et al. (2016) 
# See also chapter 6 of Tutz and Schmid (2016)

###### INPUT

# data = data frame in person-period (i.e. augmented, one line per subject per period) format.
# timefixedcov = vector of names of the time fixed covariates.
# timevaryingcov = vector of names of the timevarying covariates.
# y = name of the target (the 0-1 event indicator).
# period = name of the variable representing the period.
# id = name of subject identifier.
# currentperiod = moment (period) at which we want to make the predictions. We can use the covariates information up to that period. Possible values are 0,1,...,maxperiod-1.
# lagtouse = vector indicating which lag to use for the time-varying covariates, same for each subject
#            (0 means only the current values, (0,1) means the current and lag 1 values etc.).
#            Default NULL indicating only the current values; if it is larger the number of available lags, 
#            then all the available lags are used (baseline to the current period).
# missinglag = What to do with lags that are before baseline. (1 = input the baseline value ; 0 = put a NA)
# hazard = name of the variable representing the true hazard (if available).

###### OUTPUT

# A data frame to be used for the pool over u given t approach


createpoolugivent <- function(data,
                              id,
                              timefixedcov = NULL,
                              timevaryingcov = NULL,
                              y,
                              period,
                              hazard = NULL,
                              currentperiod, 
                              baseline = FALSE,
                              lagtouse = 0,
                              missinglag = 1)
{
  maxperiod <- max(data[, period])
  # t = 0, ..., T-1, u = t+1, ..., T => u > t+1
  if ((currentperiod + 1) > maxperiod) stop(" currentperiod + 1 must be <= maximum observed period ")
  

  datfinal <- lapply(1:(maxperiod - currentperiod), function(i){
    datfinali <- createseparate(data = data, id = id, timefixedcov = timefixedcov, timevaryingcov = timevaryingcov, y = y, 
                                period = period, currentperiod = currentperiod, predperiod = currentperiod + i, 
                                lagtouse = lagtouse, missinglag = missinglag, hazard = hazard)
    datfinali
  })
  datfinal <- plyr::rbind.fill(datfinal) # library(plyr)

  # Sort by vector name id then period
  # datfinal[with(datfinal, order(id, "predict")),]
  # Does not work if id != "id"
  
  # Sort by column index [1] then [3]
  datfinal <- datfinal[order(datfinal[, id], datfinal[, period]), ]
  
  datfinal$I <- 1:nrow(datfinal)
  
  if (currentperiod == 0 & baseline) {
    finalI <- sapply(unique(datfinal[, id]), function(jj){
      max(datfinal[datfinal[, id] == jj, ]$I)
    })
    datfinal <- datfinal[finalI, ]
    datfinal$I <- 1:nrow(datfinal)
  }
  rownames(datfinal) <- NULL
  datfinal
}
