###### GOAL
# Create the data set for the superpp approach -- pooling all the information for all combinations of (t, u) at once

###### DETAILS
# stacking the data sets from all values of t in the Bou-Hamad et al. (2009) or the Schmid et al. (2016) method

###### INPUT:

# data = data frame in person-period (i.e. augmented, one line per subject per period) format.
# timefixedcov = vector of names of the time fixed covariates.
# timevaryingcov = vector of names of the timevarying covariates.
# y = name of the target (the 0-1 event indicator).
# period = name of the variable representing the period.
# id = name of subject identifier.
# currentperiod = moment (period) at which we want to make the predictions. We can use the covariates information up to that period. Possible values are 0,1,...,maxperiod-1.
# lagtouse = vector indicating which lag to use for the time-varying covariates. (0 means only the current values, (0,1) means the current and lag 1 values etc.).
# missinglag = What to do with lags that are before baseline. (1 = input the baseline value ; 0 = put a NA)
# hazard = name of the variable representing the true hazard (if available).


## OUTPUT

# A data frame to be used for the superpp method


createsuperpp <- function(data, 
                          id,
                          timefixedcov = NULL,
                          timevaryingcov = NULL,
                          y,
                          period,
                          hazard = NULL,
                          lagtouse = 0,
                          missinglag = 1)
{
  
  maxperiod <- max(data[, period])
 
  datfinal <- lapply(0:(maxperiod - 1), function(j){
    datfinalj <- createpoolugivent(data = data, id = id, 
                                   timefixedcov = timefixedcov, timevaryingcov = timevaryingcov, 
                                   y = y,
                                   period = period, currentperiod = j, 
                                   lagtouse = lagtouse, missinglag = missinglag, hazard = hazard)
    datfinalj
  })
  
  datfinal <- plyr::rbind.fill(datfinal) # library(plyr)
  
  # Sort by vector name id then period
  # datfinal[with(datfinal, order(id, period)),]
  # Does not work if id != "id"
  
  # Sort by column index [1] then [3]
  datfinal[order(datfinal[, id], datfinal[, "current"], datfinal[, period]), ]
  
  datfinal
  
}
