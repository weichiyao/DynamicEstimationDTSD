# 2011 paper:
# The sample has 1143 firms, 189 of them went bankrupt during the study period. 
# However, since 174 of the 189 bankruptcies occurred between years 3 and 8 after 
# the IPO, only this 6-year period is retained for the final analysis. The 15 
# remaining bankruptcies, which are scattered among the eight remaining years, 
# do not convey enough information to allow accurate estimations.

# Code to prepare the data for the analysis

# See the Excel file also.
precreate_realdata <- function(folder, dnum = FALSE){
  ########### Set your own path here!!!!!!!!!!!
  
  
  dat <- data.frame(read.table(paste(folder,"bank1.txt", sep = ""), header= TRUE, sep = ";"))
  # summary(dat)
  
  # datCIC <- data.frame(read.table(paste(pa,"cicdata.txt", sep = ""), header = TRUE, sep = " "))
  
  # Important. In this data set:
  # year = the year for the covariates values
  # period = the period for the covariates values
  # yearfory = the year for the response y
  # periodfory = the period for the response y
  
  # we don't need those variables
  dat$d1=NULL
  dat$d2=NULL
  dat$d3=NULL
  dat$d4=NULL
  dat$d5=NULL
  dat$d6=NULL
  
  # create the time independent covariate ipoyear
  
  dat$ipoyear=0
  
  for(i in 1:nrow(dat))
  {
    if(dat$period[i]==1){dat$ipoyear[i]=dat$year[i]}
    else{dat$ipoyear[i]=dat$ipoyear[i-1]}
  }
  
  # create a variable that says what is the maximum value of u for that company
  
  dat$maxu=0
  
  for(i in 1:nrow(dat))
  {
    if(dat$period[i]==1)
    {
      if(i!=1)
      {
        dat$maxu[ind]=dat$u[i-1]
      }
      ind=i
    }
    if(dat$period[i]>1)
    {
      ind=c(ind,i)
    }
    if(i==nrow(dat))
    {
      dat$maxu[ind]=dat$u[nrow(dat)]
    }
  }
  
  # the 1143 company names
  # uniqueid <- unique(dat$Company_Name)
  
  if (dnum) {
    datdnum <- data.frame(read.table(paste(pa,"dnumdata.txt", sep = ""), header= TRUE, sep = " "))
    
    uniqueid <- unique(dat$Company_Name)
    dat$x11raw <- 0
    dat$x11 <- 0
    for (id in uniqueid){
      dnum <- floor(datdnum[datdnum$Company_Name == id, ]$DNUM / 100)
      dat[dat$Company_Name == id, ]$x11raw <- dnum
      if (dnum <= 9) {
        dat[dat$Company_Name == id, ]$x11 <- "A"
      } else {
        if (dnum <= 14) {
          dat[dat$Company_Name == id, ]$x11 <- "B"
        } else {
          if (dnum <= 17) {
            dat[dat$Company_Name == id, ]$x11 <- "C"
          } else {
            if (dnum <= 39) {
              dat[dat$Company_Name == id, ]$x11 <- "D"
            } else {
              if (dnum <= 49) {
                dat[dat$Company_Name == id, ]$x11 <- "E"
              } else {
                if (dnum <= 51) {
                  dat[dat$Company_Name == id, ]$x11 <- "F"
                } else {
                  if (dnum <= 59) {
                    dat[dat$Company_Name == id, ]$x11 <- "G"
                  } else {
                    if (dnum <= 67) {
                      dat[dat$Company_Name == id, ]$x11 <- "H"
                    } else {
                      # if (dnum <= 89) {
                      #   dat[dat$Company_Name == id, ]$x11 <- "I"
                      # } else {
                        dat[dat$Company_Name == id, ]$x11 <- "I"
                        
                      #}
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return(dat)
}

# This is the function to use to create the training and test data sets with the separate method.
# The SKHT and superpp taining data sets are obtained by stacking them ion the right way as explained below.
# This is the function to use to create the training and test data sets with the separate method.
# The SKHT and superpp taining data sets are obtained by stacking them ion the right way as explained below.

## Denis's original version
createseparatebank <- function(dat,uniqueid,year,t,u){
  # dat = data set created above
  # uniqueid = vector created above
  # year = must be 1997, 1998 or 1999 (see the setup in the paper)
  # t = must be 0, 1, 2 or 3 (see table 7 in the paper)
  # u = must be 3, 4, 5, 6 (see table 7 in the paper)
  
  # OUTPUT: the train data set and the test data set for the separate method for these values of t and u and for that year
  #  only the current values of the covariates are present (no lags)
  
  out=dat[1,]
  
  # create the data for all the (t,u) pairs
  for(i in 1:length(uniqueid))
  {
    dati=dat[dat$Company_Name==uniqueid[i],]
    if(dati$maxu[1]<u){next}
    else
    {
      out=rbind(out,dati[dati$t==t,])
      out[nrow(out),c("y","yearfory","periodfory","u")]=dati[dati$u==u,c("y","yearfory","periodfory","u")]
    }
  }
  out[-1,]
  
  #print(summary(out))
  #print(dim(out))
  
  # create the train and test data sets
  
  print(year)
  datrain=out[out$yearfory<=year,]
  datest=out[out$year==year,]
  
  out1=vector("list",2)
  out1[[1]]=datrain
  out1[[2]]=datest
  names(out1)=c("datrain","datest")
  out1
  
}
# ## Weichi's cleaner version
# createseparatebank <- function(dat,uniqueid,year,t,u){
#   # dat = data set created above
#   # uniqueid = vector created above
#   # year = must be 1997, 1998 or 1999 (see the setup in the paper)
#   # t = must be 0, 1, 2 or 3 (see table 7 in the paper)
#   # u = must be 3, 4, 5, 6 (see table 7 in the paper)
#   
#   # OUTPUT: the train data set and the test data set for the separate method for these values of t and u and for that year
#   #  only the current values of the covariates are present (no lags)
#   
#   
#   coltake <- c("Company_Name", "y" , 
#                "x1", "x2", "x3", "x4", "x5", "x6", 
#                "x7", "x8", "x9", "x10",
#                "year", "yearfory", "t", "u", "maxu")
#   dat <- dat[, coltake]
#   out <- dat[1, ]
#   
#   # create the data for all the (t,u) pairs
#   for(i in 1:length(uniqueid))
#   {
#     dati <- dat[dat$Company_Name == uniqueid[i], ]
#     if(dati$maxu[1] < u){next}
#     else
#     {
#       out <- rbind(out, dati[dati$t == t, ])
#       out[nrow(out), c("y", "yearfory", "u")] <- dati[dati$u == u, c("y", "yearfory", "u")]
#     }
#   }
#   
#   out[-1, ]
#   
#   #print(summary(out))
#   #print(dim(out))
#   coltake <- c("Company_Name", "y" , 
#                "x1", "x2", "x3", "x4", "x5", 
#                "x6", "x7", "x8", "x9", "x10",
#                "t", "u")
#   # create the train and test data sets
#   datrain <- out[out$yearfory<=year, coltake]
#   datest <- out[out$year==year, coltake] 
#   
#   out1 <- vector("list",2)
#   out1[[1]] <-datrain
#   out1[[2]] <- datest
#   names(out1) <- c("datrain", "datest")
#   out1
#   
# }