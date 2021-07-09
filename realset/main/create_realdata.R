########### Set your own path here
pa <- "./codes/realset/"

dat <- data.frame(read.table(paste(pa,"bank1.txt", sep = ""), header= TRUE, sep = ";"))
# summary(dat)

datCIC <- data.frame(read.table(paste(pa,"cicdata.txt", sep = ""), header = TRUE, sep = " "))

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
uniqueid <- unique(dat$Company_Name)

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

## Weichi's cleaner version
createseparatebank <- function(dat,uniqueid,year,t,u){
  # dat = data set created above
  # uniqueid = vector created above
  # year = must be 1997, 1998 or 1999 (see the setup in the paper)
  # t = must be 0, 1, 2 or 3 (see table 7 in the paper)
  # u = must be 3, 4, 5, 6 (see table 7 in the paper)
  
  # OUTPUT: the train data set and the test data set for the separate method for these values of t and u and for that year
  #  only the current values of the covariates are present (no lags)
  
  
  coltake <- c("Company_Name", "y" , 
               "x1", "x2", "x3", "x4", "x5", "x6", 
               "x7", "x8", "x9", "x10",
               "year", "yearfory", "t", "u", "maxu")
  dat <- dat[, coltake]
  out <- dat[1, ]
  
  # create the data for all the (t,u) pairs
  for(i in 1:length(uniqueid))
  {
    dati <- dat[dat$Company_Name == uniqueid[i], ]
    if(dati$maxu[1] < u){next}
    else
    {
      out <- rbind(out, dati[dati$t == t, ])
      out[nrow(out), c("y", "yearfory", "u")] <- dati[dati$u == u, c("y", "yearfory", "u")]
    }
  }

  out[-1, ]
  
  #print(summary(out))
  #print(dim(out))
  coltake <- c("Company_Name", "y" , 
               "x1", "x2", "x3", "x4", "x5", 
               "x6", "x7", "x8", "x9", "x10",
               "t", "u")
  # create the train and test data sets
  datrain <- out[out$yearfory<=year, coltake]
  datest <- out[out$year==year, coltake] 
  
  out1 <- vector("list",2)
  out1[[1]] <-datrain
  out1[[2]] <- datest
  names(out1) <- c("datrain", "datest")
  out1
  
}
# try the function
t=1
u=5
year=1997
out=createseparatebank(dat,uniqueid,year,t,u)
summary(out$datrain)
summary(out$datest)


# Check if the data are correct by reproducing the results in table 8 of the paper

checkdata=function(dat,uniqueid)
{
  out=data.frame(matrix(1,ncol=7,nrow=1))
  names(out)=c("year","t","u","natrisktrain","nbankrupttrain","natrisktest","nbankrupttest")
  
  for(year in 1997:1999)
  {
    for(t in 0:3)
    {
      for(u in (t+3):6)
      {
        dat1=createseparatebank(dat,uniqueid,year,t,u)
        datrain=dat1$datrain
        datest=dat1$datest
        out=rbind(out,c(year,t,u,nrow(datrain),sum(datrain$y),nrow(datest),sum(datest$y)))
      }
    }
  }
  out[-1,]
}


out=checkdata(dat,uniqueid)
summary(out)




createskhtbank <- function(sepdat, uniqueid, year, t){
  # dat = data set created for separate
  # uniqueid = vector created above
  # year = must be 1997, 1998 or 1999 (see the setup in the paper)
  # t = must be 0, 1, 2 or 3 (see table 7 in the paper)
  
  # OUTPUT: the train data set for the SKHT method for these t values and for that year
  #  only the current values of the covariates are present (no lags)
  
  datfinal <- lapply((t + 3):6, function(i){
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
  
  # Sort by vector name id then period
  # datfinal[with(datfinal, order(id, "predict")),]
  # Does not work if id != "id"
  
  # Sort by column index [1] then [3]
  datfinal <- out[order(out[, uniqueid], out[, period]), ]
  
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
}

