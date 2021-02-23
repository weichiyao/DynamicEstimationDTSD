library(truncdist)
library(survival)
library(expm)
library(ggplot2)
library(tidyverse)
library(ggpubr) # for ggarrange
library(ranger)


source("./utils.R")
source("./testdtv_gnrt.R")
source("./genvar.R")
source("./Timevarying_gnrt.R")
source("./traindtv_autocorr_gnrt.R")
source("./createsuperpp.R")
source("./createseparate.R")
source("./createpoolugivent.R")
source("./create_tblplot.R")
source("./predfunc.R")
source("./evalfunc.R")
source("./boxplotfunc.R")



## == INPUT (see boxplotfunc.R):
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
## Boxplots 
boxplotfunc(scenario = "2TI1TV",
            SNR = "high", 
            num.train = 50,
            num.test = 10,
            censor.rate = "50%",
            autocorrtype = "strong",
            distribution = "WI",
            relationship = "nonlinear",
            measurement = "alor",
            num.period = 4, 
            nsim = 3)

## == INPUT (see maineffectsfunc.R):
## nsim         -- number of simulations
## horizon      -- estimation horizon: "=1", ">1"
## num.test     -- number of test samples
##
## == OUTPUT:
## main effects plots
maineffectsfunc(nsim = 3, 
                num.test = 10,
                horizon = "=1")
