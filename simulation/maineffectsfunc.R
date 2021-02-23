maineffectsfunc <- function(nsim = 100, num.test = 1000, horizon = c("=1", ">1")){
  ################################################################################
  ## == GOAL:
  ## This function is create the maineffect plots when number of period equals to 4
  ## from predfunc.R by comparing across all the following factors:
  ##    scenario     -- different TI/TV combinations: c("0TI2TV", "0TI4TV", "1TI4TV", "2TI1TV", "2TI4TV")
  ##    SNR          -- signal to noise ratio: c("low", "high")
  ##    censor.rate  -- censoring rate: c("10%", "50%"),
  ##    autocorrtype -- autocorrelation type: c("weak", "strong"),
  ##    distribution -- survival distribution: c("Exp", "WI", "Gtz"),
  ##    relationship -- survival relationship: c("linear", "nonlinear", "interaction"),
  ##    num.period   -- number of periods: c(4, 8) 
  ##    num.train    -- number of training samples
  ## one for each measurement ("adist", "alor", "cindex")
  ##
  ## == INPUT:
  ## nsim         -- number of simulations
  ## horizon      -- estimation horizon: "=1", ">1"
  ## num.test     -- number of test samples
  ##
  ## == OUTPUT:
  ## main effects plots
  ################################################################################
  nntrain <- c(200, 1000, 5000)
  ddist <- c("Exp", "WI", "Gtz")
  ccrate <- c("10%", "50%")
  mmodel <- c("linear", "nonlinear", "interaction")
  aautocorr <- c("strong", "weak")
  mmeasr <- c("alor", "adist", "cindex")
  ssnr <- c("high", "low")
  sscen <- c("0TI4TV", "0TI2TV", "1TI4TV", "2TI1TV", "2TI4TV")
  # To compare all scenarios, T = 4
  dblfinal <- lapply(1:5, function(ss){ # different scenario
    lapply(1:2, function(cc){ # different censoring rate
      lapply(1:2, function(bb){ # boosted = TRUE vs. boosted = FALSE
        lapply(1:3, function(dd){ # different distributions
          lapply(1:3, function(nn){ # different n
            lapply(1:3, function(mm){ # different survival relationships
              lapply(1:2, function(aa){ # different autocorrelations
                lapply(1:3, function(ee){ # different measurements: "alor", "adist", "cindex"
                  lapply(1:nsim, function(ll){
                    predfunc(scenario = sscen[ss],
                             SNR = ssnr[bb], 
                             num.train = nntrain[nn],
                             num.test = num.test,
                             censor.rate = ccrate[cc],
                             autocorrtype = aautocorr[aa],
                             distribution = ddist[dd],
                             relationship = mmodel[mm],
                             measurement = mmeasr[ee],
                             num.period = 4) %>% 
                      select(sephlgr, supphlgr, t, u) %>%
                      mutate(utdiff = u - t,
                             measurement = mmeasr[ee],
                             distribution = ddist[dd],
                             model = mmodel[mm],
                             autocorrtype = aautocorr[aa],
                             crate = ccrate[cc],
                             ntrain = nntrain[nn],
                             scenario = sscen[ss],
                             signal = ssnr[bb]) %>%
                      select(sephlgr, supphlgr, utdiff, measurement, distribution, model, autocorrtype, crate, ntrain, scenario, signal)
                  }) %>% bind_rows(.id = "simulation")
                }) %>% bind_rows()
              }) %>% bind_rows()
            }) %>% bind_rows()
          }) %>% bind_rows()
        }) %>% bind_rows()
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows() 
  
  dblfinal <- dblfinal %>% 
    pivot_longer(-c(measurement, scenario, model, distribution, signal, autocorrtype, crate, ntrain, utdiff, simulation), 
                 names_to = "method", values_to = "value") %>%
    select(c(method, measurement, scenario, model, distribution, signal, autocorrtype, crate, ntrain, utdiff, simulation, value))%>%
    group_by(method, measurement, scenario, model, distribution, signal, autocorrtype, crate, ntrain, utdiff) %>%
    summarise(medoversim = median(value)) %>%
    ungroup() %>%
    group_by(measurement, scenario, model, distribution, signal, autocorrtype, crate, ntrain, utdiff) %>%
    mutate(sepvssupp = c(0, -diff(medoversim))) %>%
    filter(method == "supphlgr") %>%
    select(-c(method, medoversim)) %>%
    mutate(sepvssupp = round(sepvssupp, 4)) %>%
    pivot_wider(names_from = c(measurement, utdiff), values_from = sepvssupp) %>% 
    ungroup()
  
  dblfinal$ntrain = factor(dblfinal$ntrain, levels = c("200", "1000", "5000"))
  dblfinal$autocorrtype = factor(dblfinal$autocorrtype, levels = c("weak", "strong"))
  dblfinal$signal = factor(dblfinal$signal, levels = c("low", "high"))
  

  
  create_tblplot(dblfinal = dblfinal, 
                 horizon = horizon, 
                 measurement = c("adist", "alor", "cindex"))
  
}

