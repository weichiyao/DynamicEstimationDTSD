boxplotfunc <- function(scenario = c("0TI2TV", "0TI4TV", "1TI4TV", "2TI1TV", "2TI4TV"),
                        SNR = c("low", "high"), 
                        num.train = c(200, 1000, 5000),
                        num.test = 1000,
                        censor.rate = c("10%", "50%"),
                        autocorrtype = c("weak", "strong"),
                        distribution = c("Exp", "WI", "Gtz"),
                        relationship = c("linear", "nonlinear", "interaction"),
                        measurement = c("adist", "alor", "cindex"),
                        num.period = c(4, 8), 
                        nsim = 100){
  ################################################################################
  ## == GOAL:
  ## This function is create the boxplot after getting measurment values from
  ## predfunc.R
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
  ## Boxplots 
  ################################################################################
  df_gp <- lapply(1:nsim, function(ll){
    acc <- predfunc(scenario = scenario,
                    SNR = SNR, 
                    num.train = num.train,
                    num.test = num.test,
                    censor.rate = censor.rate,
                    autocorrtype = autocorrtype,
                    distribution = distribution,
                    relationship = relationship,
                    measurement = measurement,
                    num.period = num.period)
  }) %>% 
    bind_rows(.id = "simulation") %>% 
    pivot_longer(-c(t, u, simulation), names_to = "method", values_to = "value")
    
    
  if (measurement == "cindex") {
    df_gp <- df_gp %>% filter(method != "bench")
    nmethod <- 5
  } else {
    nmethod <- 6
  }
  
  ## Given n, weak v.s. strong: Organize the plot by t: In each plot, each block shows a specific u orded by method.
  p <- lapply(0:(num.period - 1), function(tt){
    df_gp %>% 
      dplyr::filter(t == tt) %>%
      mutate(method = factor(method)) %>%
      mutate(newgroup = interaction(u, method)) %>%
      mutate(newgroup = fct_reorder(newgroup, u)) %>%
      ggplot() + 
      geom_boxplot(aes(x = newgroup, 
                       y = value, fill = method)) + 
      labs(x = "method", y = "measurement") +
      geom_vline(xintercept = nmethod + 0.5 + nmethod * c(0 : (num.period - tt - 1)),
                 linetype = 4, colour = "black") +
      theme_classic() + 
      ggtitle(sprintf("t = %1.0f", tt)) +
      theme(plot.title = element_text(size = 16, face = "bold"),
            axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0),
            legend.position = "top")
  })
  
  p1 <- ggarrange(p[[1]], p[[2]], p[[3]], p[[4]], nrow = 2, ncol = 2)
  print(p1)
  if (num.period == 8){
    p2 <- ggarrange(p[[5]], p[[6]], p[[7]], p[[8]], nrow = 2, ncol =2)
    print(p2)
  } 
  
}

