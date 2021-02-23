create_tblplot <- function(dblfinal, horizon, measurement){
  if (horizon == "=1") {
    utdiff <- 1
    tblplot <- lapply(measurement, function(meas){
      meascolname <- paste(c(meas, utdiff), collapse = "_")
      Formula <- adist_1 ~ scenario + model + distribution + signal + autocorrtype + crate + ntrain
      Formula <- as.formula(paste(c(meascolname, "~", Formula[[3]]), collapse = " "))
      mod  <- lm(Formula, data = dblfinal)
      
      # mean of all
      MEall <- mean(dblfinal[[meascolname]])
      
      MEdf <- as.data.frame(effects::allEffects(mod))
      # "scenario"     "model"        "distribution" "signal"       "autocorrtype" "crate"        "ntrain"
      namMEdf <- names(MEdf)
      namMEdfnew <- c("Scenario", "Relationship", "Distribution", "SNR", "Autocorrelation", "Censor rate", "Sample size")
      
      lapply(1:length(MEdf), function(mm){
        MEdf[[mm]][, 1:2] %>%
          mutate(cat = namMEdfnew[mm]) %>%
          rename("subcat" = namMEdf[mm],
                 "Mean" = fit)
      }) %>%
        bind_rows() %>%
        select(cat, subcat, Mean) %>%
        as_tibble() %>%
        mutate(measurement = meas) %>%
        mutate(meanall = MEall)
    }) %>%
      bind_rows()  %>%
      mutate(measurement = recode(measurement,
                                  "adist" = "ADIST",
                                  "alor" = "ALOR",
                                  "cindex" = "C-index"))
  } else if (horizon == ">1") {
    tblplot <- lapply(measurement, function(meas){
      meascolname <- sprintf("%s_%1.0f", meas, 2:4)
      measvalue <- rowMeans(as.data.frame(dblfinal)[meascolname])
      dblres <- dblfinal %>%
        mutate(measvalue = measvalue)
      
      Formula <- measvalue ~ scenario + model + distribution + signal + autocorrtype + crate + ntrain
      mod  <- lm(Formula, data = dblres)
      
      # mean of all
      MEall <- mean(dblres$measvalue)
      
      MEdf <- as.data.frame(effects::allEffects(mod))
      # "scenario"     "model"        "distribution" "signal"       "autocorrtype" "crate"        "ntrain"
      namMEdf <- names(MEdf)
      namMEdfnew <- c("Scenario", "Relationship", "Distribution", "SNR", "Autocorrelation", "Censor rate", "Sample size")
      
      lapply(1:length(MEdf), function(mm){
        MEdf[[mm]][, 1:2] %>%
          mutate(cat = namMEdfnew[mm]) %>%
          rename("subcat" = namMEdf[mm],
                 "Mean" = fit)
      }) %>%
        bind_rows() %>%
        select(cat, subcat, Mean) %>%
        as_tibble() %>%
        mutate(measurement = meas) %>%
        mutate(meanall = MEall)
    }) %>%
      bind_rows()  %>%
      mutate(measurement = recode(measurement,
                                  "adist" = "ADIST",
                                  "alor" = "ALOR",
                                  "cindex" = "C-index"))
  } else{
    stop("Estimation horizon can only be '=1' or '>1'!")
  }
  
  p1 <- tblplot %>%
    mutate(subcat = factor(subcat, levels = unique(subcat))) %>% # this is important to made ntrain = 200, 1000, 5000 in order
    # arrange(match(subcat, c("200", "1000", "5000"))) %>%
    ggplot(aes(subcat, Mean)) +
    geom_point(group = 1, color = "steelblue") +
    geom_line(group = 1, color = "steelblue") +
    facet_grid(cols = vars(cat), rows = vars(measurement), scales = "free") +
    geom_hline(aes(yintercept = meanall), linetype = "dashed", color = "#999999") + 
    geom_hline(aes(yintercept = 0), linetype = "solid", color = "gray40") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank()) 
  
  print(p1) # has to print, or do nothing 
}