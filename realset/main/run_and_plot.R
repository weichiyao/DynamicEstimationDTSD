## This file gives slightly different pred_realdata function compared to the one
## given in temp.R  -- also compute vimp for superpp w/ RF hellinger method

library(tidyverse)
library(ranger)
setwd("/Users/wyao/Dropbox/RESEARCH/discreteRSF/codes/realset/source")
source("computeAUC_new.R")
source("plotgainchart_new.R")
source("plotroccurve_new.R")
source("computeroc.R")
source("computecumgain.R")
source("precreate_realdata.R")
create_realdata <- function(dat, savefilepath) {
  uniqueid <- unique(dat$Company_Name)
  out <- lapply(1997:1999, function(yyear){
    lapply(0:3, function(tt){
      lapply((tt + 3):6, function(uu){
        createseparatebank(dat, uniqueid, year = yyear, t = tt, u = uu) %>%
          bind_rows(.id = "type") 
      }) %>%
        bind_rows() 
    }) %>%
      bind_rows() 
  }) %>% 
    bind_rows(.id = "year")
  saveRDS(out, savefilepath)
  return(out)
}
pred_realdata <- function(dat, cov, savefilepath){
  ## Set the formula
  CFsep <- as.formula(paste("y~", paste(sprintf("x%1.0f", cov), collapse = "+")))
  CFSKHT <- as.formula(paste("y~", paste(sprintf("x%1.0f", cov), collapse = "+"), "+t"))
  CFsupp <- as.formula(paste("y~", paste(sprintf("x%1.0f", cov), collapse = "+"), "+t+u"))
  ntree <- 500L
  
  vimp <- vector(mode = "list", length = length(unique(dat$year)))
  saveRDS(vimp, savefilepath)
  allRES <- lapply(unique(dat$year), function(yyear){
    outyear <- dat %>%
      dplyr::filter(year == yyear)
    
    
    # Create the training set for superpp
    supptrain <- outyear %>%
      dplyr::filter(type == "datrain")
    
    # Train -- bench
    benchtrain <- supptrain %>%
      filter(t == 0)
    
    hhatbench <- sapply(unique(benchtrain$u), function(uu){
      dat <- benchtrain[benchtrain$u == uu, ]
      h_hati <- sum(dat$y) / NROW(dat)
      h_hati
    })
    
    # Train -- superpp with Hellinger splitting rule
    RFsupphlgr <- ranger(formula = CFsupp, 
                         data = supptrain, 
                         replace = FALSE, 
                         probability = TRUE,
                         splitrule = "hellinger", 
                         num.trees = ntree, 
                         num.threads = 1, 
                         # oob.error = FALSE, 
                         verbose = FALSE,
                         importance = "permutation")
    vimptemp <- importance(RFsupphlgr)
    vimptemp <- c(yyear, vimptemp)
    names(vimptemp)[1] <- "year"
    vimp <- readRDS(savefilepath)
    vimp[[as.numeric(yyear)]] <- vimptemp
    saveRDS(vimp, savefilepath)
    # Train -- superpp with DTPO
    benchDTPO <- glm(formula = CFsep, data = supptrain, family = "binomial",
                     na.action = "na.omit")
    
    # Train -- superpp with baseline
    supp0train <- supptrain
    uniqueid <- unique(supptrain$Company_Name)
    for (id in uniqueid){
      tm <- supptrain[supptrain$Company_Name == id, ]
      tm <- tm[tm$t == 0, ]
      # alltc <- which(supptrain[supptrain$Company_Name == id, ]$t > 0)
      supp0train[supp0train$Company_Name == id, sprintf("x%1.0f", cov)] <- 
        tm[1, sprintf("x%1.0f", cov)]
    }
    # Train -- superpp with baseline 
    RFsupp0hlgr <- ranger(formula = CFsupp, 
                          data = supp0train, 
                          replace = FALSE, 
                          probability = TRUE,
                          splitrule = "hellinger", 
                          num.trees = ntree, 
                          num.threads = 1, 
                          oob.error = FALSE, 
                          verbose = FALSE)
    
    ## for each t
    lapply(unique(outyear$t), function(tt){
      outyeart <- outyear %>% 
        dplyr::filter(t == tt)
      
      # Create the training set for skhtrain
      skhtrain <- outyeart %>%
        dplyr::filter(type == "datrain")
      skhtrain$y <- factor(skhtrain$y, levels = c("0", "1"))
      # Train -- SKHT with Hellinger splitting rule
      RFSKHThlgr <- ranger(formula = CFSKHT, 
                           data = skhtrain, 
                           replace = FALSE, 
                           probability = TRUE, 
                           splitrule = "hellinger", 
                           num.trees = ntree, 
                           num.threads = 1, 
                           oob.error = FALSE, 
                           verbose = FALSE)
      
      ## for each u
      lapply(unique(outyeart$u), function(uu){
        outyeartu <- outyeart %>%
          dplyr::filter(u == uu)
        
        # Create the training set for separate
        septrain <- outyeartu %>%
          dplyr::filter(type == "datrain")
        
        # Create the test set
        testData <- outyeartu %>%
          dplyr::filter(type == "datest")
        
        # Prediction for bench
        hbench <- rep(hhatbench[uu - 3 + 1], nrow(testData))
        
        ntest <- nrow(testData)
        # Prediction for separate
        if (sum(septrain$y) == 0){
          hsephlgr <- rep(0, ntest)
        } else if (all(septrain$y == 1)) {
          hsephlgr <- rep(1, ntest)
        } else {
          septrain$y <- factor(septrain$y, levels = c("0", "1"))
          
          # Train the model with Hellinger splitting rule
          RFsephlgr <- ranger(formula = CFsep, data = septrain, replace = FALSE, 
                              probability = TRUE,
                              splitrule = "hellinger", 
                              num.trees = ntree, 
                              num.threads = 1, 
                              oob.error = FALSE, 
                              verbose = FALSE)
          # Prediction
          hsephlgr <- predict(RFsephlgr, 
                              data = testData, 
                              num.threads = 1, 
                              verbose = FALSE)$predictions[, 2]
        } 
        
        # Prediction for SKHT
        hSKHThlgr <- predict(RFSKHThlgr, 
                             data = testData, 
                             num.threads = 1, 
                             verbose = FALSE)$predictions[, 2]
        
        # Prediction for superpp
        hsupphlgr <- predict(RFsupphlgr, 
                             data = testData,
                             num.threads = 1, 
                             verbose = FALSE)$predictions[, 2]
        
        # Prediction for superpp with baseline
        hsupp0hlgr <- predict(RFsupp0hlgr, 
                              data = testData,
                              num.threads = 1, 
                              verbose = FALSE)$predictions[, 2]
        
        # Prediction for superppDTPO
        hsuppDTPO <- as.numeric(predict(benchDTPO, 
                                        newdata = testData, 
                                        type = "response"))
        
        list(hsuppDTPO = hsuppDTPO, 
             hsupphlgr = hsupphlgr, 
             hsupp0hlgr = hsupp0hlgr,
             hSKHThlgr = hSKHThlgr,
             hsephlgr = hsephlgr,
             hbench = hbench,
             t = testData$t,
             u = testData$u,
             y = testData$y,
             year = testData$year) %>%
          as_tibble() 
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()
  
  vimp <- readRDS(savefilepath) %>% 
    bind_rows() #%>%
    #as.data.frame(as.numeric)
  RES <- list(vimp_supphlgr = vimp, pred = allRES)
  # write.csv(allRES, "./allRES.csv")
  saveRDS(RES, savefilepath)
  return(RES)
}
#########################################################################################################
pa <- "/Users/wyao/Dropbox/RESEARCH/discreteRSF/codes/realset/main/"
dat_pre_add_DNUM <- precreate_realdata(folder = pa, dnum = TRUE)

## Create the function for Separate
out_add_DNUM <- create_realdata(dat = dat_pre_add_DNUM, savefilepath = "out_add_DNUM.rds")

out_add_DNUM <- readRDS("out_add_DNUM.rds")
out_add_DNUM <- out_add_DNUM %>% # readRDS("out_add_dnum.rds") %>%
  mutate(x11 = factor(x11)) %>%
  as_tibble()

set.seed(41)
allRES_add_DNUM <- pred_realdata(dat = out_add_DNUM, cov = 1:11, savefilepath = "allRES_add_DNUM_wVIMP.rds")

allRES_add_DNUM <- readRDS("allRES_add_DNUM_wVIMP.rds")

################# ============= VIMP results =============== ##############################
# Show exact values
allRES_add_DNUM$vimp_supphlgr %>%
  pivot_longer(-year, names_to = "var", values_to = "vimp") %>%
  mutate(vimp = as.numeric(vimp)) %>%
  group_by(year) %>%
  mutate(rank = rank(desc(vimp))) %>%
  pivot_wider(names_from = var, values_from = c(vimp, rank)) %>%
  select(year, starts_with("rank"))

# Plot the results
names(allRES_add_DNUM$vimp_supphlgr) <- c("year", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11", "t", "u")
allRES_add_DNUM$vimp_supphlgr %>%
  pivot_longer(-year, names_to = "var", values_to = "vimp") %>%
  mutate(vimp = as.numeric(vimp)) %>%
  group_by(var) %>%
  summarise(vimp = mean(vimp)) %>%
  arrange(vimp) %>%
  mutate(var = factor(var, levels = var)) %>%
  ggplot(mapping = aes(x = var, y = vimp)) +
  geom_bar(stat = "identity", width = 0.75, fill = "gray50") + 
  geom_text(aes(label = format(round(vimp,4), 4)), 
            position = position_dodge(width = 1),
            hjust = -0.15) +
  ylim(0, 0.00875) +
  coord_flip() + 
  labs(y = "VIMP", x = "Covariates") +
  theme_classic()
  


################# ============= BS results =============== ##############################
pred <- allRES_add_DNUM$pred
## compute Brier score
BS <- pred %>% 
  mutate(horizon3 = ifelse(u - t == 3, 1, 0)) %>%
  group_by(t, u, year, horizon3) %>%
  summarise(across(c("hsuppDTPO", "hsupphlgr", "hsupp0hlgr", "hSKHThlgr", "hsephlgr", "hbench"), #.names = "{.col}.{.fn}"))
                   ~ mean((.x - y)^2))) %>%
  group_by(horizon3) %>%
  summarise(across(c("hsuppDTPO", "hsupphlgr", "hsupp0hlgr", "hSKHThlgr", "hsephlgr", "hbench"), #.names = "{.col}.{.fn}"))
                   ~ mean(.x)))

print(BS)
## plot cumulative gain chart
plotgainchartmany(dat = as.data.frame(pred), ngroups = 100, horizon3 = 1)
plotgainchartmany(dat = as.data.frame(pred), ngroups = 100, horizon3 = 0)

################# ========= Compute area under the curve ============= ####################
AUGC <- lapply(c(0, 1), function(horizon3){
  computeAUGC(dat = as.data.frame(pred), ngroups = 1000, horizon3 = horizon3) 
}) %>% 
  bind_rows()
print(AUGC)

