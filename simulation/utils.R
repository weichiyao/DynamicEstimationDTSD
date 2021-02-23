###### ======== autocorrelation: matsigma in traindtv_autocorr_gnrt.R ======== ######
create_matsigma <- function(scenario, autocorrtype){
  if (scenario == "0TI2TV") {
    if (autocorrtype == "strong"){
      matsigma <- .7 * diag(6) + matrix(.2, 6, 6)
      # to make the first two variables time-independent
      matsigma[, c(1:3, 5)] <- 0
      matsigma[c(1:3, 5), ] <- 0
    } else {
      matsigma <- .3 * diag(6) + matrix(.1, 6, 6)
      # to make the first two variables time-independent
      matsigma[, c(1:3, 5)] <- 0
      matsigma[c(1:3, 5), ] <- 0
    }
    
  } else if (scenario == "0TI4TV"){
    if (autocorrtype == "strong"){
      matsigma <- .7 * diag(6) + matrix(.2, 6, 6)
      # to make the first two variables time-independent
      matsigma[1, ] <- c(1, 0, 0, 0, 0, 0)
      matsigma[2, ] <- c(0, 1, 0, 0, 0, 0)
      matsigma[3:6, 1:2] <- 0
    } else {
      matsigma <- .3 * diag(6) + matrix(.1, 6, 6)
      # to make the first two variables time-independent
      matsigma[1, ] <- c(1, 0, 0, 0, 0, 0)
      matsigma[2, ] <- c(0, 1, 0, 0, 0, 0)
      matsigma[3:6, 1:2] <- 0
    }
  } else if (scenario == "1TI4TV"){
    if (autocorrtype == "strong"){
      matsigma <- .7 * diag(6) + matrix(.2, 6, 6)
      # to make the first two variables time-independent
      matsigma[1, ] <- c(1, 0, 0, 0, 0, 0)
      matsigma[2, ] <- c(0, 1, 0, 0, 0, 0)
      matsigma[3:6, 1] <- 0
    } else {
      matsigma <- .3 * diag(6) + matrix(.1, 6, 6)
      # to make the first two variables time-independent
      matsigma[1, ] <- c(1, 0, 0, 0, 0, 0)
      matsigma[2, ] <- c(0, 1, 0, 0, 0, 0)
      matsigma[3:6, 1] <- 0
    }
  } else if (scenario == "2TI1TV"){
    if (autocorrtype == "strong"){
      matsigma <- diag(6)
      matsigma[6, ] <- c(0.2, 0.2, 0, 0, 0, 0.9)
      
    } else {
      matsigma <- diag(6)
      matsigma[6, ] <- c(0.1, 0.1, 0, 0, 0, 0.4)
    }
  } else if (scenario == "2TI4TV"){
    if (autocorrtype == "strong"){
      matsigma <- .7 * diag(6) + matrix(.2, 6, 6)
      # to make the first two variables time-independent
      matsigma[1, ] <- c(1, 0, 0, 0, 0, 0)
      matsigma[2, ] <- c(0, 1, 0, 0, 0, 0)
    } else {
      matsigma <- .3 * diag(6) + matrix(.1, 6, 6)
      # to make the first two variables time-independent
      matsigma[1, ] <- c(1, 0, 0, 0, 0, 0)
      matsigma[2, ] <- c(0, 1, 0, 0, 0, 0)
    }
  } 
  return(matsigma)
}


###### ======== censoring rate: Censor.time in traindtv_autocorr_gnrt.R ======== ######
create_ctime <- function(censor.rate, nsub, scenario, SNR, model, distribution, nperiod){
  if (nperiod == 8){
    if (SNR == "high"){
      if (model == "linear"){
        if (censor.rate == "50%") {
          if (distribution == "Exp") {
            Censor.time = rexp(nsub, rate = 1/200)
          } else if (distribution == "WI") {
            Censor.time = rexp(nsub, rate = 1/150) 
          } else if (distribution == "Gtz") {
            Censor.time = rexp(nsub, rate = 1/20) 
          }
        } else if (censor.rate == "10%"){
          Censor.time <- rep(Inf, nsub)
        } 
      } else if (model == "nonlinear"){
        if (censor.rate == "50%") {
          if (distribution == "Exp") {
            Censor.time = rexp(nsub, rate = 1/7.2)
          } else if (distribution == "WI") {
            Censor.time = rexp(nsub, rate = 1/27)
          } else if (distribution == "Gtz") {
            Censor.time = rexp(nsub, rate = 1/200)
          }
        } else if (censor.rate == "10%"){
          Censor.time <- rep(Inf, nsub)
        }
      } else if (model == "interaction") {
        if (censor.rate == "50%") {
          if (distribution == "Exp") {
            Censor.time = rexp(nsub, rate = 1/120)
          } else if (distribution == "WI") {
            Censor.time = rexp(nsub, rate = 1/54)
          } else if (distribution == "Gtz") {
            Censor.time = rexp(nsub, rate = 1/29)
          }
        } else if (censor.rate == "10%"){
          Censor.time <- rep(Inf, nsub)
        } 
      } 
    } else if (SNR == "low") { # not boosted
      if (model == "linear"){
        if (censor.rate == "50%") {
          if (distribution == "Exp") {
            Censor.time = rexp(nsub, rate = 1/500) 
          } else if (distribution == "WI") {
            Censor.time = rexp(nsub, rate = 1/53)
          } else if (distribution == "Gtz") {
            Censor.time = rexp(nsub, rate = 1/10) / 16
          }
        } else if (censor.rate == "10%"){
          Censor.time <- rep(Inf, nsub)
        } 
      } else if (model == "nonlinear"){
        if (censor.rate == "50%") {
          if (distribution == "Exp") {
            Censor.time = rexp(nsub, rate = 1/27) 
          } else if (distribution == "WI") {
            Censor.time = rexp(nsub, rate = 1/212)
          } else if (distribution == "Gtz") {
            Censor.time = rexp(nsub, rate = 1/806)
          }
        } else if (censor.rate == "10%"){
          Censor.time <- rep(Inf, nsub)
        } 
      } else if (model == "interaction") {
        if (censor.rate == "50%") {
          if (distribution == "Exp") {
            Censor.time = rexp(nsub, rate = 1/295)
          } else if (distribution == "WI") {
            Censor.time = rexp(nsub, rate = 1/22)
          } else if (distribution == "Gtz") {
            Censor.time = rexp(nsub, rate = 1/65)
          }
        } else if (censor.rate == "10%"){
          Censor.time <- rep(Inf, nsub)
        } 
      } 
    }
  } else if (nperiod == 4) { # when nperiod != 8 we need to see other parameters
    if (scenario == "0TI2TV"){
      if (SNR == "high"){
        if (model == "linear"){
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/165)
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/225) 
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/32) 
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } else if (model == "nonlinear"){
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/16)
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/90)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/286)
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } else if (model == "interaction") {
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/570)
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/100)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/46)
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          }
        } 
      } else if (SNR == "low") { # not boosted
        if (model == "linear"){
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/585) 
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/53)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/10) / 14
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } else if (model == "nonlinear"){
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/28) 
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/155)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/780)
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } else if (model == "interaction") {
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/170)
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/22)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/72)
            } # end of distribution
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } # end of censoring
        }  # end of model 
      } # end of boosted 
    } else if (scenario == "0TI4TV") {
      if (SNR == "high"){
        if (model == "linear"){
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/150)
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/155) 
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/21) 
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } else if (model == "nonlinear"){
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/11)
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/32.5)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/195)
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } else if (model == "interaction") {
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/220)
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/79)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/46)
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } 
      } else if (SNR == "low") { # not boosted
        if (model == "linear"){
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/550) 
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/49)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/10) / 15
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } else if (model == "nonlinear"){
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/28) 
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/170)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/830)
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } else if (model == "interaction") {
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/120)
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/20)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/65)
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } 
      }
    } else if (scenario == "1TI4TV") {
      if (SNR == "high"){
        if (model == "linear"){
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/240)
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/170) 
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/24.5) 
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } else if (model == "nonlinear"){
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/15.2)
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/38)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/226)
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } else if (model == "interaction") {
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/255)
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/74)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/42)
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } 
      } else if (SNR == "low") { # not boosted
        if (model == "linear"){
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/605) 
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/53)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/10) / 14
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          }
        } else if (model == "nonlinear"){
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/36) 
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/192)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/880)
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } else if (model == "interaction") {
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/250)
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/21)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/75)
            } # end of distribution
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } # end of censoring
        }  # end of model 
      } # end of boosted 
    } else if (scenario == "2TI1TV") {
      if (SNR == "high"){
        if (model == "linear"){
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/260)
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/250) 
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/33) 
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } else if (model == "nonlinear"){
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/50)
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/65)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/242)
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } else if (model == "interaction") {
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/390)
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/77)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/44)
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } 
      } else if (SNR == "low") {
        if (model == "linear"){
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/605) 
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/57)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/10) / 14
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          }
        } else if (model == "nonlinear"){
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/45) 
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/182)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/820)
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } else if (model == "interaction") {
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/205)
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/19)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/69)
            } # end of distribution
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } # end of censoring
        } # end of model 
      } # end of boosted 
    } else if (scenario == "2TI4TV") {
      if (SNR == "high"){
        if (model == "linear"){
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/2)/50
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/200) 
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/27.5) 
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } else if (model == "nonlinear"){
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/10)
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/32)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/210)
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } else if (model == "interaction") {
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/190)
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/78)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/39)
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } 
      } else if (SNR == "low") { # not boosted
        if (model == "linear"){
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/670) 
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/61)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/10) / 13
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } else if (model == "nonlinear"){
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/33) 
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/222)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/910)
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } else if (model == "interaction") {
          if (censor.rate == "50%") {
            if (distribution == "Exp") {
              Censor.time = rexp(nsub, rate = 1/295)
            } else if (distribution == "WI") {
              Censor.time = rexp(nsub, rate = 1/28)
            } else if (distribution == "Gtz") {
              Censor.time = rexp(nsub, rate = 1/78)
            }
          } else if (censor.rate == "10%"){
            Censor.time <- rep(Inf, nsub)
          } 
        } 
      }
    } # end of scenarios
  } # end of nperiods
  return(Censor.time)
}

###### ======== create interaction model: itct_term in Timevarying_gnrt.R ======== ######
# INPUT #
# data is of size (nsub x nperiod) by ncov
# beta is a vector of length ncov 
#
# OUTPUT #
# a vector of length (nsub x nperiod)
itct_term <- function(data, coeff, scenario, SNR){
  if (scenario == "0TI4TV") {
    if (SNR == "high") {
      R3 <- 3 * log(abs(data %*% coeff$Beta2 + coeff$Beta0[3])) + coeff$Beta0[2]
      R1 <- 3 * (data[, 3] - 1) * (data[, 5] / 3) ^ 2 - 0.5 * data[, 4] ^ 0.5
      R2 <- data %*% coeff$Beta1 + coeff$Beta0[1]
      R0 <- (data[, 4] >= 0.5)  * R1 + 
        (data[, 4] < 0.5) * (data[, 6] >= 1) * R2 + 
        (data[, 4] < 0.5) * (data[, 6] < 1) * R3
      
    } else if (SNR == "low") {
      R3 <- 1.5 * log(abs(data %*% coeff$Beta2 + coeff$Beta0[3])) + coeff$Beta0[2]
      R1 <- (data[, 3] - 1) * (data[, 5] / 3) ^ 2 - 0.5 * data[, 4] ^ 0.5
      R2 <- data %*% coeff$Beta1 + coeff$Beta0[1]
      R0 <- (data[, 4] >= 0.5)  * R1 + 
        (data[, 4] < 0.5) * (data[, 6] >= 1) * R2 + 
        (data[, 4] < 0.5) * (data[, 6] < 1) * R3
    }

  } else if (scenario == "0TI2TV") {
    if (SNR == "high"){
      R3 <- 3 * log(abs(data %*% coeff$Beta2 + coeff$Beta0[3])) + coeff$Beta0[2]
      R1 <- - 1.042 - 0.5 * data[, 4] ^ 0.5
      R2 <- data %*% coeff$Beta1 + coeff$Beta0[1]
      R0 <- (data[, 4] >= 0.5)  * R1 + 
        (data[, 4] < 0.5) * (data[, 6] >= 1) * R2 + 
        (data[, 4] < 0.5) * (data[, 6] < 1) * R3
    } else if (SNR ==" low"){
      R3 <- 1.5 * log(abs(data %*% coeff$Beta2 + coeff$Beta0[3])) + coeff$Beta0[2]
      R1 <- - 0.347 - 0.5 * data[, 4] ^ 0.5
      R2 <- data %*% coeff$Beta1 + coeff$Beta0[1]
      R0 <- (data[, 4] >= 0.5)  * R1 + 
        (data[, 4] < 0.5) * (data[, 6] >= 1) * R2 + 
        (data[, 4] < 0.5) * (data[, 6] < 1) * R3
    }

  } else if (scenario == "1TI4TV") {
    if (SNR == "high"){
      R3 <- 3 * log(abs(data %*% coeff$Beta2 + coeff$Beta0[3])) + coeff$Beta0[2]
      R1 <- 3 * (data[, 3] - 1) * (data[, 5] / 3) ^ 2 - 0.5 * data[, 4] ^ data[, 2]
      R2 <- data %*% coeff$Beta1 + coeff$Beta0[1]
      R0 <- (data[, 2] >= 0.5)  * R1 + 
        (data[, 2] < 0.5) * (data[, 6] >= 1) * R2 + 
        (data[, 2] < 0.5) * (data[, 6] < 1) * R3

    } else if (SNR == "low"){
      R3 <- 1.5 * log(abs(data %*% coeff$Beta2 + coeff$Beta0[3])) + coeff$Beta0[2]
      R1 <- (data[, 3] - 1) * (data[, 5] / 3) ^ 2 - 0.5 * data[, 4] ^ data[, 2]
      R2 <- data %*% coeff$Beta1 + coeff$Beta0[1]
      R0 <- (data[, 2] >= 0.5)  * R1 + 
        (data[, 2] < 0.5) * (data[, 6] >= 1) * R2 + 
        (data[, 2] < 0.5) * (data[, 6] < 1) * R3
    }
    
  } else if (scenario == "2TI1TV") {
    if (SNR == "high"){
      R3 <- 3 * log(abs(data %*% coeff$Beta2 + coeff$Beta0[3])) + coeff$Beta0[2]
      R1 <- -1.042 + (data[, 1] - 1) * 0.5 ^ data[, 2]
      R2 <- data %*% coeff$Beta1 + coeff$Beta0[1]
      R0 <- (data[, 2] >= 0.5)  * R1 + 
        (data[, 2] < 0.5) * (data[, 6] >= 1) * R2 + 
        (data[, 2] < 0.5) * (data[, 6] < 1) * R3
    } else if (SNR == "low") {
      R3 <- 1.5 * log(abs(data %*% coeff$Beta2 + coeff$Beta0[3])) + coeff$Beta0[2]
      R1 <- -0.347 + (data[, 1] - 1) * 0.5 ^ data[, 2]
      R2 <- data %*% coeff$Beta1 + coeff$Beta0[1]
      R0 <- (data[, 2] >= 0.5)  * R1 + 
        (data[, 2] < 0.5) * (data[, 6] >= 1) * R2 + 
        (data[, 2] < 0.5) * (data[, 6] < 1) * R3
    }
    
  } else if (scenario == "2TI4TV") {
    if (SNR == "high"){
      R3 <- 3 * log(abs(data %*% coeff$Beta2)) + coeff$Beta0[2]
      R1 <- 3 * (data[, 3] - 1) * (data[, 5] / 3) ^ 2 + (data[, 1] - 1) * data[, 4] ^ data[, 2]
      R2 <- data %*% coeff$Beta1 + coeff$Beta0[1]
      R0 <- (data[, 2] >= 0.5)  * R1 + 
        (data[, 2] < 0.5) * (data[, 6] >= 1) * R2 + 
        (data[, 2] < 0.5) * (data[, 6] < 1) * R3
       
    } else if (SNR == "low"){
      R3 <- 1.5 * log(abs(data %*% coeff$Beta2)) + coeff$Beta0 
      R1 <- (data[, 3] - 1) * (data[, 5] / 3) ^ 2 + (data[, 1] - 1) * data[, 4] ^ data[, 2]
      R2 <- data %*% coeff$Beta1 + coeff$Beta0
      R0 <- (data[, 2] >= 0.5)  * R1 + 
        (data[, 2] < 0.5) * (data[, 6] >= 1) * R2 + 
        (data[, 2] < 0.5) * (data[, 6] < 1) * R3
    }
  }
  
  return(R0)
}

###### ======== create coefficient lists: Coeff in Timevarying_gnrt.R ======== ######
create_coeff <- function(nsub, scenario, model, distribution, SNR, nperiod){
  if (nperiod == 8) {
    if (SNR == "high") {
      Beta1 <- 3 * c(1, -1, 1, -1, 0.25, -0.5) 
      Beta2 <- c(2, 2, 2, 2, 0.2, 0.2) 
      if (model == "linear") {
        if(distribution == "Exp"){
          Beta1 <- 3 * c(1, -1, 1, -1, -0.25, 0.5) 
          Lambda = 2
          V = 0
          Alpha = 0
          Beta0 = -5
          TS <- as.vector(replicate(nsub, 
                                    c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 0.1, shape1 = 2, shape2 = 10)) * 900)
          ))
        } else if(distribution == "WI"){
          Lambda = 0.006
          V = 2
          Alpha = 0
          Beta0 = -5
          TS <- as.vector(replicate(nsub, 
                                    c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 2, shape1 = 2, shape2 = 2)) * 900)
          ))
        } else if (distribution == "Gtz"){
          Alpha <- 0.1
          Lambda <- 0.01
          V <- 0
          Beta0 <- 0
          TS <- as.vector(replicate(nsub, 
                                    c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.002, b = 5, shape1 = 0.5, shape2 = 2)) * 800)
          ))
        } else {
          stop("Wrong distribution is specified.")
        }
      } else if (model == "nonlinear") {      
        Beta2 <- 2 * c(-2, 2, 2, -2, 0.2, 0.2) 
        if(distribution == "Exp") {
          Beta2 <- 2 * c(2, -2, 2, 2, -0.2, -0.2) 
          Lambda <- 0.001
          V <- 0
          Alpha <- 0
          Beta0 <- c(2, 0)
          TS <- as.vector(replicate(nsub, 
                                    c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.2, shape2 = 2)) * 300)
          ))
        } else if (distribution == "WI") {
          V <- 1.8
          Lambda <- 10
          Alpha <- 0
          Beta0 <- c(-10, 0) 
          TS <- as.vector(replicate(nsub, 
                                    c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.005, b = 1, shape1 = 2, shape2 = 2)) * 300)
          ))
        } else if (distribution == "Gtz") {
          Alpha <- 0.05
          Lambda <- 0.001
          V <- 0
          Beta0 <- c(-5, 0)
          TS <- as.vector(replicate(nsub, 
                                    c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 2, shape1 = 5, shape2 = 10)) * 800)
          ))
        } else {
          stop("Wrong Distribution is specified.")
        }
      } else if (model == "interaction") {
        if(distribution == "Exp"){
          Beta2 <- c(2, -2, 2, 2, -0.5, -0.5) 
          Lambda = 0.03
          V = 0
          Alpha = 0
          Beta0 <- c(-3, -2, 0)
          TS <- as.vector(replicate(nsub, 
                                    c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.005, shape1 = 0.8, shape2 = 0.8)) * 900)
          ))
        } else if(distribution == "WI"){
          Beta2 <- c(2, -2, 2, 2, -0.2, 0.2) 
          Lambda = 0.005
          V = 1.8
          Alpha = 0
          Beta0 <- c(-3, -3, 0)
          TS <- as.vector(replicate(nsub, 
                                    c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.1, b = 2, shape1 = 0.1, shape2 = 0.5)) * 500) 
          ))
        } else if(distribution == "Gtz"){
          Beta2 <- c(2, -2, 2, 2, -0.2, 0.2) 
          Alpha = 0.1
          Lambda = 0.01
          V = 0
          Beta0 <- c(-3, 2, 0)
          TS <- as.vector(replicate(nsub, 
                                    c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.1, shape2 = 10)) * 900)
          ))
        } else {
          stop("Wrong distribution is specified.")
        }# end of distribution
      } else {
        stop("Wrong model is given")
      } # end of model
    } else if (SNR == "low") { # not boosted
      Beta1 <- c(1, 1, 1, 1, 0.25, 0.5) 
      Beta2 <- c(2, 2, 2, 2, 0.2, 0.2) 
      if (model == "linear") {
        if (distribution == "Exp"){ # GOOD
          Beta1 <- c(1, -1, 1, -1, -0.25, 0.5) 
          Lambda = 0.5
          Alpha = 0
          V = 0
          Beta0 = -5
          TS <- as.vector(replicate(nsub, 
                                    c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.0001, b = 1, shape1 = 0.1, shape2 = 2)) * 900)
          ))
        } else if(distribution == "WI"){
          Lambda = 0.006
          V = 2
          Alpha = 0
          Beta0 = -5
          TS <- as.vector(replicate(nsub, 
                                    c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.005, b = 1, shape1 = 0.01, shape2 = 2)) * 400)
          ))
        } else if(distribution == "Gtz"){
          Beta1 <- c(1, 1, -1, 1, -0.25, 0.5) 
          Alpha = 0.1
          Lambda = 1
          V = 0
          Beta0 = 0
          TS <- as.vector(replicate(nsub, 
                                    c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.006, b = 5, shape1 = 0.05, shape2 = 5))*80)
          ))
        } else {
          stop("Wrong Distribution is specified.")
        }
      } else if (model == "nonlinear") {
        if(distribution == "Exp") {
          Beta2 <- c(2, -2, 2, 2, -0.2, 0.2) 
          Lambda <- 0.03
          V <- 0.2
          Alpha <- 0
          Beta0 <- c(0, 0)
          TS <- as.vector(replicate(nsub, 
                                    c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.2, shape2 = 2)) * 300)
          ))
        } else if (distribution == "WI") {
          V <- 1.8
          Lambda <- 0.008
          Alpha <- 0
          Beta0 <- c(-6, 0)
          TS <- as.vector(replicate(nsub, 
                                    c(0, sort(rtrunc(nperiod - 1, a = 0.001, spec = "beta", shape1 = 0.9, shape2 = 4)) * 800)
          ))
        } else if (distribution == "Gtz") {
          Alpha <- 0.01
          Lambda <- 0.001
          V <- 0
          Beta0 <- c(-5, 0)
          TS <- as.vector(replicate(nsub, 
                                    c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.1, b = 2, shape1 = 0.5, shape2 = 0.5)) * 800)
          ))
        } else {
          stop("Wrong Distribution is specified.")
        }
      } else if (model == "interaction") {
        Beta0 <- c(-5, -5, 0)
        if(distribution == "Exp"){
          Beta1 <- c(1, -1, 1, -1, -0.25, 0.5)
          Beta2 <- c(2, -2, 2, 2, -0.2, 0.2) 
          Lambda = 0.14
          V = 0.5
          Alpha = 0
          TS <- as.vector(replicate(nsub, 
                                    c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.0001, b = 1, shape1 = 0.01, shape2 = 2)) * 900)
          ))
        } else if(distribution == "WI"){
          Beta1 <- c(1, -1, 1, 1, -0.25, 0.5)
          Beta2 <- c(2, -2, 2, 2, -0.2, 0.2) 
          Beta0 <- c(-4, -4, 0)
          Lambda = 0.08
          V = 1.8
          Alpha = 0
          TS <- as.vector(replicate(nsub, 
                                    c(0, sort(rtrunc(nperiod - 1, spec = "beta",b = 0.8, shape1 = 0.6, shape2 = 1)) * 500) # * 800 
          ))
        } else if(distribution == "Gtz"){
          Beta1 <- c(1, -1, 1, 1, -0.25, 0.5)
          Beta2 <- c(2, -2, 2, 2, -0.2, 0.2) 
          Alpha = 0.1
          Lambda = 0.015
          V = 0
          TS <- as.vector(replicate(nsub, 
                                    c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.2, shape2 = 2)) * 900)
          ))
        } else{
          stop("Wrong Distribution is specified.")
        } # end of distribution
      } else {
        stop("Wrong model is given")
      }# end of model
    } # end of boosted
  } else if (nperiod == 4) {
    if (scenario == "0TI2TV") {
      if (SNR == "high") {
        Beta1 <- 3 * c(0, 0, 0, -1, 0, -0.5) 
        Beta2 <- c(0, 0, 0, 2, 0, 0.2) 
        if (model == "linear") {
          if(distribution == "Exp"){
            Beta1 <- 3 * c(0, 0, 0, -1, 0, 0.5) 
            Lambda = 2
            V = 0
            Alpha = 0
            Beta0 = -5.375 # - 5;  3 / 2 - 3 * 5 / 8 = -0.375
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 0.1, shape1 = 2, shape2 = 10)) * 900)
            ))
          } else if(distribution == "WI"){
            Lambda = 0.006
            V = 2
            Alpha = 0
            Beta0 = -1.625 # - 5; 3 / 2 + 15 / 8 = 3.375
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 2, shape1 = 2, shape2 = 2)) * 900)
            ))
          } else if (distribution == "Gtz"){
            Alpha <- 0.1
            Lambda <- 0.01
            V <- 0
            Beta0 <- 3.375 # 0; 3 / 2 + 15 / 8 = 3.375
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.002, b = 5, shape1 = 0.5, shape2 = 2)) * 800)
            ))
          } else {
            stop("Wrong distribution is specified.")
          }
        } else if (model == "nonlinear") {      
          Beta2 <- 2 * c(0, 0, 0, -2, 0, 0.2) 
          if(distribution == "Exp") {
            Beta2 <- 2 * c(0, 0, 0, 2, 0, -0.2) 
            Lambda <- 0.001
            V <- 0
            Alpha <- 0
            Beta0 <- c(2, 1) # 4 / 2 - 2.5 * 0.2 * 2 = 2 - 1 = 1 
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.2, shape2 = 2)) * 300)
            ))
          } else if (distribution == "WI") {
            V <- 1.8
            Lambda <- 10
            Alpha <- 0
            Beta0 <- c(-10, 3) # 4 / 2 + 2.5 * 0.2 * 2 = 2 + 1 = 3 
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.005, b = 1, shape1 = 2, shape2 = 2)) * 300)
            ))
          } else if (distribution == "Gtz") {
            Alpha <- 0.05
            Lambda <- 0.001
            V <- 0
            Beta0 <- c(-5, 3) # 4 / 2 + 2.5 * 0.2 * 2 = 2 + 1 = 3 
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 2, shape1 = 5, shape2 = 10)) * 800)
            ))
          } else {
            stop("Wrong Distribution is specified.")
          }
        } else if (model == "interaction") {
          if(distribution == "Exp"){
            Beta2 <- c(0, 0, 0, 2, 0, -0.5) 
            Lambda = 0.03
            V = 0
            Alpha = 0
            Beta0 <- c(0.375, -2, -0.25) # -3; 3 / 2 + 15 / 8 = 3.375   2 / 2 - 2.5 * 0.5 = 1 - 1.25 = -0.25
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.005, shape1 = 0.8, shape2 = 0.8)) * 900)
            ))
          } else if(distribution == "WI"){
            Beta2 <- c(0, 0, 0, 2, 0, 0.2) 
            Lambda = 0.005
            V = 1.8
            Alpha = 0
            Beta0 <- c(0.375, -3, 0.5) # -3; 3 / 2 + 15 / 8 = 3.375   2 / 2 - 2.5 * 0.2 = 1 - 0.5 = 0.5
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.1, b = 2, shape1 = 0.1, shape2 = 0.5)) * 500) 
            ))
          } else if(distribution == "Gtz"){
            Beta2 <- c(0, 0, 0, 2, 0, 0.2) 
            Alpha = 0.1
            Lambda = 0.01
            V = 0
            Beta0 <- c(0.375, 2, 0.5) # -3; 3 / 2 + 15 / 8 = 3.375   2 / 2 - 2.5 * 0.2 = 1 - 0.5 = 0.5
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.1, shape2 = 10)) * 900)
            ))
          } else {
            stop("Wrong distribution is specified.")
          }# end of distribution
        } else {
          stop("Wrong model is given")
        } # end of model
      } else if (SNR == "low") { # not boosted
        Beta1 <- c(0, 0, 0, 1, 0, 0.5) 
        Beta2 <- c(0, 0, 0, 2, 0, 0.2) 
        if (model == "linear") {
          if (distribution == "Exp"){ # GOOD
            Beta1 <- c(0, 0, 0, -1, 0, 0.5) 
            Lambda = 0.5
            Alpha = 0
            V = 0
            Beta0 = -5.125 # -5; 0.5 - 5 / 8
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.0001, b = 1, shape1 = 0.1, shape2 = 2)) * 900)
            ))
          } else if(distribution == "WI"){
            Lambda = 0.006
            V = 2
            Alpha = 0
            Beta0 = -2.875 # -5; 1.5 + 5 / 8
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.005, b = 1, shape1 = 0.01, shape2 = 2)) * 400)
            ))
          } else if(distribution == "Gtz"){
            Beta1 <- c(0, 0, 0, 1, 0, 0.5) 
            Alpha = 0.1
            Lambda = 1
            V = 0
            Beta0 = -0.125 # 0; 0.5 - 5 / 8
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.006, b = 5, shape1 = 0.05, shape2 = 5))*80)
            ))
          } else {
            stop("Wrong Distribution is specified.")
          }
        } else if (model == "nonlinear") {
          if(distribution == "Exp") {
            Beta2 <- c(0, 0, 0, 2, 0, 0.2) 
            Lambda <- 0.03
            V <- 0.2
            Alpha <- 0
            Beta0 <- c(0, 0.5) # 2 / 2 - 0.2 * 2.5 = 1 - 0.5
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.2, shape2 = 2)) * 300)
            ))
          } else if (distribution == "WI") {
            V <- 1.8
            Lambda <- 0.008
            Alpha <- 0
            Beta0 <- c(-6, 4.5) # 8 / 2 + 0.2 * 2.5 = 4 + 0.5 = 4.5
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, a = 0.001, spec = "beta", shape1 = 0.9, shape2 = 4)) * 800)
            ))
          } else if (distribution == "Gtz") {
            Alpha <- 0.01
            Lambda <- 0.001
            V <- 0
            Beta0 <- c(-5, 4.5) # 8 / 2 + 0.2 * 2.5 = 4 + 0.5 = 4.5
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.1, b = 2, shape1 = 0.5, shape2 = 0.5)) * 800)
            ))
          } else {
            stop("Wrong Distribution is specified.")
          }
        } else if (model == "interaction") {
          if(distribution == "Exp"){
            Beta1 <- c(0, 0, 0, -1, 0, 0.5)
            Beta2 <- c(0, 0, 0, 2, 0, 0.2) 
            Lambda = 0.14
            V = 0.5
            Alpha = 0
            Beta0 <- c(-5.125, -5, 0.5) # -5; 1 / 2 - 5 / 8 = -0.125  2 / 2 - 0.2 * 2.5 = 1 - 0.5 = 0.5
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.0001, b = 1, shape1 = 0.01, shape2 = 2)) * 900)
            ))
          } else if(distribution == "WI"){
            Beta1 <- c(0, 0, 0, 1, 0, 0.5)
            Beta2 <- c(0, 0, 0, 2, 0, 0.2) 
            Beta0 <- c(-4.125, -4, 0.5) # -4; 1 / 2 - 5 / 8 = -0.125  2 / 2 - 0.2 * 2.5 = 1 - 0.5 = 0.5
            Lambda = 0.08
            V = 1.8
            Alpha = 0
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta",b = 0.8, shape1 = 0.6, shape2 = 1)) * 500) # * 800 
            ))
          } else if(distribution == "Gtz"){
            Beta1 <- c(0, 0, 0, 1, 0, 0.5)
            Beta2 <- c(0, 0, 0, 2, 0, 0.2) 
            Alpha = 0.1
            Lambda = 0.015
            V = 0
            Beta0 <- c(-5.125, -5, 0.5) # -5; 1 / 2 - 5 / 8 = -0.125  2 / 2 - 0.2 * 2.5 = 1 - 0.5 = 0.5
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.2, shape2 = 2)) * 900)
            ))
          } else{
            stop("Wrong Distribution is specified.")
          } # end of distribution
        } else {
          stop("Wrong model is given")
        }# end of model
      } # end of boosted
    } else if (scenario == "0TI4TV") {
      if (SNR == "high") {
        Beta1 <- 3 * c(0, 0, 1, -1, 0.25, -0.5) 
        Beta2 <- c(0, 0, 2, 2, 0.2, 0.2) 
        if (model == "linear") {
          if(distribution == "Exp"){
            Beta1 <- 3 * c(0, 0, 1, -1, -0.25, 0.5) 
            Lambda = 2
            V = 0
            Alpha = 0
            Beta0 = -5
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 0.1, shape1 = 2, shape2 = 10)) * 900)
            ))
          } else if(distribution == "WI"){
            Lambda = 0.006
            V = 2
            Alpha = 0
            Beta0 = -5
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 2, shape1 = 2, shape2 = 2)) * 900)
            ))
          } else if (distribution == "Gtz"){
            Alpha <- 0.1
            Lambda <- 0.01
            V <- 0
            Beta0 <- 0
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.002, b = 5, shape1 = 0.5, shape2 = 2)) * 800)
            ))
          } else {
            stop("Wrong distribution is specified.")
          }
        } else if (model == "nonlinear") {      
          Beta2 <- 2 * c(0, 0, 2, -2, 0.2, 0.2) 
          if(distribution == "Exp") {
            Beta2 <- 2 * c(0, 0, 2, 2, -0.2, -0.2) 
            Lambda <- 0.001
            V <- 0
            Alpha <- 0
            Beta0 <- c(2, 0)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.2, shape2 = 2)) * 300)
            ))
          } else if (distribution == "WI") {
            V <- 1.8
            Lambda <- 10
            Alpha <- 0
            Beta0 <- c(-10, 0)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.005, b = 1, shape1 = 2, shape2 = 2)) * 300)
            ))
          } else if (distribution == "Gtz") {
            Alpha <- 0.05
            Lambda <- 0.001
            V <- 0
            Beta0 <- c(-5, 0)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 2, shape1 = 5, shape2 = 10)) * 800)
            ))
          } else {
            stop("Wrong Distribution is specified.")
          }
        } else if (model == "interaction") {
          if(distribution == "Exp"){
            Beta2 <- c(0, 0, 2, 2, -0.5, -0.5) 
            Lambda = 0.03
            V = 0
            Alpha = 0
            Beta0 <- c(-3, -2, 0)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.005, shape1 = 0.8, shape2 = 0.8)) * 900)
            ))
          } else if(distribution == "WI"){
            Beta2 <- c(0, 0, 2, 2, -0.2, 0.2) 
            Lambda = 0.005
            V = 1.8
            Alpha = 0
            Beta0 <- c(-3, -3, 0)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.1, b = 2, shape1 = 0.1, shape2 = 0.5)) * 500) 
            ))
          } else if(distribution == "Gtz"){
            Beta2 <- c(0, 0, 2, 2, -0.2, 0.2) 
            Alpha = 0.1
            Lambda = 0.01
            V = 0
            Beta0 <- c(-3, 2, 0)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.1, shape2 = 10)) * 900)
            ))
          } else {
            stop("Wrong distribution is specified.")
          }# end of distribution
        } else {
          stop("Wrong model is given")
        } # end of model
      } else if (SNR == "low") { # not boosted
        Beta1 <- c(0, 0, 1, 1, 0.25, 0.5) 
        Beta2 <- c(0, 0, 2, 2, 0.2, 0.2) 
        if (model == "linear") {
          if (distribution == "Exp"){ # GOOD
            Beta1 <- c(0, 0, 1, -1, -0.25, 0.5) 
            Lambda = 0.5
            Alpha = 0
            V = 0
            Beta0 = -5
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.0001, b = 1, shape1 = 0.1, shape2 = 2)) * 900)
            ))
          } else if(distribution == "WI"){
            Lambda = 0.006
            V = 2
            Alpha = 0
            Beta0 = -4
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.005, b = 1, shape1 = 0.01, shape2 = 2)) * 400)
            ))
          } else if(distribution == "Gtz"){
            Beta1 <- c(0, 0, -1, 1, -0.25, 0.5) 
            Alpha = 0.1
            Lambda = 1
            V = 0
            Beta0 = 1
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.006, b = 5, shape1 = 0.05, shape2 = 5))*80)
            ))
          } else {
            stop("Wrong Distribution is specified.")
          }
        } else if (model == "nonlinear") {
          if(distribution == "Exp") {
            Beta2 <- c(0, 0, 2, 2, -0.2, 0.2) 
            Lambda <- 0.03
            V <- 0.2
            Alpha <- 0
            Beta0 <- c(0, 0)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.2, shape2 = 2)) * 300)
            ))
          } else if (distribution == "WI") {
            V <- 1.8
            Lambda <- 0.008
            Alpha <- 0
            Beta0 <- c(-6, 2)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, a = 0.001, spec = "beta", shape1 = 0.9, shape2 = 4)) * 800)
            ))
          } else if (distribution == "Gtz") {
            Alpha <- 0.01
            Lambda <- 0.001
            V <- 0
            Beta0 <- c(-5, 2)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.1, b = 2, shape1 = 0.5, shape2 = 0.5)) * 800)
            ))
          } else {
            stop("Wrong Distribution is specified.")
          }
        } else if (model == "interaction") {
          
          if(distribution == "Exp"){
            Beta1 <- c(0, 0, 1, -1, -0.25, 0.5)
            Beta2 <- c(0, 0, 2, 2, -0.2, 0.2) 
            Lambda = 0.14
            V = 0.5
            Alpha = 0
            Beta0 <- c(-5, -5, 0)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.0001, b = 1, shape1 = 0.01, shape2 = 2)) * 900)
            ))
          } else if(distribution == "WI"){
            Beta1 <- c(0, 0, 1, 1, -0.25, 0.5)
            Beta2 <- c(0, 0, 2, 2, -0.2, 0.2) 
            Beta0 <- c(-4, -4, 0)
            Lambda = 0.08
            V = 1.8
            Alpha = 0
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta",b = 0.8, shape1 = 0.6, shape2 = 1)) * 500) # * 800 
            ))
          } else if(distribution == "Gtz"){
            Beta1 <- c(0, 0, 1, 1, -0.25, 0.5)
            Beta2 <- c(0, 0, 2, 2, -0.2, 0.2) 
            Alpha = 0.1
            Lambda = 0.015
            V = 0
            Beta0 <- c(-5, -5, 0)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.2, shape2 = 2)) * 900)
            ))
          } else{
            stop("Wrong Distribution is specified.")
          } # end of distribution
        } else {
          stop("Wrong model is given")
        }# end of model
      } # end of boosted
    } else if (scenario == "1TI4TV") {
      if (SNR == "high") {
        Beta1 <- 3 * c(0, -1, 1, -1, 0.25, -0.5) 
        Beta2 <- c(2, 2, 2, 2, 0.2, 0.2) 
        if (model == "linear") {
          if(distribution == "Exp"){
            Beta1 <- 3 * c(0, -1, 1, -1, -0.25, 0.5) 
            Lambda = 2
            V = 0
            Alpha = 0
            Beta0 = -3.5
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 0.1, shape1 = 2, shape2 = 10)) * 900)
            ))
          } else if(distribution == "WI"){
            Lambda = 0.006
            V = 2
            Alpha = 0
            Beta0 = -3.5
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 2, shape1 = 2, shape2 = 2)) * 900)
            ))
          } else if (distribution == "Gtz"){
            Alpha <- 0.1
            Lambda <- 0.01
            V <- 0
            Beta0 <- 1.5
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.002, b = 5, shape1 = 0.5, shape2 = 2)) * 800)
            ))
          } else {
            stop("Wrong distribution is specified.")
          }
        } else if (model == "nonlinear") {      
          Beta2 <- 2 * c(0, 2, 2, -2, 0.2, 0.2) 
          if(distribution == "Exp") {
            Beta2 <- 2 * c(0, -2, 2, 2, -0.2, -0.2) 
            Lambda <- 0.001
            V <- 0
            Alpha <- 0
            Beta0 <- c(2, 2)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.2, shape2 = 2)) * 300)
            ))
          } else if (distribution == "WI") {
            V <- 1.8
            Lambda <- 10
            Alpha <- 0
            Beta0 <- c(-10, -2)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.005, b = 1, shape1 = 2, shape2 = 2)) * 300)
            ))
          } else if (distribution == "Gtz") {
            Alpha <- 0.05
            Lambda <- 0.001
            V <- 0
            Beta0 <- c(-5, -2)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 2, shape1 = 5, shape2 = 10)) * 800)
            ))
          } else {
            stop("Wrong Distribution is specified.")
          }
        } else if (model == "interaction") {
          if(distribution == "Exp"){
            Beta2 <- c(0, -2, 2, 2, -0.5, -0.5) 
            Lambda = 0.03
            V = 0
            Alpha = 0
            Beta0 <- c(-1.5, -2, 1)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.005, shape1 = 0.8, shape2 = 0.8)) * 900)
            ))
          } else if(distribution == "WI"){
            Beta2 <- c(0, -2, 2, 2, -0.2, 0.2) 
            Lambda = 0.005
            V = 1.8
            Alpha = 0
            Beta0 <- c(-1.5, -3, 1)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.1, b = 2, shape1 = 0.1, shape2 = 0.5)) * 500) 
            ))
          } else if(distribution == "Gtz"){
            Beta2 <- c(0, -2, 2, 2, -0.2, 0.2) 
            Alpha = 0.1
            Lambda = 0.01
            V = 0
            Beta0 <- c(-1.5, 2, 1)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.1, shape2 = 10)) * 900)
            ))
          } else {
            stop("Wrong distribution is specified.")
          }# end of distribution
        } else {
          stop("Wrong model is given")
        } # end of model
      } else if (SNR == "low") { # not boosted
        Beta1 <- c(0, 1, 1, 1, 0.25, 0.5) 
        Beta2 <- c(0, 2, 2, 2, 0.2, 0.2) 
        if (model == "linear") {
          if (distribution == "Exp"){ # GOOD
            Beta1 <- c(0, -1, 1, -1, -0.25, 0.5) 
            Lambda = 0.5
            Alpha = 0
            V = 0
            Beta0 = -4.5
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.0001, b = 1, shape1 = 0.1, shape2 = 2)) * 900)
            ))
          } else if(distribution == "WI"){
            Lambda = 0.006
            V = 2
            Alpha = 0
            Beta0 = -4.5
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.005, b = 1, shape1 = 0.01, shape2 = 2)) * 400)
            ))
          } else if(distribution == "Gtz"){
            Beta1 <- c(0, 1, -1, 1, -0.25, 0.5) 
            Alpha = 0.1
            Lambda = 1
            V = 0
            Beta0 = 0.5
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.006, b = 5, shape1 = 0.05, shape2 = 5))*80)
            ))
          } else {
            stop("Wrong Distribution is specified.")
          }
        } else if (model == "nonlinear") {
          if(distribution == "Exp") {
            Beta2 <- c(0, -2, 2, 2, -0.2, 0.2) 
            Lambda <- 0.03
            V <- 0.2
            Alpha <- 0
            Beta0 <- c(0, 1)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.2, shape2 = 2)) * 300)
            ))
          } else if (distribution == "WI") {
            V <- 1.8
            Lambda <- 0.008
            Alpha <- 0
            Beta0 <- c(-6, 1)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, a = 0.001, spec = "beta", shape1 = 0.9, shape2 = 4)) * 800)
            ))
          } else if (distribution == "Gtz") {
            Alpha <- 0.01
            Lambda <- 0.001
            V <- 0
            Beta0 <- c(-5, 1)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.1, b = 2, shape1 = 0.5, shape2 = 0.5)) * 800)
            ))
          } else {
            stop("Wrong Distribution is specified.")
          }
        } else if (model == "interaction") {
          if(distribution == "Exp"){
            Beta1 <- c(0, -1, 1, -1, -0.25, 0.5)
            Beta2 <- c(0, -2, 2, 2, -0.2, 0.2) 
            Lambda = 0.14
            V = 0.5
            Alpha = 0
            Beta0 <- c(-4.5, -5, 1)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 2, shape1 = 1, shape2 = 2)) * 900)
            ))
          } else if(distribution == "WI"){
            Beta1 <- c(0, -1, 1, 1, -0.25, 0.5)
            Beta2 <- c(0, -2, 2, 2, -0.2, 0.2) 
            Beta0 <- c(-3.5, -4, 1)
            Lambda = 0.08
            V = 1.8
            Alpha = 0
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta",b = 0.8, shape1 = 0.6, shape2 = 1)) * 500) # * 800 
            ))
          } else if(distribution == "Gtz"){
            Beta1 <- c(0, -1, 1, 1, -0.25, 0.5)
            Beta2 <- c(0, -2, 2, 2, -0.2, 0.2) 
            Alpha = 0.1
            Lambda = 0.015
            V = 0
            Beta0 <- c(-4.5, -5, 1)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.2, shape2 = 2)) * 900)
            ))
          } else{
            stop("Wrong Distribution is specified.")
          } # end of distribution
        } else {
          stop("Wrong model is given")
        }# end of model
      } # end of boosted
    } else if (scenario == "2TI1TV") {
      if (SNR == "high") {
        Beta1 <- 3 * c(1, -1, 0, 0, 0, -0.5) 
        Beta2 <- c(2, 2, 0, 0, 0, 0.2) 
        if (model == "linear") {
          if(distribution == "Exp"){
            Beta1 <- 3 * c(1, -1, 0, 0, 0, 0.5) 
            Lambda = 2
            V = 0
            Alpha = 0
            Beta0 = -6.875 # -5 - 3*0.625
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 0.1, shape1 = 2, shape2 = 10)) * 900)
            ))
          } else if(distribution == "WI"){
            Lambda = 0.006
            V = 2
            Alpha = 0
            Beta0 = -3.125 # -5 + 3*0.625
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 2, shape1 = 2, shape2 = 2)) * 900)
            ))
          } else if (distribution == "Gtz"){
            Alpha <- 0.1
            Lambda <- 0.01
            V <- 0
            Beta0 <- 1.875 # 0.625+0
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.002, b = 5, shape1 = 0.5, shape2 = 2)) * 800)
            ))
          } else {
            stop("Wrong distribution is specified.")
          }
        } else if (model == "nonlinear") {      
          Beta2 <- 2 * c(-2, 2, 0, 0, 0, 0.2) 
          if(distribution == "Exp") {
            Beta2 <- 2 * c(2, -2, 0, 0, 0, -0.2) 
            Lambda <- 0.001
            V <- 0
            Alpha <- 0
            Beta0 <- c(2, 1.5)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.2, shape2 = 2)) * 300)
            ))
          } else if (distribution == "WI") {
            V <- 1.8
            Lambda <- 10
            Alpha <- 0
            Beta0 <- c(-10, 0.5)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.005, b = 1, shape1 = 2, shape2 = 2)) * 300)
            ))
          } else if (distribution == "Gtz") {
            Alpha <- 0.05
            Lambda <- 0.001
            V <- 0
            Beta0 <- c(-5, 0.5)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 2, shape1 = 5, shape2 = 10)) * 800)
            ))
          } else {
            stop("Wrong Distribution is specified.")
          }
        } else if (model == "interaction") {
          if(distribution == "Exp"){
            Beta2 <- c(2, -2, 0, 0, 0, -0.5) 
            Lambda = 0.03
            V = 0
            Alpha = 0
            Beta0 <- c(-1.125, -2, 0.75)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.005, shape1 = 0.8, shape2 = 0.8)) * 900)
            ))
          } else if(distribution == "WI"){
            Beta2 <- c(2, -2, 0, 0, 0, 0.2) 
            Lambda = 0.005
            V = 1.8
            Alpha = 0
            Beta0 <- c(-1.125, -3, 1.5)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.1, b = 2, shape1 = 0.1, shape2 = 0.5)) * 500) 
            ))
          } else if(distribution == "Gtz"){
            Beta2 <- c(2, -2, 0, 0, 0, 0.2) 
            Alpha = 0.1
            Lambda = 0.01
            V = 0
            Beta0 <- c(-1.125, 2, 1.5)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.1, shape2 = 10)) * 900)
            ))
          } else {
            stop("Wrong distribution is specified.")
          }# end of distribution
        } else {
          stop("Wrong model is given")
        } # end of model
      } else if (SNR == "low") { # not boosted
        Beta1 <- c(1, 1, 0, 0, 0, 0.5) 
        Beta2 <- c(2, 2, 0, 0, 0, 0.2) 
        if (model == "linear") {
          if (distribution == "Exp"){ 
            Beta1 <- c(1, -1, 0, 0, 0, 0.5) 
            Lambda = 0.5
            Alpha = 0
            V = 0
            Beta0 = -5.625
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.0001, b = 1, shape1 = 0.1, shape2 = 2)) * 900)
            ))
          } else if(distribution == "WI"){
            Lambda = 0.006
            V = 2
            Alpha = 0
            Beta0 = -3.375
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.005, b = 1, shape1 = 0.01, shape2 = 2)) * 400)
            ))
          } else if(distribution == "Gtz"){
            Beta1 <- c(1, 1, 0, 0, 0, 0.5) 
            Alpha = 0.1
            Lambda = 1
            V = 0
            Beta0 = -0.625
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.006, b = 5, shape1 = 0.05, shape2 = 5))*80)
            ))
          } else {
            stop("Wrong Distribution is specified.")
          }
        } else if (model == "nonlinear") {
          if(distribution == "Exp") {
            Beta2 <- c(2, -2, 0, 0, 0, 0.2) 
            Lambda <- 0.03
            V <- 0.2
            Alpha <- 0
            Beta0 <- c(0, 1.5)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.2, shape2 = 2)) * 300)
            ))
          } else if (distribution == "WI") {
            V <- 1.8
            Lambda <- 0.008
            Alpha <- 0
            Beta0 <- c(-6, 2.5)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, a = 0.001, spec = "beta", shape1 = 0.9, shape2 = 4)) * 800)
            ))
          } else if (distribution == "Gtz") {
            Alpha <- 0.01
            Lambda <- 0.001
            V <- 0
            Beta0 <- c(-5, 2.5)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.1, b = 2, shape1 = 0.5, shape2 = 0.5)) * 800)
            ))
          } else {
            stop("Wrong Distribution is specified.")
          }
        } else if (model == "interaction") {
          if(distribution == "Exp"){
            Beta1 <- c(1, -1, 0, 0, 0, 0.5)
            Beta2 <- c(2, -2, 0, 0, 0, 0.2) 
            Lambda = 0.14
            V = 0
            Alpha = 0
            Beta0 <- c(-5.625, -5, 1.5)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.0001, b = 1, shape1 = 0.01, shape2 = 2)) * 900)
            ))
          } else if(distribution == "WI"){
            Beta1 <- c(1, -1, 0, 0, 0, 0.5)
            Beta2 <- c(2, -2, 0, 0, 0, 0.2) 
            Beta0 <- c(-3.625, -4, 1.5)
            Lambda = 0.08
            V = 1.8
            Alpha = 0
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta",b = 0.8, shape1 = 0.6, shape2 = 1)) * 500) # * 800 
            ))
          } else if(distribution == "Gtz"){
            Beta1 <- c(1, -1, 0, 0, 0, 0.5)
            Beta2 <- c(2, -2, 0, 0, 0, 0.2) 
            Alpha = 0.1
            Lambda = 0.015
            V = 0
            Beta0 <- c(-4.625, -5, 1.5)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.2, shape2 = 2)) * 900)
            ))
          } else{
            stop("Wrong Distribution is specified.")
          } # end of distribution
        } else {
          stop("Wrong model is given")
        }# end of model
      } # end of boosted
    } else if (scenario == "2TI4TV"){
      if (SNR == "high") {
        Beta1 <- 3 * c(1, -1, 1, -1, 0.25, -0.5) 
        Beta2 <- c(2, 2, 2, 2, 0.2, 0.2) 
        if (model == "linear") {
          if(distribution == "Exp"){
            Beta1 <- 3 * c(1, 1, -1, 1, 0.25, 0.5) 
            Lambda = 10
            V = 0
            Alpha = 0
            Beta0 = -5
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.2, shape2 = 0.2)) * 1000)
            ))
          } else if(distribution == "WI"){
            Lambda = 0.006
            V = 2
            Alpha = 0
            Beta0 = -5
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 2, shape1 = 2, shape2 = 2)) * 900)
            ))
          } else if (distribution == "Gtz"){
            Alpha <- 0.1
            Lambda <- 0.01
            V <- 0
            Beta0 <- 0
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.002, b = 5, shape1 = 0.5, shape2 = 2)) * 800)
            ))
          } else {
            stop("Wrong distribution is specified.")
          }
        } else if (model == "nonlinear") {      
          Beta2 <- 2 * c(-2, 2, 2, -2, 0.2, 0.2) 
          if(distribution == "Exp") {
            Beta2 <- 2 * c(2, -2, 2, 2, -0.2, -0.2) 
            Lambda <- 0.001
            V <- 0
            Alpha <- 0
            Beta0 <- c(2, 0)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.2, shape2 = 2)) * 300)
            ))
          } else if (distribution == "WI") {
            V <- 1.8
            Lambda <- 10
            Alpha <- 0
            Beta0 <- c(-10, 0)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.005, b = 1, shape1 = 2, shape2 = 2)) * 300)
            ))
          } else if (distribution == "Gtz") {
            Alpha <- 0.05
            Lambda <- 0.001
            V <- 0
            Beta0 <- c(-5, 0)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 2, shape1 = 5, shape2 = 10)) * 800)
            ))
          } else {
            stop("Wrong Distribution is specified.")
          }
        } else if (model == "interaction") {
          if(distribution == "Exp"){
            Beta2 <- c(2, -2, 2, 2, -0.5, -0.5) 
            Lambda = 0.03
            V = 0
            Alpha = 0
            Beta0 <- c(-3, -2, 0)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.005, shape1 = 0.8, shape2 = 0.8)) * 900)
            ))
          } else if(distribution == "WI"){
            Beta2 <- c(2, -2, 2, 2, -0.2, 0.2) 
            Lambda = 0.005
            V = 1.8
            Alpha = 0
            Beta0 <- c(-3, -3, 0)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.1, b = 2, shape1 = 0.1, shape2 = 0.5)) * 500) 
            ))
          } else if(distribution == "Gtz"){
            Beta2 <- c(2, -2, 2, 2, -0.2, 0.2) 
            Alpha = 0.1
            Lambda = 0.01
            V = 0
            Beta0 <- c(-3, 2, 0)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.1, shape2 = 10)) * 900)
            ))
          } else {
            stop("Wrong distribution is specified.")
          }# end of distribution
        } else {
          stop("Wrong model is given")
        } # end of model
      } else if (SNR == "low") { # not boosted
        Beta1 <- c(1, 1, 1, 1, 0.25, 0.5) 
        Beta2 <- c(2, 2, 2, 2, 0.2, 0.2) 
        if (model == "linear") {
          if (distribution == "Exp"){ # GOOD
            Beta1 <- c(1, -1, 1, -1, -0.25, 0.5) 
            Lambda = 0.5
            Alpha = 0
            V = 0
            Beta0 = -5
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.0001, b = 1, shape1 = 0.1, shape2 = 2)) * 900)
            ))
          } else if(distribution == "WI"){
            Lambda = 0.006
            V = 2
            Alpha = 0
            Beta0 = -5
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.005, b = 1, shape1 = 0.01, shape2 = 2)) * 400)
            ))
          } else if(distribution == "Gtz"){
            Beta1 <- c(1, 1, -1, 1, -0.25, 0.5) 
            Alpha = 0.1
            Lambda = 1
            V = 0
            Beta0 = 0
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.006, b = 5, shape1 = 0.05, shape2 = 5))*80)
            ))
          } else {
            stop("Wrong Distribution is specified.")
          }
        } else if (model == "nonlinear") {
          if(distribution == "Exp") {
            Beta2 <- c(2, -2, 2, 2, -0.2, 0.2) 
            Lambda <- 0.03
            V <- 0.2
            Alpha <- 0
            Beta0 <- c(0, 0)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.2, shape2 = 2)) * 300)
            ))
          } else if (distribution == "WI") {
            V <- 1.8
            Lambda <- 0.008
            Alpha <- 0
            Beta0 <- c(-6, 0)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, a = 0.001, spec = "beta", shape1 = 0.9, shape2 = 4)) * 800)
            ))
          } else if (distribution == "Gtz") {
            Alpha <- 0.01
            Lambda <- 0.001
            V <- 0
            Beta0 <- c(-5, 0)
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.1, b = 2, shape1 = 0.5, shape2 = 0.5)) * 800)
            ))
          } else {
            stop("Wrong Distribution is specified.")
          }
        } else if (model == "interaction") {
          Beta0 <- c(-5, -5, 0)
          if(distribution == "Exp"){
            Beta1 <- c(1, -1, 1, -1, -0.25, 0.5)
            Beta2 <- c(2, -2, 2, 2, -0.2, 0.2) 
            Lambda = 0.14
            V = 0.5
            Alpha = 0
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.0001, b = 1, shape1 = 0.01, shape2 = 2)) * 900)
            ))
          } else if(distribution == "WI"){
            Beta1 <- c(1, -1, 1, 1, -0.25, 0.5)
            Beta2 <- c(2, -2, 2, 2, -0.2, 0.2) 
            Beta0 <- c(-5, -4, 0)
            Lambda = 0.08
            V = 1.8
            Alpha = 0
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta",b = 0.8, shape1 = 0.6, shape2 = 1)) * 500) # * 800 
            ))
          } else if(distribution == "Gtz"){
            Beta1 <- c(1, -1, 1, 1, -0.25, 0.5)
            Beta2 <- c(2, -2, 2, 2, -0.2, 0.2) 
            Alpha = 0.1
            Lambda = 0.015
            V = 0
            TS <- as.vector(replicate(nsub, 
                                      c(0, sort(rtrunc(nperiod - 1, spec = "beta", a = 0.001, b = 1, shape1 = 0.2, shape2 = 2)) * 900)
            ))
          } else{
            stop("Wrong Distribution is specified.")
          } # end of distribution
        } else {
          stop("Wrong model is given")
        }# end of model
      } # end of boosted
    }
  }
  
  Coeff <- list(Lambda = Lambda, Alpha = Alpha, V = V, 
                Beta1 = Beta1, Beta2 = Beta2, Beta0 = Beta0)
  
  return(list(TS = TS, Coeff = Coeff))
}

