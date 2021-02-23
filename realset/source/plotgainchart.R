plotgainchart=function(dat, method, ngroups, horizon3)
{
  
  # dat=read.table(file="d:/11000347/Desktop/landmarking_dynamic_prediction/discrete_time/bankruptcy_hoora/Allresultsforr.txt",header=TRUE)   # portable
  
  # dat=read.table(file="d:/11000347/Desktop/landmarking_dynamic_prediction/discrete_time/bankruptcy_hoora/Allresultsforr.txt",header=TRUE)    # maison
  
  # dat=read.table(file="d:/11000347/Bureau/bankruptcy_hoora/Allresultsforr.txt",header=TRUE)
  # 
  # dat=data.frame(dat)
  # 
  # dat[,"average"]=(dat[,"superpp"]+dat[,"schmid"]+dat[,"frtrn"]+dat[,"tu"])/4
  
  out=seq(1:ngroups)/ngroups*100
  
  
  
  if(horizon3==1)
  {
    for(year in unique(dat[, "year"]))
    {
      for(t in unique(dat[, "t"]))
      {
        
        dat1 <- dat[dat[, "year"] == year & dat[,"t"] == t & dat[, "u"] == t + 3, ]
        out <- cbind(out, computecumgain(pred = dat1[, method], y=dat1[, "y"],ngroups = ngroups)[, 2])
        
        
      }}
  }
  
  if(horizon3 == 0)
  {
    for(year in 1:3)
    {
      for(t in 0:3)
      {
        if(t==0)
        {
          for(u in 4:6)
          {
            dat1=dat[dat[, "year"] == year & dat[,"t"] == t & dat[,"u"] == u, ]
            out=cbind(out,computecumgain(pred=dat1[,method],y=dat1[,"y"],ngroups=ngroups)[,2])
          }
        }
        if(t==1)
        {
          for(u in 5:6)
          {
            dat1=dat[dat[, "year"] == year & dat[,"t"] == t & dat[,"u"] == u, ]
            out=cbind(out,computecumgain(pred=dat1[, method],y = dat1[, "y"], ngroups = ngroups)[, 2])
          }
        }
        
        if(t==2)
        {
          for(u in 6:6)
          {
            dat1=dat[dat[,"year"]==year & dat[,"t"]==t & dat[,"u"]==u,]
            out=cbind(out,computecumgain(pred=dat1[,method],y=dat1[,"y"],ngroups=ngroups)[,2])
          }
        }
        
        
      }}
  }
  
  
  # print(out)
  
  out=cbind(out,apply(out[,-1],1,mean))
  
  out=out[,c(1,NCOL(out))]
  
  out=data.frame(out)
  
  names(out)=c("base",method)
  
  # print(out)
  out
  
  
}

plotgainchartmany=function(dat, ngroups, horizon3)
{
  out <- cbind(c(0, seq(1:ngroups) / ngroups * 100), 
               sapply(names(dat[, c(5,4,3,2)]), function(namk){
                 c(0, plotgainchart(dat, namk, ngroups, horizon3)[, namk])
                 })) %>% as.data.frame()
  
  names(out)[1] = "base"
  names(out)[names(out) == "hsuppDTPO"] <- "SuperppDTPO"
  names(out)[names(out) == "hsephlgr"] <- "Separate"
  names(out)[names(out) == "hSKHThlgr"] <- "Poolt"
  names(out)[names(out) == "hsupphlgr"] <- "Superpp"
  names(out)[names(out) == "hsupp0hlgr"] <- "Superpp0"
  # cumulative gains chart
  
  colorset <- c("red", "blue", "forestgreen", "chocolate1")
  # linetypeset <- c("dotted", "dashed", "dotdash", "F1")

  matplot(matrix(rep(out[, 1], ncol(out) - 1), ncol = ncol(out) - 1),
          out[, -1],
          type = "l",
          xlim = c(0, 100),
          ylim = c(0, 100),
          col  = colorset[1:(ncol(out) - 1)],
          lty  = 2:ncol(out), # linetypeset[1:(ncol(out) - 1)],# 
          lwd = 1.5,
          axes = FALSE,
          ylab = "% bankruptcies",
          xlab = "% firms targeted at risk")
  
  axis(side=1,at=seq(0, 100, 10))
  axis(side=2,at=seq(0, 100, 10))
  box()
  abline(h=0,col = "lightgray")
  abline(h=10,col = "lightgray")
  abline(h=20,col = "lightgray")
  abline(h=30,col = "lightgray")
  abline(h=40,col = "lightgray")
  abline(h=50,col = "lightgray")
  abline(h=60,col = "lightgray")
  abline(h=70,col = "lightgray")
  abline(h=80,col = "lightgray")
  abline(h=90,col = "lightgray")
  abline(h=100,col = "lightgray")
  abline(v=0,col = "lightgray")
  abline(v=10,col = "lightgray")
  abline(v=20,col = "lightgray")
  abline(v=30,col = "lightgray")
  abline(v=40,col = "lightgray")
  abline(v=50,col = "lightgray")
  abline(v=60,col = "lightgray")
  abline(v=70,col = "lightgray")
  abline(v=80,col = "lightgray")
  abline(v=90,col = "lightgray")
  abline(v=100,col = "lightgray")
  
  leg = legend(70, 30,
               legend = names(out)[-1],
               col    = colorset[1:(ncol(out) - 1)],
               lty    = 2:ncol(out), #linetypeset[1:(ncol(out) - 1)], # 2:ncol(out)) 
               plot = FALSE)
  # adjust as desired
  leftx <- leg$rect$left
  rightx <- (leg$rect$left + leg$rect$w) * 1
  topy <- leg$rect$top
  bottomy <- (leg$rect$top - leg$rect$h) * 1
  
  legend(x = c(leftx, rightx), y = c(topy, bottomy), 
         legend = names(out)[-1],
         col    = colorset[1:(ncol(out) - 1)],
         lty    = 2:ncol(out)) #linetypeset[1:(ncol(out) - 1)])# 2:ncol(out))
  
  abline(0, 1)
  
  
  # lift chart
  
  #matplot(cbind(out[,1],out[,1],out[,1],out[,1],out[,1]),out[,-1]/cbind(out[,1],out[,1],out[,1],out[,1],out[,1]),type="l",xlim=c(0,100),ylim=c(0,10),col=1:5)
  #legend(60,3,legend=c("tu","schmid","frtrn","superpp","average"),col=1:5,lty=1)
  #abline(1,0)
  
  # print(apply(out,2,mean))
  
  # out
  
  
  
}



# Compute the area under the cumulative gain curve

computeAUGC <- function(dat, ngroups, horizon3) {
  out <- cbind(sapply(names(dat[, 1:5]), function(namk){
                 c(0, plotgainchart(dat, namk, ngroups, horizon3)[, namk])
               })) 
  
  # Compute area under the curve
  out %>% 
    as_tibble() %>%
    summarise(across(.cols = everything(),
                     ~ sum((diff(.x) / 2 + .x[-length(.x)])) / ngroups)) %>%
    mutate(horizon3 = horizon3)
}



