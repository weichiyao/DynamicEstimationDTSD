plotroccurve=function(dat, method, ngroups, horizon3)
{
  
  # dat=read.table(file="d:/11000347/Desktop/landmarking_dynamic_prediction/discrete_time/bankruptcy_hoora/Allresultsforr.txt",header=TRUE)   # portable
  
  out1 = seq(1:(ngroups + 1))
  out2 = seq(1:(ngroups + 1))
  
  
  
  if(horizon3 == 1)
  {
    for(year in unique(dat[, "year"]))
    {
      for(t in unique(dat[, "t"]))
      {
        
        dat1 = dat[dat[, "year"] == year & dat[, "t"] == t & dat[, "u"] == t + 3, ]
        ro = computeroc(pred=dat1[,method], y=dat1[, "y"],ngroups = ngroups)
        out1 = cbind(out1,ro[, 1])
        out2 = cbind(out1,ro[, 2])
        
        
      }}
  }
  
  if(horizon3 == 0)
  {
    for(year in 1:3)
    {
      for(t in 0:3)
      {
        if(t == 0)
        {
          for(u in 4:6)
          {
            dat1=dat[dat[,"year"]==year & dat[,"t"]==t & dat[,"u"]==u,]
            ro=computeroc(pred=dat1[,method],y=dat1[,"y"],ngroups=ngroups)
            out1=cbind(out1,ro[,1])
            out2=cbind(out1,ro[,2])
          }
        }
        if(t == 1)
        {
          for(u in 5:6)
          {
            dat1=dat[dat[,"year"]==year & dat[,"t"]==t & dat[,"u"]==u,]
            ro=computeroc(pred=dat1[,method],y=dat1[,"y"],ngroups=ngroups)
            out1=cbind(out1,ro[,1])
            out2=cbind(out1,ro[,2])
          }
        }
        
        if(t == 2)
        {
          for(u in 6:6)
          {
            dat1=dat[dat[,"year"]==year & dat[,"t"]==t & dat[,"u"]==u,]
            ro=computeroc(pred=dat1[,method],y=dat1[,"y"],ngroups=ngroups)
            out1=cbind(out1,ro[,1])
            out2=cbind(out1,ro[,2])
          }
        }
        
        
      }}
  }
  
  
  
  #print(out1)
  #print(out2)
  
  out1=apply(out1[, -1], 1, mean)
  out2=apply(out2[, -1], 1, mean)
  
  out=cbind(out1,out2)
  out=data.frame(out)
  
  names(out)=c(method,method)
  
  out
  
  
}

plotroccurvemany=function(dat, ngroups, horizon3)
{
  
  out1=plotroccurve(dat, "hbench",ngroups,horizon3)
  out2=plotroccurve(dat, "hsuppDTPO",ngroups,horizon3)
  out3=plotroccurve(dat, "hsupp0hlgr",ngroups,horizon3)
  out4=plotroccurve(dat, "hsephlgr",ngroups,horizon3)
  out5=plotroccurve(dat, "hSKHThlgr",ngroups,horizon3)
  out6=plotroccurve(dat, "hsupphlgr",ngroups,horizon3)
  
  allx=cbind(out1[,1],out2[,1],out3[,1],out4[,1],out5[,1],out6[,1])
  ally=cbind(out1[,2],out2[,2],out3[,2],out4[,2],out5[,2],out6[,2])
  
  # out=data.frame(out)
  
  # names(out)=c("base","tu","schmid","frtrn","superpp","average")
  
  # ROC curves
  
  matplot(allx, ally, type="l", xlim=c(0, 100), ylim=c(0, 100), col = 1:6)
  legend(60, 30,
         legend = c("bench","suppDTPO","supp0hlgr","sephlgr","SKHThlgr","supphlgr"),
         col = 1:6,
         lty = 1,
         lwd = 1.5)
  abline(1, 1)
  
  
  # out
  
  
  
}
